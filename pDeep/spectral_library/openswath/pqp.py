import sqlite3
import numpy as np
import time

from ...utils.mass_calc import PeptideIonCalculator as ion_calc
from ..library_base import LibraryBase

__mod_dict = {
    "Carbamidomethyl[C]": "(UniMod:4)",
    "Oxidation[M]": "(UniMod:35)",
    "Phospho[S]": "(UniMod:21)",
    "Phospho[T]": "(UniMod:21)",
    "Phospho[Y]": "(UniMod:21)",
    "SILACnoLabel_13C(6)15N(2)[K]": "(8.014199)",
    "SILACnoLabel_13C(6)15N(4)[R]": "(10.008269)",
}

def pDeepFormat2PeptideModSeq(seq, modinfo):
    if not modinfo: return seq
    moditems = modinfo.strip(";").split(";")
    modlist = []
    for moditem in moditems:
        site, mod = moditem.split(",")
        modlist.append((int(site), mod))
    modlist.sort(reverse=True)
    for site, mod in modlist:
        if not mod in __mod_dict:
            print('[E] No {} in __mod_dict in openswath/pqp.py'.format(mod))
            return None
        seq = seq[:site] + __mod_dict[mod] + seq[site:]
    return seq

def PeptideModSeq2pDeepFormat(PeptideModSeq):
    site = PeptideModSeq.find('(')
    modlist = []
    while site != -1:
        if PeptideModSeq[site-1] == 'C': modlist.append('%d,%s'%(site, 'Carbamidomethyl[C]'))
        elif PeptideModSeq[site-1] == 'M': modlist.append('%d,%s'%(site, 'Oxidation[M]'))
        elif PeptideModSeq[site-1] == 'S': modlist.append('%d,%s'%(site, 'Phospho[S]'))
        elif PeptideModSeq[site-1] == 'T': modlist.append('%d,%s'%(site, 'Phospho[T]'))
        elif PeptideModSeq[site-1] == 'Y': modlist.append('%d,%s'%(site, 'Phospho[Y]'))
        else: return None, None
        PeptideModSeq = PeptideModSeq[:site] + PeptideModSeq[PeptideModSeq.find(')')+1:]
        site = PeptideModSeq.find('(', site)
    return PeptideModSeq, ";".join(modlist)

# not support insertion, only create from empty.pqp
class PQP(LibraryBase):
    def __init__(self):
        self.sql_conn = None
        self._ion_calc = ion_calc()
        self.peptide_dict = {}
        self.peptide_list = []
        self.precursor_peptide_dict = {}
        self.peptide_protein_dict = {}
        
    def Open(self, pqp_file):
        self.sql_conn = sqlite3.connect(pqp_file)
        self.cursor = self.sql_conn.cursor()
        self.pqp_file = pqp_file
        
    def Close(self):
        self.sql_conn.close()
        self.pqp_file = None
    
    def GetAllPeptides(self):
        start = time.perf_counter()
        self.peptide_dict = {}
        peptide_list = []
        cursor = self.cursor.execute("SELECT c.MODIFIED_SEQUENCE, a.CHARGE, a.LIBRARY_RT, b.precursor_id, b.peptide_id FROM precursor a INNER JOIN precursor_peptide_mapping b ON a.ID=b.precursor_id INNER JOIN peptide c on b.peptide_id=c.ID")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            charge = int(row[1])
            RT = float(row[2])*60
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = [row[0], charge, RT, row[3], row[4]] #items = [PeptideModSeq, CHARGE, RT, precursor_id, peptide_id]
        print("reading pqp time = %.3fs"%(time.perf_counter() - start))
        return peptide_list
        
    def _clear_transition(self):
        self.cursor.excute("DELETE FROM transition")
        self.cursor.excute("DELETE FROM transition_precursor_mapping")
        self.sql_conn.commit()
        
    # CREATE TABLE TRANSITION(ID INT PRIMARY KEY NOT NULL,TRAML_ID TEXT NULL,PRODUCT_MZ REAL NOT NULL,CHARGE INT NULL,TYPE CHAR(1) NULL,ANNOTATION TEXT NULL,ORDINAL INT NULL,DETECTING INT NOT NULL,IDENTIFYING INT NOT NULL,QUANTIFYING INT NOT NULL,LIBRARY_INTENSITY REAL NULL,DECOY INT NOT NULL)
    def _insert_transition(self):
        
        
    # peak_selection = "topK"? or "intensity"
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}, peak_selection = "intensity", threshold = 0.05, mass_upper = 2000):
        print("updating dlib ...")
        count = 0
        start = time.perf_counter()
        self._clear_transition()
        
        self.transition_update_list = []
        self.transition_count = 0
        insert_list = []
        insert_pep2pro_list = []
        update_sql = "UPDATE entries SET MassEncodedLength = ?, MassArray = ?, IntensityEncodedLength = ?, IntensityArray = ?, RTInSeconds = ? WHERE PeptideModSeq = ? AND PrecursorCharge = ?"
        
        insert_sql = "INSERT INTO entries(PeptideModSeq, PeptideSeq, PrecursorMz, PrecursorCharge, RTInSeconds, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, SourceFile, Copies, Score) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, 'pDeep', 1, 0.001)"
        
        insert_pep2pro_sql = "INSERT INTO peptidetoprotein(PeptideSeq, isDecoy, ProteinAccession) VALUES (?, '0', ?)"
        
        def _encode(masses, intens):
            return len(masses)*8, EncodeMassList(masses), len(intens)*4, EncodeIntensityList(intens)
            
        def _check_protein(PeptideSeq, ProteinACs, insert_pep2pro_list):
            if not ProteinACs: return insert_pep2pro_list
            cursor = self.cursor.execute("SELECT isDecoy FROM peptidetoprotein WHERE PeptideSeq = '%s'"%PeptideSeq)
            if not cursor.fetchone():
                for ProteinAccession in ProteinACs.split(";"):
                    insert_pep2pro_list.append((PeptideSeq, ProteinAccession))
            return insert_pep2pro_list
        
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            count += 1
            seq, mod, charge = pepinfo.split("|")
            charge = int(charge)
            max_ion_charge = 2
            masses, pepmass = self._ion_calc.calc_by_and_pepmass(seq, mod, max_ion_charge)
            # print(seq, mod, masses)
            b_sites = np.tile(np.arange(1, len(seq)).reshape(-1,1), [1,max_ion_charge])
            b_types = np.array(['b']*((len(seq)-1)*max_ion_charge)).reshape(-1, max_ion_charge)
            y_sites = len(seq)-b_sites
            y_types = np.array(['y']*((len(seq)-1)*max_ion_charge)).reshape(-1, max_ion_charge)
            sites = np.concatenate((b_sites, y_sites), axis=1)
            types = np.concatenate((b_types, y_types), axis=1)
            intens = intensities[:,:masses.shape[1]]
            intens[0, 0:2] = 0 #b1+/b1++ = 0
            intens[-1, 2:] = 0 #y1+/y1++ = 0, do not consider y1/b1 in the library
            masses = masses.reshape(-1)
            intens = intens.reshape(-1)
            sites = sites.reshape(-1)
            types = types.reshape(-1)
            intens = intens/np.max(intens)
            
            if peak_selection == "intensity": 
                masses = masses[intens > threshold]
                sites = sites[intens > threshold]
                types = types[intens > threshold]
                intens = intens[intens > threshold]*10000
            else:
                intens = intens*10000
            intens = intens[masses < mass_upper]
            sites = sites[masses < mass_upper]
            types = types[masses < mass_upper]
            masses = masses[masses < mass_upper]
            RT = _prediction.GetRetentionTime(pepinfo)
            
            if pepinfo in self.peptide_dict:
                item = self.peptide_dict[pepinfo]
                update_list.append((*_encode(masses, intens), item[2], item[0], item[1])) #item[2] = RT in dlib
            else:
                pepmass = pepmass / charge + self._ion_calc.base_mass.mass_proton
                insert_pep2pro_list = _check_protein(seq, peptide_to_protein_dict[seq] if seq in peptide_to_protein_dict else "", insert_pep2pro_list)
                insert_list.append((pDeepFormat2PeptideModSeq(seq, mod), seq, pepmass, charge, RT if RT is not None else 0, *_encode(masses, intens)))
            if count%10000 == 0:
                print("[SQL UPDATE] {:.1f}%".format(100.0*count/len(_prediction.peptide_intensity_dict)), end="\r")
                if count%1000000 == 0:
                    self.cursor.executemany(update_sql, update_list)
                    self.cursor.executemany(insert_sql, insert_list)
                    self.cursor.executemany(insert_pep2pro_sql, insert_pep2pro_list)
                    update_list = []
                    insert_list = []
                    insert_pep2pro_list = []
        self.cursor.executemany(update_sql, update_list)
        self.cursor.executemany(insert_sql, insert_list)
        self.cursor.executemany(insert_pep2pro_sql, insert_pep2pro_list)
        print("[SQL UPDATE] 100%: {}".format(self.dlib_file))
        self.sql_conn.commit()
        print("updating dlib time = %.3fs"%(time.perf_counter()-start))
        
        
        