import sqlite3
import numpy as np
import time

from ...utils.mass_calc import PeptideIonCalculator as ion_calc
from ..library_base import LibraryBase
from .tsv import pDeepFormat2PeptideModSeq, PeptideModSeq2pDeepFormat

class OSW(LibraryBase):
    def __init__(self):
        self.sql_conn = None
        
    def Open(self, osw_file):
        self.sql_conn = sqlite3.connect(osw_file)
        self.cursor = self.sql_conn.cursor()
        self.osw_file = osw_file
        
    def Close(self):
        self.sql_conn.close()
        self.osw_file = None
    
    def GetAllPeptides(self):
        start = time.perf_counter()
        self.peptide_dict = {}
        peptide_list = []
        cursor = self.cursor.execute("SELECT pep.MODIFIED_SEQUENCE, pre.CHARGE, fea.EXP_RT FROM score_ms2 ms2 INNER JOIN feature fea ON ms2.FEATURE_ID=fea.ID INNER JOIN precursor pre ON pre.ID=fea.precursor_id INNER JOIN precursor_peptide_mapping prepep ON prepep.precursor_id=fea.precursor_id INNER JOIN peptide pep ON pep.ID=prepep.peptide_id where ms2.rank=1 and ms2.qvalue<0.01")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            charge = int(row[1])
            RT = float(row[2])
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = [row[0], charge, RT] #items = [PeptideModSeq, CHARGE, RT]
        print("[pDeep Info] reading osw time = %.3fs"%(time.perf_counter() - start))
        return peptide_list

# not support insertion, only create from empty.pqp
class PQP(LibraryBase):
    def __init__(self):
        self.sql_conn = None
        self._ion_calc = ion_calc()
        self.peptide_dict = {}
        self.peptide_list = []
        self.precursor_peptide_dict = {}
        self.peptide_protein_dict = {}
        
        self.decoy = "pseudo_reverse" #or reverse
        self.decoy_tag = "DECOY_"
        
        #CREATE TABLE TRANSITION(ID INT PRIMARY KEY NOT NULL,TRAML_ID TEXT NULL,PRODUCT_MZ REAL NOT NULL,CHARGE INT NULL,TYPE CHAR(1) NULL,ANNOTATION TEXT NULL,ORDINAL INT NULL,DETECTING INT NOT NULL,IDENTIFYING INT NOT NULL,QUANTIFYING INT NOT NULL,LIBRARY_INTENSITY REAL NULL,DECOY INT NOT NULL)
        self.transition_insert_sql = "INSERT INTO TRANSITION(ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ANNOTATION, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES(?, ?, ?, ?, ?, ?, ?, 1, 0, 1, ?, ?)"
        
        # CREATE TABLE PRECURSOR(ID INT PRIMARY KEY NOT NULL,TRAML_ID TEXT NULL,GROUP_LABEL TEXT NULL,PRECURSOR_MZ REAL NOT NULL,CHARGE INT NULL,LIBRARY_INTENSITY REAL NULL,LIBRARY_RT REAL NULL,LIBRARY_DRIFT_TIME REAL NULL,DECOY INT NOT NULL)
        self.precursor_insert_sql = "INSERT INTO PRECURSOR(ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, LIBRARY_RT, LIBRARY_DRIFT_TIME, DECOY) VALUES(?, ?, ?, ?, ?, '', ?, -1.0, ?)"
        
        # CREATE TABLE PEPTIDE(ID INT PRIMARY KEY NOT NULL,UNMODIFIED_SEQUENCE TEXT NOT NULL,MODIFIED_SEQUENCE TEXT NOT NULL,DECOY INT NOT NULL)
        self.peptide_insert_sql = "INSERT INTO PEPTIDE(ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES(?, ?, ?, ?)"
        
        # CREATE TABLE PROTEIN(ID INT PRIMARY KEY NOT NULL,PROTEIN_ACCESSION TEXT NOT NULL,DECOY INT NOT NULL)
        self.protein_insert_sql = "INSERT INTO protein(ID, PROTEIN_ACCESSION, DECOY) VALUES(?, ?, ?)"
        
        # CREATE TABLE PEPTIDE_GENE_MAPPING(PEPTIDE_ID INT NOT NULL,GENE_ID INT NOT NULL)
        self.peptide2gene_sql = "INSERT INTO PEPTIDE_GENE_MAPPING(PEPTIDE_ID, GENE_ID) VALUES(?, 0)"
        
        # CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID INT NOT NULL,PROTEIN_ID INT NOT NULL)
        self.peptide2protein_sql = "INSERT INTO PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID, PROTEIN_ID) VALUES(?, ?)"
        
        # CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INT NOT NULL,PEPTIDE_ID INT NOT NULL)
        self.precursor2peptide_sql = "INSERT INTO PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID, PEPTIDE_ID) VALUES(?, ?)"
        
        # CREATE TABLE TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID INT NOT NULL,PRECURSOR_ID INT NOT NULL)
        self.transition2precursor_sql = "INSERT INTO TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID, PRECURSOR_ID) VALUES(?, ?)"
        
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
        cursor = self.cursor.execute("SELECT pep.MODIFIED_SEQUENCE, pre.CHARGE, pre.LIBRARY_RT FROM precursor pre INNER JOIN precursor_peptide_mapping prepep ON pre.ID=prepep.precursor_id INNER JOIN peptide pep on prepep.peptide_id=pep.ID")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            charge = int(row[1])
            RT = float(row[2])*60
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = [row[0], charge, RT, ""] #items = [PeptideModSeq, CHARGE, RT, protein]
        print("[pDeep Info] reading pqp time = %.3fs"%(time.perf_counter() - start))
        return peptide_list
        
    # def _write_one_decoy_peptide(self, seq, mod, pre_charge, pepmass, masses, intens, charges, types, sites, RT, protein, transition_count, transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list):
        # if self.decoy == 'reverse':
            # seq = seq[::-1]
            # if mod:
                # mods = mod.strip(";").split(";")
                # modlist = []
                # for onemod in mods:
                    # site, modname = onemod.split(",")
                    # site = int(site)
                    # if site <= len(seq) and site != 0:
                        # site = len(seq)-site+1
                    # modlist.append((site, modname))
                # modlist.sort(key=lambda x: x[0])
                # mod = ";".join(['%d,%s'%(site, modname) for site, modname in modlist])
            # labeled_seq = pDeepFormat2PeptideModSeq(seq, mod)
            
            # changed_y_idx = (types=='y')
            # changed_b_idx = (types=='b')
            # changed_idx = np.logical_or(changed_y_idx, changed_b_idx)
            # sites[changed_idx] = len(seq)-sites[changed_idx]
            # peptide_mass = (pepmass - self._ion_calc.base_mass.mass_proton)*pre_charge
            # peptide_charged_masses = peptide_mass*charges + self._ion_calc.base_mass.mass_proton
            # masses[changed_y_idx] = peptide_charged_masses[changed_y_idx] - masses[changed_y_idx] + self._ion_calc.get_aamass(seq[-1])/charges[changed_y_idx]
            # masses[changed_b_idx] = peptide_charged_masses[changed_b_idx] - masses[changed_b_idx] - self._ion_calc.get_aamass(seq[-1])/charges[changed_b_idx]
            # types[changed_y_idx] = 'b'
            # types[changed_b_idx] = 'y'
        # else:
            # seq = seq[:-1][::-1]+seq[-1]
            # if mod:
                # mods = mod.strip(";").split(";")
                # modlist = []
                # for onemod in mods:
                    # site, modname = onemod.split(",")
                    # site = int(site)
                    # if site < len(seq) and site != 0:
                        # site = len(seq)-site
                    # modlist.append((site, modname))
                # modlist.sort(key=lambda x: x[0])
                # mod = ";".join(['%d,%s'%(site, modname) for site, modname in modlist])
            # labeled_seq = pDeepFormat2PeptideModSeq(seq, mod)
            
            # changed_y_idx = np.logical_and(sites > 1, types=='y')
            # changed_b_idx = np.logical_and(sites < len(seq)-1, types=='b')
            # changed_idx = np.logical_or(changed_y_idx, changed_b_idx)
            # sites[changed_idx] = (len(seq)-1)-sites[changed_idx]
            # peptide_mass = (pepmass - self._ion_calc.base_mass.mass_proton)*pre_charge
            # peptide_charged_masses = peptide_mass*charges + self._ion_calc.base_mass.mass_proton
            # masses[changed_y_idx] = peptide_charged_masses[changed_y_idx] - masses[changed_y_idx] + self._ion_calc.get_aamass(seq[-1])/charges[changed_y_idx]
            # masses[changed_b_idx] = peptide_charged_masses[changed_b_idx] - masses[changed_b_idx] - self._ion_calc.get_aamass(seq[-1])/charges[changed_b_idx]
            # types[changed_y_idx] = 'b'
            # types[changed_b_idx] = 'y'
            
        # return self._write_one_peptide(seq, mod, pre_charge, pepmass, masses, intens, charges, types, sites, RT, protein, transition_count, transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list, 1)
        
    def _write_one_peptide(self, seq, mod, pre_charge, pepmass, masses, intens, charges, types, sites, RT, protein, transition_count, transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list, decoy=0):
    
        labeled_seq = pDeepFormat2PeptideModSeq(seq, mod)
        
        if decoy: protein = self.decoy_tag+protein
        
        precursor_row = []
        # "INSERT INTO PRECURSOR(ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, LIBRARY_RT, LIBRARY_DRIFT_TIME, DECOY) VALUES(?, ?, ?, ?, ?, NULL, ?, -1.0, ?)"
        precursor_id = len(precursor2peptide)
        label = "%s%d_%s_%d"%("" if not decoy else self.decoy_tag, precursor_id, labeled_seq, pre_charge)
        precursor_row.append(precursor_id)
        precursor_row.append(precursor_id)
        precursor_row.append(label)
        precursor_row.append(pepmass)
        precursor_row.append(pre_charge)
        precursor_row.append(RT/60)
        precursor_row.append(decoy)
        precursor_list.append(precursor_row)
        
        if labeled_seq in peptide_dict:
            precursor2peptide[precursor_id] = peptide_dict[labeled_seq]
            precursor2peptide_list.append((precursor_id, peptide_dict[labeled_seq]))
        else:
            precursor2peptide[precursor_id] = len(peptide_dict)
            precursor2peptide_list.append((precursor_id, len(peptide_dict)))
            peptide_list.append((len(peptide_dict), seq, labeled_seq, decoy))
            
            peptide_dict[labeled_seq] = len(peptide_dict)
            
        if protein in protein_dict:
            peptide2protein[peptide_dict[labeled_seq]] = protein_dict[protein]
            peptide2protein_list.append((peptide_dict[labeled_seq],protein_dict[protein]))
        else:
            peptide2protein[peptide_dict[labeled_seq]] = len(protein_dict)
            peptide2protein_list.append((peptide_dict[labeled_seq],len(protein_dict)))
            protein_list.append((len(protein_dict), protein, decoy))
            
            protein_dict[protein] = len(protein_dict)
            
        for mz, inten, charge, ion_type, site in zip(masses, intens, charges, types, sites):
            # "INSERT INTO TRANSITION(ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ANNOTATION, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES(?, ?, ?, ?, ?, ?, ?, 1, 0, 1, ?, ?)"
            transition_row = []
            transition_row.append(transition_count)
            transition_row.append("%s%d_%s_%d"%("" if not decoy else self.decoy_tag, transition_count, labeled_seq, pre_charge))
            transition_row.append(mz)
            transition_row.append(int(charge))
            transition_row.append(ion_type)
            transition_row.append('%s%d^%d/0'%(ion_type, site, charge) if charge > 1 else '%s%d/0'%(ion_type, site))
            transition_row.append(int(site))
            transition_row.append(float(inten))
            transition_row.append(decoy)
            transition_list.append(transition_row)
            transition2precursor_list.append((transition_count, precursor_id))
            transition_count += 1
            
        return transition_count
        
    def _get_decoy_peptide(self, seq, mod):
        if self.decoy == 'reverse':
            seq = seq[::-1]
            if mod:
                mods = mod.strip(";").split(";")
                modlist = []
                for onemod in mods:
                    site, modname = onemod.split(",")
                    site = int(site)
                    if site <= len(seq) and site != 0:
                        site = len(seq)-site+1
                    modlist.append((site, modname))
                modlist.sort(key=lambda x: x[0])
                mod = ";".join(['%d,%s'%(site, modname) for site, modname in modlist])
        else:
            seq = seq[:-1][::-1]+seq[-1]
            if mod:
                mods = mod.strip(";").split(";")
                modlist = []
                for onemod in mods:
                    site, modname = onemod.split(",")
                    site = int(site)
                    if site < len(seq) and site != 0:
                        site = len(seq)-site
                    modlist.append((site, modname))
                modlist.sort(key=lambda x: x[0])
                mod = ";".join(['%d,%s'%(site, modname) for site, modname in modlist])
        return seq, mod
        
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}, min_intensity = 0.1, least_n_peaks = 6, max_mz = 2000):
        print("[pDeep Info] updating pqp ...")
        transition_count = 0
        count = 0
        transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list = [], [], [], {}, [], {}, [], {}, [], {}, []
        start = time.perf_counter()
        
        if self.decoy: total = len(_prediction.peptide_intensity_dict)*2
        else: total = len(_prediction.peptide_intensity_dict)
        
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            seq, mod, charge = pepinfo.split("|")
            
            pre_charge = int(charge)
            max_ion_charge = 2
            masses, pepmass = self._ion_calc.calc_by_and_pepmass(seq, mod, max_ion_charge)
            
            if self.decoy: 
                decoy_seq, decoy_mod = self._get_decoy_peptide(seq, mod)
                decoy_masses, _ = self._ion_calc.calc_by_and_pepmass(decoy_seq, decoy_mod, max_ion_charge)
            pepmass = pepmass / pre_charge + self._ion_calc.base_mass.mass_proton
            # print(seq, mod, masses)
            
            b_sites = np.tile(np.arange(1, len(seq)).reshape(-1,1), [1,max_ion_charge])
            y_sites = len(seq)-b_sites
            sites = np.concatenate((b_sites, y_sites), axis=1)
            
            b_types = np.array(['b']*((len(seq)-1)*max_ion_charge)).reshape(-1, max_ion_charge)
            y_types = np.array(['y']*((len(seq)-1)*max_ion_charge)).reshape(-1, max_ion_charge)
            types = np.concatenate((b_types, y_types), axis=1)
            
            charges = np.concatenate([np.full(len(seq)-1, i, dtype=int).reshape(-1,1) for i in range(1, max_ion_charge+1)], axis=1)
            charges = np.concatenate((charges, charges), axis=1)
            
            intens = intensities[:,:masses.shape[1]]
            intens[0, 0:2] = 0 #b1+/b1++ = 0
            intens[-1, 2:] = 0 #y1+/y1++ = 0, do not consider y1/b1 in the library
            masses = masses.reshape(-1)
            intens = intens.reshape(-1)
            sites = sites.reshape(-1)
            types = types.reshape(-1)
            charges = charges.reshape(-1)
            if self.decoy: decoy_masses = decoy_masses.reshape(-1)
            
            intens[np.abs(masses - pepmass) < 10] = 0 #delete ions around precursor m/z
            intens = intens/np.max(intens)
                
            intens = intens[masses <= max_mz]
            sites = sites[masses <= max_mz]
            types = types[masses <= max_mz]
            charges = charges[masses <= max_mz]
            if self.decoy: decoy_masses = decoy_masses[masses <= max_mz]
            masses = masses[masses <= max_mz]
            
            
            if len(masses[intens > min_intensity]) >= least_n_peaks:
                masses = masses[intens > min_intensity]
                sites = sites[intens > min_intensity]
                types = types[intens > min_intensity]
                charges = charges[intens > min_intensity]
                if self.decoy: decoy_masses = decoy_masses[intens > min_intensity]
                intens = intens[intens > min_intensity] * 10000
            else:
                indices = np.argsort(intens)[::-1]
                masses = masses[indices[:least_n_peaks]]
                sites = sites[indices[:least_n_peaks]]
                types = types[indices[:least_n_peaks]]
                charges = charges[indices[:least_n_peaks]]
                if self.decoy: decoy_masses = decoy_masses[indices[:least_n_peaks]]
                intens = intens[indices[:least_n_peaks]]*10000
            
            
            RT = _prediction.GetRetentionTime(pepinfo)
            
            if seq in peptide_to_protein_dict:
                protein = peptide_to_protein_dict[seq]
                protein = str(protein.count("/")+1)+"/"+protein
            elif pepinfo in self.peptide_dict:
                protein = self.peptide_dict[pepinfo][-1]
            else:
                protein = "1/pDeep"
                
            transition_count = self._write_one_peptide(seq, mod, pre_charge, pepmass, masses, intens, charges, types, sites, RT, protein, transition_count, transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list)
            count += 1
            if self.decoy:
                transition_count = self._write_one_peptide(decoy_seq, decoy_mod, pre_charge, pepmass, decoy_masses, intens, charges, types, sites, RT, protein, transition_count, transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list, 1)
                count += 1
            if count%10000 == 0:
                print("[PQP UPDATE] {:.1f}%".format(100.0*count/total), end="\r")
                if count%100000 == 0:
                    self.cursor.executemany(self.transition_insert_sql, transition_list)
                    self.cursor.executemany(self.precursor_insert_sql, precursor_list)
                    self.cursor.executemany(self.peptide_insert_sql, peptide_list)
                    self.cursor.executemany(self.protein_insert_sql, protein_list)
                    self.cursor.executemany(self.transition2precursor_sql, transition2precursor_list)
                    self.cursor.executemany(self.precursor2peptide_sql, precursor2peptide_list)
                    self.cursor.executemany(self.peptide2protein_sql, peptide2protein_list)
                    self.sql_conn.commit()
                    transition_list, transition2precursor_list, precursor_list, precursor2peptide_list, peptide_list, peptide2protein_list, protein_list = [], [], [], [], [], [], []
                    
        self.cursor.executemany(self.transition_insert_sql, transition_list)
        self.cursor.executemany(self.precursor_insert_sql, precursor_list)
        self.cursor.executemany(self.peptide_insert_sql, peptide_list)
        self.cursor.executemany(self.protein_insert_sql, protein_list)
        self.cursor.executemany(self.transition2precursor_sql, transition2precursor_list)
        self.cursor.executemany(self.precursor2peptide_sql, precursor2peptide_list)
        self.cursor.executemany(self.peptide2protein_sql, peptide2protein_list)
        self.cursor.executemany(self.peptide2gene_sql, [(peptide_id, ) for peptide_id in peptide_dict.values()])
        self.sql_conn.commit()
        print("[PQP UPDATE] 100%: {}".format(self.pqp_file))
        print("[pDeep Info] updating pqp time = %.3fs"%(time.perf_counter()-start))
        
        
        