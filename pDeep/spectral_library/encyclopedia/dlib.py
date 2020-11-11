import sqlite3
import zlib
import struct
import numpy as np
import time

from ...utils.mass_calc import PeptideIonCalculator as ion_calc
from ..library_base import LibraryBase

__mod_dict = {
    "Carbamidomethyl[C]": "[+57.021464]",
    "Oxidation[M]": "[+15.994915]",
    "Phospho[S]": "[+79.966331]",
    "Phospho[T]": "[+79.966331]",
    "Phospho[Y]": "[+79.966331]",
    "Acetyl[ProteinN-term]": "[+42.010565]",
    "Acetyl[AnyN-term]": "[+42.010565]",
    "Label_13C(6)15N(2)[K]": "[+8.014199]",
    "Label_13C(6)15N(4)[R]": "[+10.008269]",
}
# from ...prediction import pDeepPrediction as prediction

# pack/unpack(fmt,...), fmt=">" means big-endian

class DLIB(LibraryBase):
    def __init__(self, pDeepParam=None):
        super(self.__class__, self).__init__(pDeepParam)
        self.sql_conn = None
        self.peptide_dict = {}
        self.peptide_list = []
        self.lib_version = "0.1.14"
        
    def Open(self, dlib_file):
        self.sql_conn = sqlite3.connect(dlib_file)
        self.cursor = self.sql_conn.cursor()
        self.dlib_file = dlib_file
        
    def Close(self):
        self.sql_conn.close()
        self.dlib_file = None
            
    def CreateTables(self):
        self.cursor.execute("CREATE TABLE peptidetoprotein( PeptideSeq string not null, IsDecoy int not null, ProteinAccession string not null)") #, PRIMARY KEY (PeptideSeq) )")
        self.cursor.execute("CREATE TABLE entries ( PrecursorMz double not null, PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, Copies int not null, RTInSeconds double not null, Score double not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, CorrelationEncodedLength int, CorrelationArray blob, RTInSecondsStart double, RTInSecondsStop double, MedianChromatogramEncodedLength int, MedianChromatogramArray blob, SourceFile string not null )")
        self.cursor.execute("CREATE TABLE metadata ( Key string not null, Value string not null, PRIMARY KEY (Key) )")
        self.cursor.execute("CREATE TABLE fragmentquants ( PrecursorCharge int not null, PeptideModSeq string not null, SourceFile string not null, IonType string not null, FragmentMass double not null, Correlation double not null, DeltaMassPPM double not null, Intensity double not null, FOREIGN KEY (PrecursorCharge, PeptideModSeq, SourceFile) REFERENCES entries (PrecursorCharge, PeptideModSeq, SourceFile) )")
        self.cursor.execute("CREATE TABLE peptidequants ( PrecursorCharge int not null, PeptideModSeq string not null, SourceFile string not null, RTInSecondsStart double not null, RTInSecondsStop double not null, TotalIntensity double not null, NumberOfQuantIons int not null, BestFragmentCorrelation double not null, BestFragmentDeltaMassPPM double not null, MedianChromatogramEncodedLength int not null, MedianChromatogramArray blob not null,PRIMARY KEY (PrecursorCharge, PeptideModSeq, SourceFile), FOREIGN KEY (PrecursorCharge, PeptideModSeq, SourceFile) REFERENCES entries (PrecursorCharge, PeptideModSeq, SourceFile) )")
        
        self.cursor.execute("INSERT INTO metadata(key,value) VALUES('version', '%s')"%self.lib_version)
        
    def GetAllPeptides(self):
        start = time.perf_counter()
        self.peptide_dict = {}
        peptide_list = []
        cursor = self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='peptidetoprotein'")
        if cursor.fetchone():
            cursor = self.cursor.execute("SELECT entries.PeptideModSeq, entries.PrecursorCharge, entries.RTInSeconds, peptidetoprotein.ProteinAccession FROM entries INNER JOIN peptidetoprotein ON entries.PeptideSeq=peptidetoprotein.PeptideSeq")
        else:
            cursor = self.cursor.execute("SELECT entries.PeptideModSeq, entries.PrecursorCharge, entries.RTInSeconds, proteins.ProteinAccessions FROM entries INNER JOIN proteins ON entries.PeptideSeq=proteins.PeptideSeq")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            if not seq:
                print("[pDeep Warning] unknown modseq '%s'"%row[0])
                continue
            charge = int(row[1])
            RT = float(row[2])
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = (row[0], charge, RT, '', -1, row[3]) #items = [PeptideModSeq, PrecursorCharge, RT, raw, scan, protein]
        print("reading dlib time = %.3fs"%(time.perf_counter() - start))
        return peptide_list
        
    # peak_selection = "topK" or "intensity"
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}):
        self.decoy = None
        print("[pDeep Info] updating dlib ...")
        count = 0
        start = time.perf_counter()
        
        update_list = []
        insert_list = []
        insert_pep2pro_list = set()
        update_sql = "UPDATE entries SET MassEncodedLength = ?, MassArray = ?, IntensityEncodedLength = ?, IntensityArray = ?, RTInSeconds = ? WHERE PeptideModSeq = ? AND PrecursorCharge = ?"
        
        insert_sql = "INSERT INTO entries(PeptideModSeq, PeptideSeq, PrecursorMz, PrecursorCharge, RTInSeconds, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, SourceFile, Copies, Score) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, 'pDeep', 1, 0.001)"
        
        insert_pep2pro_sql = "INSERT INTO peptidetoprotein(PeptideSeq, isDecoy, ProteinAccession) VALUES (?, '0', ?)"
        
        def _encode(masses, intens):
            return len(masses)*8, EncodeMassList(masses), len(intens)*4, EncodeIntensityList(intens)
            
        def _check_protein(PeptideSeq, ProteinACs, insert_pep2pro_list):
            if not ProteinACs: return insert_pep2pro_list
            cursor = self.cursor.execute("SELECT isDecoy FROM peptidetoprotein WHERE PeptideSeq = '%s'"%PeptideSeq)
            if not cursor.fetchone():
                for ProteinAccession in ProteinACs.split("/"):
                    insert_pep2pro_list.add((PeptideSeq, ProteinAccession))
            return insert_pep2pro_list
        
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            seq, mod, charge = pepinfo.split("|")
            charge = int(charge)
            
            pepmass, masses, intens, sites, types, charges, decoy_seq, decoy_mod, decoy_masses = self._calc_ions(seq, mod, charge, intensities)
                
            RT = _prediction.GetRetentionTime(pepinfo)
            
            # if pepinfo in self.peptide_dict:
                # item = self.peptide_dict[pepinfo]
                # update_list.append((*_encode(masses, intens), item[2], item[0], item[1])) #item[2] = RT in dlib
            
            insert_pep2pro_list = _check_protein(seq, peptide_to_protein_dict[seq] if seq in peptide_to_protein_dict else "", insert_pep2pro_list)
            modseq = pDeepFormat2PeptideModSeq(seq, mod)
            if modseq: insert_list.append((modseq, seq, pepmass, charge, RT if RT is not None else 0, *_encode(masses, intens)))
            count += 1
            if count%10000 == 0:
                print("[SQL UPDATE] {:.1f}%".format(100.0*count/len(_prediction.peptide_intensity_dict)), end="\r")
                if count%1000000 == 0:
                    # self.cursor.executemany(update_sql, update_list)
                    self.cursor.executemany(insert_sql, insert_list)
                    self.cursor.executemany(insert_pep2pro_sql, insert_pep2pro_list)
                    # update_list = []
                    insert_list = []
                    insert_pep2pro_list = set()
        # self.cursor.executemany(update_sql, update_list)
        self.cursor.executemany(insert_sql, insert_list)
        self.cursor.executemany(insert_pep2pro_sql, insert_pep2pro_list)
        print("[SQL UPDATE] 100%: {}".format(self.dlib_file))
        self.sql_conn.commit()
        print("[pDeep Info] updating dlib time = %.3fs"%(time.perf_counter()-start))

def pDeepFormat2PeptideModSeq(seq, modinfo):
    if not modinfo: return seq
    moditems = modinfo.strip(";").split(";")
    modlist = []
    for moditem in moditems:
        site, mod = moditem.split(",")
        modlist.append((int(site), mod))
    modlist.sort(reverse=True)
    for site, mod in modlist:
        if not mod in __mod_dict: return None
        seq = seq[:site] + __mod_dict[mod] + seq[site:]
    return seq

def PeptideModSeq2pDeepFormat(PeptideModSeq):
    site = PeptideModSeq.find('[')
    modlist = []
    while site != -1:
        right = PeptideModSeq.find(']')
        int_mass = round(float(PeptideModSeq[site+1:right]))
        if PeptideModSeq[site-1] == 'C' and int_mass == 57: modlist.append('%d,%s'%(site, 'Carbamidomethyl[C]'))
        elif PeptideModSeq[site-1] == 'M' and int_mass == 16: modlist.append('%d,%s'%(site, 'Oxidation[M]'))
        elif PeptideModSeq[site-1] == 'S' and int_mass == 80: modlist.append('%d,%s'%(site, 'Phospho[S]'))
        elif PeptideModSeq[site-1] == 'T' and int_mass == 80: modlist.append('%d,%s'%(site, 'Phospho[T]'))
        elif PeptideModSeq[site-1] == 'Y' and int_mass == 80: modlist.append('%d,%s'%(site, 'Phospho[Y]'))
        elif site == 1 and int_mass == 42: modlist.append('0,%s'%('Acetyl[ProteinN-term]'))
        else: return None, None
        PeptideModSeq = PeptideModSeq[:site] + PeptideModSeq[right+1:]
        site = PeptideModSeq.find('[', site)
    return PeptideModSeq, ";".join(modlist)

def DecodeDoubleList(code):
    bytes = zlib.decompress(code)
    return struct.unpack(">%dd"%(len(bytes)/8), bytes)
    
def DecodeFloatList(code):
    bytes = zlib.decompress(code)
    return struct.unpack(">%df"%(len(bytes)/4), bytes)
    
def DecodeIntensityList(lst):
    return DecodeFloatList(lst)
    
def DecodeMassList(lst):
    return DecodeDoubleList(lst)
    
def EncodeDoubleList(lst):
    bytes = struct.pack(">%dd"%len(lst), *lst)
    return zlib.compress(bytes)
    
def EncodeFloatList(lst):
    bytes = struct.pack(">%df"%len(lst), *lst)
    return zlib.compress(bytes)
    
def EncodeIntensityList(lst):
    return EncodeFloatList(np.array(lst, dtype=np.float32))
    
def EncodeMassList(lst):
    return EncodeDoubleList(np.array(lst, dtype=np.float64))