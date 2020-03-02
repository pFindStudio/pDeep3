import sqlite3
import zlib
import struct
import numpy as np

from ...utils.mass_calc import PeptideIonCalculator as ion_calc

mod_dict = {
    "Carbamidomethyl[C]": "[+57.021464]",
    "Oxidation[M]": "[+15.994915]",
    "Phospho[S]": "[+79.966331]",
    "Phospho[T]": "[79.966331]",
    "Phospho[Y]": "[79.966331]",
}
# from ...prediction import pDeepPrediction as prediction

# pack/unpack(fmt,...), fmt=">" means big-endian

class DLIB(object):
    def __init__(self):
        self.sql_conn = None
        self._ion_calc = ion_calc()
        
    def Open(self, dlib_file):
        self.sql_conn = sqlite3.connect(dlib_file)
        self.dlib_file = dlib_file
        
    def Close(self):
        self.sql_conn.close()
        self.dlib_file = None
        
    def GetAllPeptides(self):
        self.peptide_dict = {}
        peptide_list = []
        cursor = self.sql_conn.execute("SELECT PeptideModSeq, PrecursorCharge, RTInSeconds FROM entries")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            charge = int(row[1])
            RT = float(row[2])
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = [row[0], charge, RT, None, None] #items = [PeptideModSeq, PrecursorCharge, PredictedMassList, PredictedIntensityList]
        return peptide_list
        
    def UpdateMassIntensity(self, PeptideModSeq, PrecursorCharge, MassList, IntensityList, RT, commit_now = False):
        sql = "UPDATE entries SET MassEncodedLength = ?, MassArray = ?, IntensityEncodedLength = ?, IntensityArray = ?, RTInSeconds = ? WHERE PeptideModSeq = '%s' AND PrecursorCharge = %d"%(PeptideModSeq, PrecursorCharge)
        self.sql_conn.execute(sql, (len(MassList)*8, EncodeMassList(MassList), len(IntensityList)*4, EncodeIntensityList(IntensityList), RT))
        if commit_now: sql_conn.commit()
    
    def InsertNewPeptide(self, PeptideModSeq, PeptideSeq, PrecursorMz, PrecursorCharge, MassList, IntensityList, RT, ProteinACs = "", commit_now = False):
        sql = "INSERT INTO entries(PeptideModSeq, PeptideSeq, PrecursorMz, PrecursorCharge, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, RTInSeconds, SourceFile, Copies, Score) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, 'pDeep', 1, 0.001)"
        self.sql_conn.execute(sql, (
                PeptideModSeq, 
                PeptideSeq, 
                PrecursorMz, 
                PrecursorCharge, 
                len(MassList)*8, 
                EncodeMassList(MassList), 
                len(IntensityList)*4, 
                EncodeIntensityList(IntensityList), 
                RT,
            ))
        if ProteinACs:
            for ProteinAccession in ProteinACs.split(";"):
                cursor = self.sql_conn.execute("SELECT isDecoy FROM peptidetoprotein WHERE PeptideSeq = '%s'"%PeptideSeq)
                if not cursor.fetchone():
                    self.sql_conn.execute("INSERT INTO peptidetoprotein(PeptideSeq, isDecoy, ProteinAccession) VALUES ('%s', '0', '%s')"%(PeptideSeq, ProteinAccession))
        if commit_now: sql_conn.commit()
    
    def UpdateByPeptideDict(self, peptide_dict):
        for pepinfo, (modseq, charge, mass, intensity) in peptide_dict.items():
            if not mass: continue
            self.UpdateMassIntensity(modseq, charge, mass, intensity)
        self.sql_conn.commit()
    
    # peak_selection = "topK" or "intensity"
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}, peak_selection = "intensity", threshold = 0.01, mass_upper = 1500):
        count = 0
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            count += 1
            seq, mod, charge = pepinfo.split("|")
            masses, pepmass = self._ion_calc.calc_by_and_pepmass(seq, mod, 2)
            # print(seq, mod, masses)
            intens = intensities[:,:masses.shape[1]]
            intens[0, 0:2] = 0 #b1+/b1++ = 0
            intens[-1, 2:] = 0 #y1+/y1++ = 0, do not consider y1/b1 in the library
            masses = masses.reshape(-1)
            intens = intens.reshape(-1)
            intens = intens/np.max(intens)
            
            if peak_selection == "intensity": 
                masses = masses[intens > threshold]
                intens = intens[intens > threshold]*10000
            else:
                intens = intens*10000
            intens = intens[masses < mass_upper]
            masses = masses[masses < mass_upper]
            RT = _prediction.GetRetentionTime(pepinfo)
            
            if pepinfo in self.peptide_dict:
                item = self.peptide_dict[pepinfo]
                self.UpdateMassIntensity(item[0], item[1], masses, intens, RT if RT is not None else item[2])
            else:
                pepmass = pepmass / charge + self._ion_calc.base_mass.mass_proton
                self.InsertNewPeptide(pDeepFormat2PeptideModSeq(seq, mod), seq, pepmass, charge, masses, intens, RT if RT is not None else 0, peptide_to_protein_dict[seq] if seq in peptide_to_protein_dict else "")
            if count%10000 == 0:
                self.sql_conn.commit()
                print("[SQL UPDATE] {:.1f}%".format(100.0*count/len(_prediction.peptide_intensity_dict)), end="\r")
        self.sql_conn.commit()
        print("[SQL UPDATE] 100.0%: {}".format(self.dlib_file))
        
def GetMassIntensity(sql_conn, PeptideSeq, PrecursorCharge):
    cursor = sql_conn.execute("SELECT MassArray, IntensityArray FROM entries WHERE PeptideSeq = '%s' AND PrecursorCharge = %d"%(PeptideSeq, PrecursorCharge))
    row = cursor.fetchone()
    return DecodeMassList(row[0]), DecodeIntensityList(row[1])

def pDeepFormat2PeptideModSeq(seq, modinfo):
    if not modinfo: return seq
    moditems = modinfo.strip(";").split(";")
    modlist = []
    for moditem in moditems:
        site, mod = moditem.split(",")
        modlist.append((int(site), mod))
    modlist.sort(reverse=True)
    for site, mod in modlist:
        seq = seq[:site] + mod_dict[mod] + seq[site:]
    return seq

def PeptideModSeq2pDeepFormat(PeptideModSeq):
    site = PeptideModSeq.find('[')
    modlist = []
    while site != -1:
        if PeptideModSeq[site-1] == 'C': modlist.append('%d,%s;'%(site, 'Carbamidomethyl[C]'))
        elif PeptideModSeq[site-1] == 'M': modlist.append('%d,%s;'%(site, 'Oxidation[M]'))
        elif PeptideModSeq[site-1] == 'S': modlist.append('%d,%s;'%(site, 'Phospho[S]'))
        elif PeptideModSeq[site-1] == 'T': modlist.append('%d,%s;'%(site, 'Phospho[T]'))
        elif PeptideModSeq[site-1] == 'Y': modlist.append('%d,%s;'%(site, 'Phospho[Y]'))
        PeptideModSeq = PeptideModSeq[:site] + PeptideModSeq[PeptideModSeq.find(']')+1:]
        site = PeptideModSeq.find('[', site)
    return PeptideModSeq, "".join(modlist)

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
    
if __name__ == "__main__":
    dlib_obj = DLIB()
    dlib_obj.Open(r'e:\DIATools\encyclopedia\pan_human_subset.dlib')
    print ("Opened database successfully")
    
    peptide = 'AAAAAAAA'
    
    dlib_obj.InsertNewPeptide(peptide, peptide, 100, 2, [0]*10, [0]*10, 11)
    dlib_obj.InsertNewPeptide(peptide, peptide, 100, 3, [0]*10, [0]*10, 11)
    dlib_obj.sql_conn.commit()

    print ("Operation done successfully")
    dlib_obj.Close()