import sqlite3
import zlib
import struct
import numpy as np

from ..mass_calc import PeptideIonCalculator as ion_calc
# from ...prediction import pDeepPrediction as prediction

# pack/unpack(fmt,...), fmt=">" means big-endian

class DLIB:
    def __init__(self):
        self.sql_conn = None
        self.mod_dict = {
            "Carbamidomethyl[C]": "[+57.021464]",
            "Oxidation[M]": "[+15.994915]",
        }
        self._ion_calc = ion_calc()
        
    def Open(self, dlib_file):
        self.sql_conn = sqlite3.connect(dlib_file)
        
    def Close(self):
        self.sql_conn.close()
        
    def GetAllPeptides(self):
        self.peptide_dict = {}
        peptide_list = []
        cursor = self.sql_conn.execute("SELECT PeptideModSeq, PrecursorCharge FROM entries")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            charge = int(row[1])
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = [row[0], charge, None, None] #items = [PeptideModSeq, PrecursorCharge, PredictedMassList, PredictedIntensityList]
        return peptide_list
        
    def UpdateMassIntensity(self, PeptideModSeq, PrecursorCharge, MassList, IntensityList, commit_now = False):
        sql = "UPDATE entries SET MassEncodedLength = ?, MassArray = ?, IntensityEncodedLength = ?, IntensityArray = ? WHERE PeptideModSeq = '%s' AND PrecursorCharge = %d"%(PeptideModSeq, PrecursorCharge)
        self.sql_conn.execute(sql, (len(MassList)*8, EncodeMassList(MassList), len(IntensityList)*4, EncodeIntensityList(IntensityList)))
        if commit_now: sql_conn.commit()
    
    def UpdateByPeptideDict(self, peptide_dict):
        for pepinfo, (modseq, charge, mass, intensity) in peptide_dict.items():
            if not mass: continue
            self.UpdateMassIntensity(modseq, charge, mass, intensity)
        self.sql_conn.commit()
    
    def UpdateByPrediction(self, _prediction, intensity_threshold = 0.1):
        count = 0
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            count += 1
            seq, mod, charge = pepinfo.split("|")
            masses = self._ion_calc.calc_b_y_ions(seq, mod, 2)
            # print(seq, mod, masses)
            intens = intensities[:,:masses.shape[1]]
            intens[0, 0:2] = 0 #b1+/b1++ = 0
            intens[-1, 2:] = 0 #y1+/y1++ = 0, do not consider y1/b1 in the library
            masses = masses.reshape(-1)
            intens = intens.reshape(-1)
            
            masses = masses[intens > intensity_threshold]
            intens = intens[intens > intensity_threshold]*10000
            
            item = self.peptide_dict[pepinfo]
            self.UpdateMassIntensity(item[0], item[1], masses, intens)
            if count%10000 == 0:
                self.sql_conn.commit()
                print("[SQL UPDATE] {:.1f}%".format(100.0*count/len(_prediction.peptide_intensity_dict)), end="\r")
        self.sql_conn.commit()
        print("[SQL UPDATE] 100.0%")
        
def GetMassIntensity(sql_conn, PeptideSeq, PrecursorCharge):
    cursor = sql_conn.execute("SELECT MassArray, IntensityArray FROM entries WHERE PeptideSeq = '%s' AND PrecursorCharge = %d"%(PeptideSeq, PrecursorCharge))
    row = cursor.fetchone()
    return DecodeMassList(row[0]), DecodeIntensityList(row[1])

def PeptideModSeq2pDeepFormat(PeptideModSeq):
    site = PeptideModSeq.find('[')
    modlist = []
    while site != -1:
        if PeptideModSeq[site-1] == 'C': modlist.append('%d,%s;'%(site, 'Carbamidomethyl[C]'))
        elif PeptideModSeq[site-1] == 'M': modlist.append('%d,%s;'%(site, 'Oxidation[M]'))
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
    dlib_obj.Open(r'c:\DataSets\DIA\DIA-Tools\pan_human_library.dlib')
    print ("Opened database successfully")
    
    peptide = 'ILITIVEEVETLR'
    charge2 = 2
    charge3 = 3
    def print_mass_inten(dlib_obj, peptide, charge):
        mass, intensity = dlib_obj.GetMassIntensity(peptide, charge)
        print(mass)
        print(intensity)
    print_mass_inten(dlib_obj, peptide, charge2)
    print_mass_inten(dlib_obj, peptide, charge3)
    # pep_dict = {
        # "a": [peptide, charge2, [0]*10, [0]*10],
        # "b": [peptide, charge3, [1]*10, [1]*10],
    # }
    pep_dict = {
        "a": [peptide, charge2, (389.250694684, 441.30714692, 518.293287778, 554.391210902, 617.361701696, 651.36919237, 746.40429479, 875.446887884, 974.515301802, 1087.59936578, 1188.64704426, 1301.73110824),(2175.5, 2418.60009765625, 2238.39990234375, 1386.699951171875, 2182.0, 1083.199951171875, 2944.5, 4886.7998046875, 9557.400390625, 3463.89990234375, 6086.39990234375, 2990.199951171875)],
        "b": [peptide, charge3, (389.250694684, 441.30714692, 505.810250713, 518.293287778, 554.391210902, 617.361701696, 653.45962482, 746.40429479, 782.502217914, 875.446887884, 911.544811008, 974.515301802, 1087.59936578),(5543.89990234375, 3531.0, 2671.800048828125, 9054.400390625, 2884.5, 4665.2001953125, 2104.300048828125, 3087.60009765625, 1850.199951171875, 3758.60009765625, 3443.800048828125, 3082.800048828125, 1543.9000244140625)],
    }
    UpdateByPeptideDict(conn, pep_dict)
    # UpdateMassIntensity(conn, peptide, charge, [0]*10, [0]*10)
    # UpdateMassIntensity(conn, peptide, charge, (389.250694684, 441.30714692, 518.293287778, 554.391210902, 617.361701696, 651.36919237, 746.40429479, 875.446887884, 974.515301802, 1087.59936578, 1188.64704426, 1301.73110824),(2175.5, 2418.60009765625, 2238.39990234375, 1386.699951171875, 2182.0, 1083.199951171875, 2944.5, 4886.7998046875, 9557.400390625, 3463.89990234375, 6086.39990234375, 2990.199951171875))

    print ("Operation done successfully")
    dlib_obj.Close()