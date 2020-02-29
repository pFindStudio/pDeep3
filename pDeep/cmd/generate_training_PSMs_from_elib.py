from ..pyRawFileReader.RawFileReader import RawFileReader
from ..spectral_library.encyclopedia.dlib import DLIB
from ..utils.mass_calc import PeptideIonCalculator

def CheckScanHasPrecursorMz(rawFile, scan, precursorMz):
    if scan < 1: return False
    elif rawFile.GetMSOrderForScanNum(scan) <= 1: return False
    iso_radius = rawFile.GetIsolationWidthForScanNum(scan)/2
    spectrumPrecursorMz = rawFile.GetPrecursorMassForScanNum(scan)
    if precursorMz >= spectrumPrecursorMz-iso_radius and precursorMz <= spectrumPrecursorMz+iso_radius:
        return True
    else:
        return False

def FindMS2ScanNumFromPrecursorMzWithRTInSeconds(rawFile, precursorMz, RT):
    scan = rawFile.ScanNumFromRTInSeconds(RT)
    if CheckScanHasPrecursorMz(rawFile, scan, precursorMz): return scan
    for i in range(1, 100):
        if CheckScanHasPrecursorMz(rawFile, scan+i, precursorMz): return scan+i
        elif CheckScanHasPrecursorMz(rawFile, scan-i, precursorMz): return scan-i
    return None
    
def GenerateCFGpsmLabel(cfg_file, input_PSM, raw_file):
    cfg_str =  "psm_type = none\n"
    cfg_str += "mode = pDeep\n"
    cfg_str += "num_psm_file = 1\n"
    cfg_str += "psm_file1 = %s\n"%input_PSM
    cfg_str += "ms2_type = raw\n"
    cfg_str += "num_ms2_file = 1\n"
    cfg_str += "ms2_file1 = %s\n"%raw_file
    cfg_str += "output_folder = %s\n"%os.path.split(raw_file)[0]
    cfg_str += "NH3_loss = true\n"
    cfg_str += "H2O_loss = true\n"
    cfg_str += "Mod_loss = true\n"
    cfg_str += "num_ion_type = 2\n"
    cfg_str += "iontype1 = b|N_term|0\n"
    cfg_str += "iontype2 = y|C_term|0\n"
    cfg_str += "num_new_aa = 0"
    with open(cfg_file, "w") as f: f.write(cfg_str)

if __name__ == "__main__":
    import sys
    import os
    rawName,_ = os.path.splitext(os.path.basename(sys.argv[1]))
    rawFile = RawFileReader(sys.argv[1])
    
    dlib = DLIB()
    dlib.Open(sys.argv[2])
    dlib.GetAllPeptides()
    
    with open(sys.argv[3], "w") as output:
        output.write("raw_name\tscan\tpeptide\tmodinfo\tcharge\n")
        
        ion_calc = PeptideIonCalculator()
        
        for pepinfo, (_, charge, RT, _, _) in dlib.peptide_dict.items():
            peptide, modinfo, _ = pepinfo.split("|")
            precursorMz = ion_calc.calc_pepmass(peptide, modinfo)
            precursorMz = precursorMz/charge + ion_calc.base_mass.mass_proton
            scan = FindMS2ScanNumFromPrecursorMzWithRTInSeconds(rawFile, precursorMz, RT)
            if scan:
                output.write("{}\t{}\t{}\t{}\t{}\n".format(rawName, scan, peptide, modinfo, charge))
            else:
                print("no scan found for {}".format(pepinfo))
                
    #run psmLabel
    raw_filename = sys.argv[1]
    cfg_file = os.path.join(os.path.split(raw_filename)[0], "psmLabel.cfg")
    input_PSM = sys.argv[3]
    print(raw_filename, cfg_file, input_PSM)
    
    GenerateCFGpsmLabel(cfg_file, input_PSM, raw_filename)
    os.chdir("psmLabel")
    os.system('psmLabel.exe "%s"'%cfg_file)
    os.chdir("..")
    
