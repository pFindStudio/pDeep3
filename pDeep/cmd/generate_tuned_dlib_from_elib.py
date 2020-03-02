from math import log

from ..pyRawFileReader.RawFileReader import RawFileReader
from ..spectral_library.encyclopedia.dlib import DLIB
from ..utils.mass_calc import PeptideIonCalculator
from . import tune_and_predict

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
    
def GenerateCFGpDeep(cfg_file, psmLabel):
    with open('tmp/predict/pDeep-tune-template.cfg') as f, open(cfg_file, "w") as out:
        cfg_str = "".join(f.readlines())
        out.write(cfg_str.format(psmLabel))

if __name__ == "__main__":
    import sys
    import os
    
    argd = {}
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i+1]
    
    elib_db = argd['-elib']
    dlib_db = argd['-dlib']
    raw_path = argd['-raw']
    
    PSMfile = elib_db+".psm.txt"
    rawName,_ = os.path.splitext(os.path.basename(raw_path))
    rawFile = RawFileReader(raw_path)
    
    elib = DLIB()
    elib.Open(elib_db)
    elib.GetAllPeptides()
    
    with open(PSMfile, "w") as output:
        output.write("raw_name\tscan\tpeptide\tmodinfo\tcharge\n")
        
        ion_calc = PeptideIonCalculator()
        
        for pepinfo, (_, charge, RT, _, _) in elib.peptide_dict.items():
            peptide, modinfo, _ = pepinfo.split("|")
            precursorMz = ion_calc.calc_pepmass(peptide, modinfo)
            precursorMz = precursorMz/charge + ion_calc.base_mass.mass_proton
            scan = FindMS2ScanNumFromPrecursorMzWithRTInSeconds(rawFile, precursorMz, RT)
            if scan:
                output.write("{}\t{}\t{}\t{}\t{}\n".format(rawName, scan, peptide, modinfo, charge))
            else:
                print("no scan found for {}".format(pepinfo))
    elib.Close()
    rawFile.Close()
                
    raw_filename = raw_path
    cfg_file = os.path.join(os.path.split(raw_filename)[0], "psmLabel.cfg")
    pDeep_cfg = os.path.join(os.path.split(raw_filename)[0], "pDeep-tune.cfg")
    print(raw_filename, cfg_file, PSMfile)
    
    GenerateCFGpsmLabel(cfg_file, PSMfile, raw_filename)
    os.chdir("psmLabel")
    os.system('psmLabel.exe "%s"'%cfg_file)
    os.chdir("..")
    
    psmLabel = os.path.splitext(raw_path)[0]+".psmLabel"
    GenerateCFGpDeep(pDeep_cfg, psmLabel)
    
    with open(psmLabel) as f:
        head_line = f.readline()
        head = head_line.strip().split("\t")
        headidx = dict(zip(head, range(len(head))))
        lines = f.readlines()
        
    def get_PSM_score(item, ion_type):
        item = item[headidx[ion_type]]
        if not item: return 0
        intens = []
        for ion_inten in item.strip(";").split(";"):
            if ion_inten[ion_inten.find('+')+1]!='0':
                intens.append(float(ion_inten.split(",")[1]))
        if len(intens) == 0: return 0
        return log(sum(intens))*len(intens)/len(item[headidx['peptide']])
        
    new_items = []
    for line in lines:
        item = line.split("\t")
        TIC = get_PSM_score(item, 'b')
        TIC += get_PSM_score(item, 'y')
        new_items.append( (TIC, line) )
    new_items.sort(key = lambda x: -x[0])
    with open(psmLabel, 'w') as f:
        f.write(head_line)
        f.writelines([item[1] for item in new_items])
    
    dlib = DLIB()
    dlib.Open(dlib_db)
    prediction = tune_and_predict.run(pDeep_cfg, dlib.GetAllPeptides())
    dlib.UpdateByPrediction(prediction)
    dlib.Close()
