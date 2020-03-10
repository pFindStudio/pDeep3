from math import log
import sys
import os
from shutil import copyfile

from ..spectral_library.encyclopedia.dlib import DLIB
from ..utils.mass_calc import PeptideIonCalculator
from ..spectral_library.library_base import SequenceLibrary 
from . import tune_and_predict
from ..sequence.protein_infer import infer_protein
try:
    from ..pyRawFileReader.RawFileReader import RawFileReader
except:
    RawFileReader = None

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
    
def GenerateCFGpsmLabel(cfg_file, input_PSM, raw_path):
    cfg_str =  "psm_type = none\n"
    cfg_str += "mode = pDeep\n"
    cfg_str += "num_psm_file = 1\n"
    cfg_str += "psm_file1 = %s\n"%input_PSM
    cfg_str += "ms2_type = raw\n"
    cfg_str += "num_ms2_file = 1\n"
    cfg_str += "ms2_file1 = %s\n"%raw_path
    cfg_str += "output_folder = %s\n"%(os.path.split(raw_path)[0])
    cfg_str += "NH3_loss = true\n"
    cfg_str += "H2O_loss = true\n"
    cfg_str += "Mod_loss = true\n"
    cfg_str += "num_ion_type = 2\n"
    cfg_str += "iontype1 = b|N_term|0\n"
    cfg_str += "iontype2 = y|C_term|0\n"
    cfg_str += "num_new_aa = 0"
    with open(cfg_file, "w") as f: f.write(cfg_str)
    
def GenerateCFGpDeep(cfg_file, psmLabel, psmRT):
    with open('tmp/predict/pDeep-tune-template.cfg') as f, open(cfg_file, "w") as out:
        lines = f.readlines()
        cfg_str = "".join(lines)
        out.write(cfg_str.format(psmLabel, psmRT))
        
def Run_psmLabel(elib_db, raw_path):
    PSMfile = elib_db+".psm.txt"
    rawName = os.path.splitext(os.path.basename(raw_path))[0]
    raw_dir = os.path.split(raw_path)[0]
    rawFile = RawFileReader(raw_path)
    
    elib = DLIB()
    elib.Open(elib_db)
    elib.GetAllPeptides()
    
    with open(PSMfile, "w") as output:
        output.write("raw_name\tscan\tpeptide\tmodinfo\tcharge\tRTInSeconds\n")
        
        ion_calc = PeptideIonCalculator()
        
        for pepinfo, (_, charge, RT, _, _) in elib.peptide_dict.items():
            peptide, modinfo, _ = pepinfo.split("|")
            precursorMz = ion_calc.calc_pepmass(peptide, modinfo)
            precursorMz = precursorMz/charge + ion_calc.base_mass.mass_proton
            scan = FindMS2ScanNumFromPrecursorMzWithRTInSeconds(rawFile, precursorMz, RT)
            if scan:
                output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(rawName, scan, peptide, modinfo, charge, RT))
            else:
                print("no scan found for {}".format(pepinfo))
    elib.Close()
    rawFile.Close()
    
    cfg_file = os.path.join(raw_dir, "psmLabel.cfg")
    print(raw_path, cfg_file, PSMfile)
    
    GenerateCFGpsmLabel(cfg_file, PSMfile, raw_path)
    
    os.chdir("psmLabel")
    os.system('psmLabel.exe "%s"'%cfg_file)
    os.chdir("..")
    return os.path.splitext(raw_path)[0]+".psmlabel", PSMfile
    
def SortPSM(psmLabel):
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

if __name__ == "__main__":
    argd = {}
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i].lower()] = sys.argv[i+1]
    
    dlib_db = argd['-dlib']
    
    if '-fasta' in argd:
        seqlib = SequenceLibrary(min_precursor_mz = 400, max_precursor_mz = 1000)
        fasta_peplist, protein_dict = seqlib.PeptideListFromFasta(argd['-fasta'])
    else:
        fasta_peplist, protein_dict = [], {}
        
    if '-elib' in argd and RawFileReader:
        raw_path = argd['-raw']
        raw_dir = os.path.split(raw_path)[0]
        psmLabel, psmRT = Run_psmLabel(argd['-elib'], raw_path)
    elif '-psmlabel' in argd:
        raw_dir = os.path.split(psmLabel)
        psmLabel = argd['-psmlabel']
        if '-psmRT' in argd: psmRT = argd['-psmRT']
        else: psmRT = ""
    else:
        raw_dir = argd['-dir']
        psmLabel = ""
        if '-psmRT' in argd: psmRT = argd['-psmRT']
        else: psmRT = ""
    
    new_dlib = os.path.join(raw_dir, "pdeep_tune.dlib")
    copyfile(dlib_db, new_dlib)
    dlib_db = new_dlib
        
    pDeep_cfg = os.path.join(raw_dir, "pDeep-tune.cfg")
    GenerateCFGpDeep(pDeep_cfg, psmLabel, psmRT)
    
    if psmLabel:
        SortPSM(psmLabel)
        
    dlib = DLIB()
    dlib.Open(dlib_db)
    
    peptide_list = dlib.GetAllPeptides()
    for pepinfo in fasta_peplist:
        if pepinfo not in dlib.peptide_dict:
            peptide_list.append(pepinfo)
    
    prediction = tune_and_predict.run(pDeep_cfg, peptide_list)
    
    pep_pro_dict = infer_protein([seq for seq, mod, charge in fasta_peplist], protein_dict)
    fasta_peptopro_dict = dict([(peptide, ";".join([pro_ac for pro_ac, site in prosites])) for peptide, prosites in pep_pro_dict.items()])
    
    dlib.UpdateByPrediction(prediction, fasta_peptopro_dict)
    dlib.Close()
