from math import log
import sys
import os
from shutil import copyfile
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

from ..spectral_library.encyclopedia.dlib import DLIB
from ..spectral_library.openswath.tsv import OSW_TSV
from ..spectral_library.openswath.pqp import PQP,OSW
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
        
def Run_psmLabel(result_file, raw_path):
    PSMfile = result_file+".psm.txt"
    rawName = os.path.splitext(os.path.basename(raw_path))[0]
    out_dir = os.path.split(raw_path)[0]
    rawFile = RawFileReader(raw_path)
    
    if result_file.endswith(".elib"): result = DLIB()
    elif result_file.endswith(".osw"): result = OSW()
    result.Open(result_file)
    result.GetAllPeptides()
    
    with open(PSMfile, "w") as output:
        output.write("raw_name\tscan\tpeptide\tmodinfo\tcharge\tRTInSeconds\n")
        
        ion_calc = PeptideIonCalculator()
        
        for pepinfo, (_, charge, RT, _) in result.peptide_dict.items():
            peptide, modinfo, _ = pepinfo.split("|")
            precursorMz = ion_calc.calc_pepmass(peptide, modinfo)
            precursorMz = precursorMz/charge + ion_calc.base_mass.mass_proton
            scan = FindMS2ScanNumFromPrecursorMzWithRTInSeconds(rawFile, precursorMz, RT)
            if scan:
                output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(rawName, scan, peptide, modinfo, charge, RT))
            else:
                print("no scan found for {}".format(pepinfo))
    result.Close()
    rawFile.Close()
    
    cfg_file = os.path.join(out_dir, "psmLabel.cfg")
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
    import argparse
    argd = {}
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i].lower()] = sys.argv[i+1]
    
    out_lib = argd['-out_library']
    out_dir = os.path.split(out_lib)[0]
    if '-decoy' in argd: decoy = argd['-decoy']
    else: decoy = "pseudo_reverse"
    
    if '-fasta' in argd:
        if '-proteins' in argd:
            protein_list = argd['-proteins'].split(",")
        else:
            protein_list = None
        seqlib = SequenceLibrary(min_precursor_mz = 400, max_precursor_mz = 1200)
        fasta_peplist, protein_dict = seqlib.PeptideListFromFasta(argd['-fasta'], protein_list)
    else:
        fasta_peplist, protein_dict = [], {}
        
    if '-tune_psm' in argd and RawFileReader:
        raw_path = argd['-raw']
        psmLabel, psmRT = Run_psmLabel(argd['-tune_psm'], raw_path)
    elif '-psmlabel' in argd:
        psmLabel = argd['-psmlabel']
        if '-psmRT' in argd: psmRT = argd['-psmRT']
        else: psmRT = ""
    else:
        psmLabel = ""
        if '-psmRT' in argd: psmRT = argd['-psmRT']
        else: psmRT = ""

    copyfile('tmp/data/library/empty'+os.path.splitext(out_lib)[-1], out_lib)
        
    pDeep_cfg = os.path.join(out_dir, "pDeep-tune.cfg")
    GenerateCFGpDeep(pDeep_cfg, psmLabel, psmRT)
    
    if psmLabel:
        SortPSM(psmLabel)
        
    if '-transition_tsv' in argd:
        speclib_file = argd['-transition_tsv']
        tsv = OSW_TSV()
    
        tsv.Open(speclib_file)
        
        peptide_list = tsv.GetAllPeptides()
        for seq, mod, charge in peptide_list:
            if "%s|%s|%d"%(seq, mod, charge) not in tsv.peptide_dict:
                peptide_list.append((seq, mod, charge))
        
        tsv.Close()
    else:
        peptide_list = fasta_peplist
    
    prediction = tune_and_predict.run(pDeep_cfg, peptide_list)
    
    pep_pro_dict = infer_protein([seq for seq, mod, charge in fasta_peplist], protein_dict)
    fasta_peptopro_dict = dict([(peptide, ";".join([pro_ac for pro_ac, site in prosites])) for peptide, prosites in pep_pro_dict.items()])
    
    if out_lib.endswith(".dlib"):
        _lib = DLIB()
    elif out_lib.endswith(".tsv"):
        _lib = OSW_TSV()
    elif out_lib.endswith(".pqp"):
        _lib = PQP()
    _lib.Open(out_lib)
    
    _lib.decoy = decoy
    _lib.UpdateByPrediction(prediction, fasta_peptopro_dict, min_intensity = 0.1, least_n_peaks = 6, max_mz = 2000)
    _lib.Close()
