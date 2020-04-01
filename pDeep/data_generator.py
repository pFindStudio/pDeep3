import os
from math import log

from .parameter import pDeepParameter
from .utils.mass_calc import PeptideIonCalculator
try:
    from .pyRawFileReader.RawFileReader import RawFileReader
except:
    RawFileReader = None

from .spectral_library.encyclopedia.dlib import DLIB
from .spectral_library.openswath.tsv import OSW_TSV
from .spectral_library.openswath.pqp import PQP,OSW
from .spectral_library.msp import MSP
from .search_engine.maxquant_reader import MaxQuantEvidenceReader as MQE
from .search_engine.pfind_reader import pFindSpectraReader as pFind
    
_library_writer_dict = {}
def _register_library_writer(ext, _class):
    _library_writer_dict[ext] = _class
    
def GetLibraryWriter(filename, pDeepParam):
    if os.path.isfile(filename):
        os.remove(filename)
    for ext, _class in _library_writer_dict.items():
        if filename.lower().endswith(ext): return _class(pDeepParam)
    return None
    
_register_library_writer('.dlib', DLIB)
_register_library_writer('.pqp', PQP)
_register_library_writer('.tsv', OSW_TSV)
_register_library_writer('.msp', MSP)
    
def GetRTInSecondsFromScanNum(rawFile, scan):
    return rawFile.RTInSecondsFromScanNum(scan)
    
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
    
def Generate_psmLabelCFG(cfg_file, input_PSM, raw_path):
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
    
def ReadModSeq(spikein_file):
    peptide_list = []
    pep_pro_dict = {}
    with open(spikein_file) as f:
        head = f.readline().strip().split("\t")
        headidx = dict(zip(head, range(len(head))))
        lines = f.readlines()
        for line in lines:
            items = line.strip().split("\t")
            seq = items[headidx['peptide']]
            mod = items[headidx['modinfo']]
            charge = items[headidx['charge']]
            protein = items[headidx['protein']]
            peptide_list.append((seq, mod, charge))
            pep_pro_dict[seq] = protein
    return peptide_list, pep_pro_dict
    
def Set_pDeepParam(param, instrument = "QE", ce = 27, psmLabel = "", psmRT = "", fixmod = "", varmod = "", n_tune=1000, psmLabel_test = "", threads = 4):
    param.predict_instrument = instrument
    param.predict_nce = ce
    if psmLabel: param.tune_psmlabels.append(psmLabel)
    param.tune_RT_psmlabel = psmRT
    if fixmod: param.fixmod.extend(fixmod.strip(",").split(","))
    if varmod: param.varmod.extend(varmod.strip(",").split(","))
    if psmLabel_test: param.test_psmlabels.append(psmLabel_test)
    param.n_tune_per_psmlabel = n_tune
    param.threads = threads
    return param
        
def GeneratePSMFile(result_file, raw_path):
    PSMfile = result_file+".psm.txt"
    rawName = os.path.splitext(os.path.basename(raw_path))[0]
    out_dir = os.path.split(raw_path)[0]
    rawFile = RawFileReader(raw_path)
    
    if result_file.lower().endswith(".elib"): 
        result = DLIB()
        multi_raw_in_result = False
    elif result_file.lower().endswith(".osw"): 
        result = OSW()
        multi_raw_in_result = False
    elif result_file.lower().endswith("evidence.txt"): 
        result = MQE()
        multi_raw_in_result = True
    elif result_file.lower().endswith(".spectra"):
        result = pFind()
        multi_raw_in_result = True
    result.Open(result_file)
    result.GetAllPeptides()
    
    with open(PSMfile, "w") as output:
        output.write("raw_name\tscan\tpeptide\tmodinfo\tcharge\tRTInSeconds\n")
        
        ion_calc = PeptideIonCalculator()
        
        for pepinfo, (_, charge, RT, raw, scan, protein) in result.peptide_dict.items():
            if multi_raw_in_result and raw != rawName: continue
            peptide, modinfo, _ = pepinfo.split("|")
            if scan != -1 and RT == -1:
                RT = GetRTInSecondsFromScanNum(rawFile, scan)
            elif scan == -1:
                precursorMz = ion_calc.calc_pepmass(peptide, modinfo)
                precursorMz = precursorMz/charge + ion_calc.base_mass.mass_proton
                scan = FindMS2ScanNumFromPrecursorMzWithRTInSeconds(rawFile, precursorMz, RT)
            if scan != -1:
                output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(rawName, scan, peptide, modinfo, charge, RT))
            else:
                print("no scan found for {}".format(pepinfo))
    result.Close()
    rawFile.Close()
    return PSMfile
        
def Run_psmLabel(PSMfile, raw_path):
    out_dir = os.path.split(raw_path)[0]
    
    cfg_file = os.path.join(out_dir, "psmLabel.cfg")
    print(raw_path, cfg_file, PSMfile)
    
    Generate_psmLabelCFG(cfg_file, PSMfile, raw_path)
    
    os.chdir("psmLabel")
    os.system('psmLabel.exe "%s"'%cfg_file)
    os.chdir("..")
    return os.path.splitext(raw_path)[0]+".psmlabel"
    
def Sort_psmLabel(psmLabel):
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
        