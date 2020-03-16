from math import log
import sys
import os
from shutil import copyfile
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

from math import log
from ..spectral_library.encyclopedia.dlib import DLIB
from ..spectral_library.openswath.tsv import OSW_TSV
from ..utils.mass_calc import PeptideIonCalculator
from ..spectral_library.library_base import SequenceLibrary 
from . import tune_and_predict
from ..sequence.protein_infer import infer_protein
from .generate_predicted_speclib import *
try:
    from ..pyRawFileReader.RawFileReader import RawFileReader
except:
    RawFileReader = None

if __name__ == "__main__":
    argd = {}
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i+1]
    
    speclib = argd['-library']
    
    seqlib = SequenceLibrary()
    peptide_list, pep_pro_dict = seqlib.PeptideListFromPeptideFile(argd['-peptide_file'])
        
    if '-raw' in argd:
        raw_path = argd['-raw']
        out_dir = os.path.split(raw_path)[0]
        if '-' in argd and RawFileReader:
            psmLabel, psmRT = Run_psmLabel(argd['-tune'], raw_path)
        else:
            psmLabel = ""
            psmRT = ""
    elif '-psmlabel' in argd:
        out_dir = os.path.split(psmLabel)
        psmLabel = argd['-psmlabel']
        if '-' in argd: psmRT = argd['-psmRT']
        else: psmRT = ""
    else:
        out_dir = argd['-dir']
        psmLabel = ""
        if '-psmRT' in argd: psmRT = argd['-psmRT']
        else: psmRT = ""
    speclib = os.path.join(out_dir, speclib)
    if speclib.endswith(".dlib"): copyfile('tmp/data/library/empty.dlib', speclib)
    elif speclib.endswith(".tsv"): copyfile('tmp/data/library/empty.tsv', speclib)
        
    pDeep_cfg = os.path.join(out_dir, "pDeep-tune.cfg")
    GenerateCFGpDeep(pDeep_cfg, psmLabel, psmRT)
    
    if psmLabel:
        SortPSM(psmLabel)
    
    if speclib.endswith(".dlib"):
        _lib = DLIB()
    elif speclib.endswith(".tsv"):
        _lib = OSW_TSV()
    _lib.Open(speclib)
    
    prediction = tune_and_predict.run(pDeep_cfg, peptide_list)
    
    _lib.UpdateByPrediction(prediction, pep_pro_dict)
    _lib.Close()
