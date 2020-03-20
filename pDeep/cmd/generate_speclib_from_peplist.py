from math import log
import sys
import os
from shutil import copyfile
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

from math import log
from ..spectral_library.encyclopedia.dlib import DLIB
from ..spectral_library.openswath.tsv import OSW_TSV
from ..spectral_library.openswath.pqp import PQP
from ..utils.mass_calc import PeptideIonCalculator
from ..spectral_library.library_base import SequenceLibrary 
from . import tune_and_predict
from ..sequence.protein_infer import infer_protein
from ..data_generator import *

if __name__ == "__main__":
    argd = {}
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i+1]
    
    speclib = argd['--output']
    out_dir = os.path.split(speclib)[0]
    
    if '--decoy' in argd: decoy = argd['--decoy']
    else: decoy = "reverse"
    
    seqlib = SequenceLibrary()
    peptide_list, pep_pro_dict = seqlib.PeptideListFromPeptideFile(argd['--input'])
    if '--raw' in argd:
        raw_path = argd['--raw']
        if '--tune_psm' in argd and RawFileReader:
            psmLabel, psmRT = Run_psmLabel(argd['--tune_psm'], raw_path)
        else:
            psmLabel = ""
            psmRT = ""
    elif '--psmlabel' in argd:
        psmLabel = argd['--psmlabel']
        if '--psmRT' in argd: psmRT = argd['--psmRT']
        else: psmRT = ""
    else:
        psmLabel = ""
        if '--psmRT' in argd: psmRT = argd['--psmRT']
        else: psmRT = ""
        
    copyfile('tmp/data/library/empty'+os.path.splitext(speclib)[-1], speclib)
        
    pDeep_cfg = os.path.join(out_dir, "pDeep-tune.cfg")
    GenerateCFGpDeep(pDeep_cfg, psmLabel, psmRT)
    
    if psmLabel:
        SortPSM(psmLabel)
    
    if speclib.endswith(".dlib"):
        _lib = DLIB()
    elif speclib.endswith(".tsv"):
        _lib = OSW_TSV()
    elif speclib.endswith(".pqp"):
        _lib = PQP()
    _lib.Open(speclib)
    
    prediction = tune_and_predict.run(pDeep_cfg, peptide_list)
    
    _lib.decoy = decoy
    _lib.UpdateByPrediction(prediction, pep_pro_dict)
    _lib.Close()
