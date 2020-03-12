from math import log
import sys
import os
from shutil import copyfile
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

from math import log
from ..spectral_library.encyclopedia.dlib import DLIB
from ..utils.mass_calc import PeptideIonCalculator
from ..spectral_library.library_base import SequenceLibrary 
from . import tune_and_predict
from ..sequence.protein_infer import infer_protein
from .generate_tuned_dlib_from_elib import *
try:
    from ..pyRawFileReader.RawFileReader import RawFileReader
except:
    RawFileReader = None

if __name__ == "__main__":
    argd = {}
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i+1]
    
    dlib_db = argd['-dlib']
    
    fasta_peplist = []
    pep_pro_dict = {}
    if '-peptide' in argd:
        seqlib = SequenceLibrary()
        peptide_list, pep_pro_dict = seqlib.PeptideListFromPeptideFile(argd['-peptide'])
        
    if '-raw' in argd:
        raw_path = argd['-raw']
        raw_dir = os.path.split(raw_path)[0]
        if '-elib' in argd and RawFileReader:
            psmLabel, psmRT = Run_psmLabel(argd['-elib'], raw_path)
        else:
            psmLabel = ""
            psmRT = ""
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
    dlib_db = os.path.join(raw_dir, dlib_db)
    copyfile('tmp/data/library/empty.dlib', dlib_db)
        
    pDeep_cfg = os.path.join(raw_dir, "pDeep-tune.cfg")
    GenerateCFGpDeep(pDeep_cfg, psmLabel, psmRT)
    
    if psmLabel:
        SortPSM(psmLabel)
    
    dlib = DLIB()
    dlib.Open(dlib_db)
    
    prediction = tune_and_predict.run(pDeep_cfg, peptide_list)
    
    dlib.UpdateByPrediction(prediction, pep_pro_dict)
    dlib.Close()
