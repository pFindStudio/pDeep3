from math import log
import sys
import os
from shutil import copyfile

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
        with open(argd['-peptide']) as f:
            head = f.readline().strip().split("\t")
            headidx = dict(zip(head, range(len(head))))
            auto_protein = None if 'protein' in headidx else "pDeep"
            lines = f.readlines()
            for line in lines:
                items = line.strip().split("\t")
                fasta_peplist.append(tuple(items[:3]))
                if not auto_protein: pep_pro_dict[items[0]] = items[headidx['protein']]
                else: pep_pro_dict[items[0]] = auto_protein
        
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
    for seq, mod, charge in fasta_peplist:
        if "%s|%s|%d"%(seq, mod, charge) not in dlib.peptide_dict:
            peptide_list.append((seq, mod, charge))
    
    prediction = tune_and_predict.run(pDeep_cfg, peptide_list)
    
    dlib.UpdateByPrediction(prediction, pep_pro_dict)
    dlib.Close()
