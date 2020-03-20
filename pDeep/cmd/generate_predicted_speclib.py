import sys
import os
from shutil import copyfile
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

from ..spectral_library.library_base import SequenceLibrary 
from ..sequence.protein_infer import infer_protein
from . import tune_and_predict
from ..data_generator import *

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generating spectral library for OpenSWATH and EncyclopDIA by pDeep.')
    
    parser.add_argument('--input', type=str, required=True, help='.peplib (text) file for peptide library; or .fasta for generating peptides from fasta; or .tsv file for transitions/assays file (for example pan-human library).')
    parser.add_argument('--output', type=str, required=True, help='The genereted library file.')
    
    parser.add_argument('--varmod', type=str, default="Oxidation[M]", required=False, help='Variable modifications, seperated by ",".')
    parser.add_argument('--fixmod', type=str, default="Carbamidomethyl[C]", required=False, help='Fixed modifications, seperated by ",".')
    
    parser.add_argument('--proteins', type=str, required=False, help='Only considering these proteins when input is fasta, seperated by "," (default: use all proteins).')
    
    
    parser.add_argument('--spikein', type=str, required=False, help='Spike-in file could be: .peplib file for peptide library; or .fasta for generating peptides from fasta; or .tsv file for transition or assays file (for example pan-human library), or .txt spike-in file with title "peptide", "modinfo", "charge", and "protein" in pDeep format, pDeep will not add modifications to peptides from this .txt file.')
    parser.add_argument('--spikein_proteins', type=str, required=False, help='Only considering these proteins when spikein is fasta, seperated by "," (default: use all proteins).')
    # parser.add_argument('--spikein_mod', type=str, default="", required=False, help='Modifications within spike-in peptides, all of them are considered as fixed, seperated by ","')
    
    
    parser.add_argument('--instrument', type=str, default="QE", required=False, help='Instrument type for prediction.')
    parser.add_argument('--ce', type=float, default=27, required=False, help='Collision energy for prediction.')
    
    parser.add_argument('--decoy', type=str, choices=['reverse','pseudo_reverse'], default='reverse', help='Decoy method when generating OSW PQP file.')
    
    parser.add_argument('--tune_psm', type=str, required=False, help='.osw (OpenSWATH), .elib (EncyclopDIA), evidence.txt (MaxQuant) or .spectra (pFind) file for tuning pDeep and pDeepRT.')
    parser.add_argument('--raw', type=str, required=False, help='Raw file for tuning pDeep and pDeepRT.')
    
    parser.add_argument('--psmlabel', type=str, default="", required=False, help='psmLabel file for tuning pDeep, could be generated by psmLabel.exe.')
    parser.add_argument('--psmRT', type=str, default="", required=False, help='psm file for tuning pDeepRT, must containing a column with title "RT" or "RTInSeconds".')
    
    args = parser.parse_args()
    
    
    out_lib = args.output
    out_dir = os.path.split(out_lib)[0]
    decoy = args.decoy
    
    copyfile('tmp/data/library/empty'+os.path.splitext(out_lib)[-1], out_lib)
    
    _lib = GetLibraryWriter(out_lib)
    _lib.Open(out_lib)
    _lib.decoy = decoy
    
    seqlib = SequenceLibrary(min_charge = 2, max_charge = 4, min_precursor_mz = 400, max_precursor_mz = 1200, varmod=args.varmod, fixmod=args.fixmod)
    
    def _from_fasta(fasta, proteins = None):
        if proteins:
            protein_list = proteins.split(",")
        else:
            protein_list = None
        peptide_list, protein_dict = seqlib.PeptideListFromFasta(args.input, protein_list)
        infer_pep_pro_dict = infer_protein([seq for seq, mod, charge in peptide_list], protein_dict)
        pep_pro_dict = dict([(peptide,"/".join([pro_ac for pro_ac, site in prosites])) for peptide, prosites in infer_pep_pro_dict.items()])
        return peptide_list, pep_pro_dict
        
    def _from_tsv(tsvfile):
        tsv = OSW_TSV()
        tsv.Open(tsvfile)
        peptide_list = tsv.GetAllPeptides()
        pep_pro_dict = dict([(pepinfo.split("|")[0], item[-1]) for pepinfo, item in tsv.peptide_dict.items()])
        tsv.Close()
        return peptide_list, pep_pro_dict
        
    if args.input.endswith('.peplib'):
        peptide_list, pep_pro_dict = seqlib.PeptideListFromPeptideFile(args.input)
    elif args.input.endswith('.fasta'):
        peptide_list, pep_pro_dict = _from_fasta(args.input, args.proteins)
    elif args.input.endswith('.tsv'):
        peptide_list, pep_pro_dict = _from_tsv(args.input)
    else:
        peptide_list, pep_pro_dict = seqlib.PeptideListFromPeptideFile(args.input)
    
    if args.spikein:
        if args.spikein.endswith('.txt'):
            spkin_list, spkin_dict = ReadSpikein(args.spikein)
        elif args.input.endswith('.peplib'):
            spkin_list, spkin_dict = seqlib.PeptideListFromPeptideFile(args.spikein)
        elif args.input.endswith('.fasta'):
            spkin_list, spkin_dict = _from_fasta(args.spikein, args.spikein_proteins)
        elif args.input.endswith('.tsv'):
            spkin_list, spkin_dict = _from_tsv(args.spikein)
        else:
            spkin_list, spkin_dict = ReadSpikein(args.spikein)
            
        peptide_list.extend(spkin_list)
        for seq, pro in spkin_dict.items():
            if seq not in pep_pro_dict:
                pep_pro_dict[seq] = pro
        for seq,mod,charge in spkin_list:
            if mod:
                for modname in mod.strip(";").split(";"):
                    modname = modname[modname.find(",")+1:]
                    if modname not in args.fixmod and modname not in args.varmod:
                        args.fixmod += ","+modname
        args.fixmod = args.fixmod.strip(',')
        print("[pDeep Info] fix modification included in spike-in peptides: '%s'"%args.fixmod)
        
    if args.tune_psm and args.raw and RawFileReader:
        raw_path = args.raw
        psmRT = GeneratePSMFile(args.tune_psm, raw_path)
        psmLabel = Run_psmLabel(psmRT, raw_path)
    else:
        psmLabel = args.psmlabel
        psmRT = args.psmRT
    if psmLabel:
        Sort_psmLabel(psmLabel)
        
    pDeep_cfg = Generate_pDeepParam(instrument=args.instrument, ce=args.ce, psmLabel=psmLabel, psmRT=psmRT, fixmod=args.fixmod, varmod=args.varmod)
    
    prediction = tune_and_predict.run(pDeep_cfg, peptide_list)
    
    _lib.UpdateByPrediction(prediction, pep_pro_dict, min_intensity = 0.1, least_n_peaks = 6, max_mz = 2000)
    _lib.Close()