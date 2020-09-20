import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from shutil import copyfile

from ..spectral_library.library_base import SequenceLibrary 
from ..sequence.protein_infer import infer_protein
from ..data_generator import *
from ..parameter import pDeepParameter
from . import tune_and_predict

def _from_fasta(_seqlib, fasta, proteins = None):
    if proteins:
        protein_list = proteins.split(",")
    else:
        protein_list = None
    peptide_list, protein_dict = _seqlib.PeptideListFromFasta(fasta, protein_list)
    infer_pep_pro_dict = infer_protein([seq for seq, mod, charge in peptide_list], protein_dict)
    pep_pro_dict = dict([(peptide,"/".join([pro_ac for pro_ac, site in prosites])) for peptide, prosites in infer_pep_pro_dict.items()])
    return peptide_list, pep_pro_dict
    
def _read_lib(lib_reader, infile):
    lib_reader.Open(infile)
    peptide_list = lib_reader.GetAllPeptides()
    pep_pro_dict = dict([(pepinfo.split("|")[0], item[-1]) for pepinfo, item in lib_reader.peptide_dict.items()])
    lib_reader.Close()
    return peptide_list, pep_pro_dict
    
def _from_tsv(tsvfile):
    tsv = OSW_TSV()
    return _read_lib(tsv, tsvfile)

def _from_pqp(pqpfile):
    pqp = PQP()
    return _read_lib(pqp, pqpfile)
    
def _from_dlib(dlibfile):
    dlib = DLIB()
    return _read_lib(dlib, dlibfile)
    
def _from_MaxQuant(txtfile):
    mq = MQE()
    return _read_lib(mq, txtfile)
    
def _from_pFind(txtfile):
    pf = pFind()
    return _read_lib(pf, txtfile)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generating spectral library for OpenSWATH and EncyclopDIA by pDeep.')
    
    parser.add_argument('--input', type=str, required=True, help='.peplib (text) file for peptide library; or .fasta for generating peptides from fasta; or .tsv file for transitions/assays file (for example pan-human library).')
    parser.add_argument('--output', type=str, required=True, help='The genereted library file. File type is determined by the file extention: .dlib for EncyclopeDIA (SQLite file), .pqp for OpenSWATH (SQLite file), .tsv for OpenSWATH (text file), .msp for generic text spectrum file.')
    
    parser.add_argument('--varmod', type=str, default="Oxidation[M]", required=False, help='Variable modifications, seperated by ",".')
    parser.add_argument('--fixmod', type=str, default="Carbamidomethyl[C]", required=False, help='Fixed modifications, seperated by ",".')
    
    parser.add_argument('--target_proteins', type=str, required=False, help='Only considering these proteins (ACs or uniprot IDs) when input is fasta, seperated by "," (default: use all proteins).')
    
    
    parser.add_argument('--spikein', type=str, required=False, help='Spike-in file could be: .peplib file for peptide library; or .fasta for generating peptides from fasta; or .tsv file for transition or assays file (for example pan-human library), or .txt spike-in file with title "peptide", "modinfo", "charge", and "protein" in pDeep format, pDeep will not add modifications to peptides from this .txt file.')
    parser.add_argument('--spikein_proteins', type=str, required=False, help='Only considering these proteins when spikein is fasta, seperated by "," (default: use all proteins).')
    parser.add_argument('--spikein_varmod', type=str, default="Oxidation[M]", required=False, help='Variable modifications within spike-in peptides, seperated by ","')
    parser.add_argument('--spikein_fixmod', type=str, default="Carbamidomethyl[C]", required=False, help='Fix modifications within spike-in peptides, seperated by ","')
    
    
    parser.add_argument('--min_intensity', type=float, default=0.1, required=False, help='Min intensity of fragments in the library')
    parser.add_argument('--least_n_peaks', type=int, default=6, required=False, help='Keep at least n fragments in the library')
    parser.add_argument('--min_mz', type=float, default=300, required=False, help='Min mz of fragments in the library')
    parser.add_argument('--max_mz', type=float, default=2000, required=False, help='Max mz of fragments in the library')
    parser.add_argument('--ion_type', type=str, default="b,y,b-ModLoss,y-ModLoss", required=False, help='Ion types in the library')
    parser.add_argument('--model', type=str, default="HCD", required=False, help='The model file or "HCD", "ETchD" ...')
    
    parser.add_argument('--min_precursor_charge', type=int, default=2, required=False, help='Min precursor charge of peptides in the library')
    parser.add_argument('--max_precursor_charge', type=int, default=4, required=False, help='Max precursor charge of peptides in the library')
    parser.add_argument('--min_precursor_mz', type=float, default=400, required=False, help='Min precursor mz of peptides in the library')
    parser.add_argument('--max_precursor_mz', type=float, default=1200, required=False, help='Max precursor mz of peptides in the library')
    parser.add_argument('--min_peptide_length', type=int, default=6, required=False, help='Min peptide length of peptides in the library')
    parser.add_argument('--max_peptide_length', type=int, default=60, required=False, help='Max peptide length of peptides in the library')
    
    
    parser.add_argument('--instrument', type=str, default="QE", required=False, help='Instrument type for prediction.')
    parser.add_argument('--ce', type=float, default=27, required=False, help='Collision energy for prediction.')
    
    parser.add_argument('--decoy', type=str, choices=['reverse','pseudo_reverse','no_decoy'], default='reverse', help='Decoy method when generating OpenSWATH PQP file.')
    
    parser.add_argument('--tune_psm', type=str, required=False, help='.osw (OpenSWATH), .elib (EncyclopDIA), evidence.txt (MaxQuant), .spectra (pFind) or .txt (tab seperated file with title "raw_name, scan, peptide, modinfo, charge, RTInSeconds" file for tuning pDeep and pDeepRT.')
    parser.add_argument('--n_tune_psm', type=int, default=1000, required=False, help='number of PSMs for tuning.')
    parser.add_argument('--raw', type=str, required=False, help='Raw file for tuning pDeep and pDeepRT.')
    
    parser.add_argument('--psmlabel', type=str, default="", required=False, help='psmLabel file for tuning pDeep, could be generated by psmLabel.exe.')
    parser.add_argument('--psmRT', type=str, default="", required=False, help='psm file for tuning pDeepRT, must containing a column with title "RT" or "RTInSeconds".')
    
    parser.add_argument('--RT_input', type=str, required=False, help='Input for RT normalization for OpenSWATH.')
    parser.add_argument('--RT_proteins', type=str, required=False, help='Protein ACs for RT normalization for OpenSWATH (if RT_input is fasta).')
    parser.add_argument('--RT_tsv', type=str, required=False, help='Output tsv for RT normalization for OpenSWATH.')
    
    args = parser.parse_args()
    
    out_lib = args.output
    decoy = args.decoy
    
    param = pDeepParameter()
    
    # copyfile('tmp/data/library/empty'+os.path.splitext(out_lib)[-1], out_lib)
    
    _lib = GetLibraryWriter(out_lib, param)
    _lib.Open(out_lib)
    _lib.CreateTables()
    _lib.decoy = decoy
    
    seqlib = SequenceLibrary(min_charge = args.min_precursor_charge, max_charge = args.max_precursor_charge, min_peptide_len=args.min_peptide_length, max_peptide_len=args.max_peptide_length, min_precursor_mz = args.min_precursor_mz, max_precursor_mz = args.max_precursor_mz, varmod=args.varmod, fixmod=args.fixmod)
    
    varmod_set = set(args.varmod.strip(',').split(",")) if args.varmod else set()
    fixmod_set = set(args.fixmod.strip(',').split(",")) if args.fixmod else set()
    def _add_mod_from_library(peptide_list):
        for seq,mod,charge in peptide_list:
            if mod:
                mods = mod.strip(";").split(";")
                for mod in mods:
                    modname = mod[mod.find(',')+1:]
                    if modname not in varmod_set and modname not in fixmod_set: 
                        fixmod_set.add(modname)
                        args.fixmod += ','+modname
        
    if args.input.endswith('.peplib') or args.input.endswith('.peplib.txt') :
        peptide_list, pep_pro_dict = seqlib.PeptideListFromPeptideFile(args.input)
    elif args.input.endswith('.modseq') or args.input.endswith('.modseq.txt'):
        peptide_list, pep_pro_dict = ReadModSeq(args.input)
        _add_mod_from_library(peptide_list)
    elif args.input.endswith('.fasta'):
        peptide_list, pep_pro_dict = _from_fasta(seqlib, args.input, args.proteins)
    elif args.input.endswith('.tsv') or args.input.endswith('.csv'):
        peptide_list, pep_pro_dict = _from_tsv(args.input)
        _add_mod_from_library(peptide_list)
    elif args.input.endswith('.pqp') or args.input.endswith('.osw'):
        peptide_list, pep_pro_dict = _from_pqp(args.input)
        _add_mod_from_library(peptide_list)
    elif args.input.endswith('.dlib') or args.input.endswith('.elib'):
        peptide_list, pep_pro_dict = _from_dlib(args.input)
        _add_mod_from_library(peptide_list)
    elif args.input.endswith('evidence.txt'):
        peptide_list, pep_pro_dict = _from_MaxQuant(args.input)
        _add_mod_from_library(peptide_list)
    elif args.input.endswith('.spectra'):
        peptide_list, pep_pro_dict = _from_pFind(args.input)
        _add_mod_from_library(peptide_list)
    else:
        peptide_list, pep_pro_dict = ReadModSeq(args.input)
        _add_mod_from_library(peptide_list)
    
    if args.RT_input:
        RT_seqlib = SequenceLibrary(min_charge = args.min_precursor_charge, max_charge = args.max_precursor_charge, min_peptide_len=args.min_peptide_length, max_peptide_len=args.max_peptide_length, min_precursor_mz = args.min_precursor_mz, max_precursor_mz = args.max_precursor_mz, varmod=args.varmod, fixmod=args.fixmod)
        if args.RT_input.endswith('.peplib') or args.RT_input.endswith('.peplib.txt') :
            RT_peptide_list, RT_pep_pro_dict = RT_seqlib.PeptideListFromPeptideFile(args.RT_input)
        elif args.RT_input.endswith('.modseq') or args.RT_input.endswith('.modseq.txt'):
            RT_peptide_list, RT_pep_pro_dict = ReadModSeq(args.RT_input)
            _add_mod_from_library(RT_peptide_list)
        elif args.RT_input.endswith('.fasta'):
            RT_peptide_list, RT_pep_pro_dict = _from_fasta(RT_seqlib, args.RT_input, args.RT_proteins)
        elif args.RT_input.endswith('.tsv') or args.RT_input.endswith('.csv'):
            RT_peptide_list, RT_pep_pro_dict = _from_tsv(args.RT_input)
            _add_mod_from_library(RT_peptide_list)
        else:
            RT_peptide_list, RT_pep_pro_dict = ReadModSeq(args.RT_input)
    
    if args.spikein:
        spkin_seqlib = SequenceLibrary(min_charge = args.min_precursor_charge, max_charge = args.max_precursor_charge, min_peptide_len=args.min_peptide_length, max_peptide_len=args.max_peptide_length, min_precursor_mz = args.min_precursor_mz, max_precursor_mz = args.max_precursor_mz, varmod=args.spikein_varmod, fixmod=args.spikein_fixmod)
        if args.spikein.endswith('.modseq') or args.spikein.endswith('.modseq.txt'):
            spkin_list, spkin_dict = ReadModSeq(args.spikein)
        elif args.spikein.endswith('.peplib') or args.spikein.endswith('.peplib.txt'):
            spkin_list, spkin_dict = spkin_seqlib.PeptideListFromPeptideFile(args.spikein)
        elif args.spikein.endswith('.fasta'):
            spkin_list, spkin_dict = _from_fasta(spkin_seqlib, args.spikein, args.spikein_proteins)
        elif args.spikein.endswith('.tsv') or args.spikein.endswith('.csv'):
            spkin_list, spkin_dict = _from_tsv(args.spikein)
        elif args.input.endswith('evidence.txt'):
            spkin_list, spkin_dict = _from_MaxQuant(args.spikein)
        elif args.input.endswith('.spectra'):
            spkin_list, spkin_dict = _from_pFind(args.spikein)
        else:
            spkin_list, spkin_dict = ReadModSeq(args.spikein)
        _add_mod_from_library(spkin_list)
            
        peptide_list.extend(spkin_list)
        for seq, pro in spkin_dict.items():
            if seq not in pep_pro_dict:
                pep_pro_dict[seq] = pro
                
    print("[pDeep Info] fix modifications included: '%s'"%args.fixmod)
    print("[pDeep Info] var modifications included: '%s'"%args.varmod)
        
    if args.tune_psm and args.raw and RawFileReader:
        raw_path = args.raw
        psmRT = GeneratePSMFile(args.tune_psm, raw_path)
        psmLabel = Run_psmLabel(psmRT, raw_path)
    else:
        psmLabel = args.psmlabel
        psmRT = args.psmRT
    if psmLabel:
        Sort_psmLabel(psmLabel)
        
    param = Set_pDeepParam(param, model=args.model, instrument=args.instrument, ce=args.ce, psmLabel=psmLabel, psmRT=psmRT, fixmod=args.fixmod, varmod=args.varmod, psmLabel_test=psmLabel, n_tune=args.n_tune_psm)
        
    
    prediction = tune_and_predict.run(param, peptide_list)
    
    _lib.min_mz = args.min_mz
    _lib.max_mz = args.max_mz
    _lib.min_intensity = args.min_intensity
    _lib.least_n_peaks = args.least_n_peaks
    _lib.ion_types = args.ion_type.split(",")
    _lib.UpdateByPrediction(prediction, pep_pro_dict)
    _lib.Close()
    
    if args.RT_input:
        RT_lib = GetLibraryWriter(args.RT_tsv, param)
        RT_lib.Open(args.RT_tsv)
        RT_lib.decoy = None
        RT_prediction = tune_and_predict.run_predict(param, RT_peptide_list)
        RT_lib.min_mz = args.min_mz
        RT_lib.max_mz = args.max_mz
        RT_lib.min_intensity = args.min_intensity
        RT_lib.least_n_peaks = args.least_n_peaks
        RT_lib.ion_types = args.ion_type.split(",")
        RT_lib.UpdateByPrediction(RT_prediction, RT_pep_pro_dict)
        RT_lib.Close()
