import sys
import os
from ..spectral_library.library_base import SequenceLibrary 
from ..sequence.protein_infer import infer_protein
from ..data_generator import *
from .generate_predicted_speclib import _from_fasta


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generating modified sequences from fasta.')
    
    parser.add_argument('--input', type=str, required=True, help='.fasta for generating peptides from fasta')
    parser.add_argument('--output', type=str, required=True, help='.modseq.txt')
    
    parser.add_argument('--varmod', type=str, default="Oxidation[M]", required=False, help='Variable modifications, seperated by ",".')
    parser.add_argument('--fixmod', type=str, default="Carbamidomethyl[C]", required=False, help='Fixed modifications, seperated by ",".')
    
    parser.add_argument('--min_precursor_charge', type=int, default=2, required=False, help='Min precursor charge of peptides in the library')
    parser.add_argument('--max_precursor_charge', type=int, default=4, required=False, help='Max precursor charge of peptides in the library')
    parser.add_argument('--min_precursor_mz', type=float, default=400, required=False, help='Min precursor mz of peptides in the library')
    parser.add_argument('--max_precursor_mz', type=float, default=1200, required=False, help='Max precursor mz of peptides in the library')
    parser.add_argument('--min_peptide_length', type=int, default=6, required=False, help='Min peptide length of peptides in the library')
    parser.add_argument('--max_peptide_length', type=int, default=60, required=False, help='Max peptide length of peptides in the library')
    parser.add_argument('--min_varmod', type=int, default=0, required=False, help='Min variable modifications of peptides in the library')
    parser.add_argument('--max_varmod', type=int, default=1, required=False, help='Max variable modifications of peptides in the library')
    
    args = parser.parse_args()
    seqlib = SequenceLibrary(min_charge = args.min_precursor_charge, max_charge = args.max_precursor_charge, min_precursor_mz = args.min_precursor_mz, min_peptide_len=args.min_peptide_length, max_peptide_len=args.max_peptide_length, max_precursor_mz = args.max_precursor_mz, varmod=args.varmod, fixmod=args.fixmod, min_varmod=args.min_varmod, max_varmod=args.max_varmod)
    
    with open(args.output,"w") as f:
        f.write("peptide\tmodinfo\tcharge\tprotein\n")
        peptide_list, pep_pro_dict = _from_fasta(seqlib, args.input)
        for seq, mod, charge in peptide_list:
            pro = pep_pro_dict[seq]
            idx = mod.find('[')
            while idx != -1:
                mod = mod[:idx+1] + mod[idx+1].upper() + mod[idx+2:]
                idx = mod.find('[', idx+2)
            f.write("\t".join([seq.upper(), mod, str(charge), pro]) + "\n")
    