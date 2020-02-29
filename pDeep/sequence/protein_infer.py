from .ahocorasick import Trie
from .fasta_reader import read_all_proteins

def infer_protein(seq_list, protein_dict):
    print("Inferring proteins ... ")

    ret_pro_acs = []
    ret_start_pos = []
    trie = Trie()
    pep_pro_dict = {}
    for peptide in seq_list:
        pep_pro_dict[peptide] = []
        trie.add_word(peptide, peptide)

    trie.make_automaton()

    for pro_ac, protein in protein_dict.items():
        for (pos, peptides) in trie.iter(protein.seq):
            for peptide in peptides:
                pep_pro_dict[peptide].append((pro_ac, pos - len(peptide) + 1))

    print("End inferring")

    return pep_pro_dict

def infer_protein_fasta(seq_list, fasta):
    protein_dict = read_all_proteins(fasta)
    return infer_protein(seq_list, protein_dict), protein_dict