from .digest import digest

anyNterm = "AnyN-term"
max_peptidoforms_per_seq = 100

class Protein:
    def __init__(self, AC, DE, seq):
        self.AC = AC
        self.DE = DE
        self.seq = seq

def read_proteins_by_AC_list(fasta, protein_list):
    '''
    For exmaple, 'sp|P08603|CFAH_HUMAN', protein in protein_list could be 'P08603' or 'CFAH_HUMAN'.
    '''
    def check_AC_in_protein_list(ac):
        for pro in protein_list:
            if pro in ac: return True
        return False
        
    ret = {}
    with open(fasta) as f:
        seq = None
        while True:
            line = f.readline()
            if line == "":
                if seq is not None:
                    ret[ac] = Protein(ac, de, seq)
                break
            elif line.startswith(">"):
                if seq is not None:
                    ret[ac] = Protein(ac, de, seq)
                line = line.strip()[1:].replace("\t", " ")
                ac, de = line.split(" ", 1)
                if check_AC_in_protein_list(ac):
                    seq = ""
                else:
                    seq = None
            else:
                if seq is not None:
                    seq += line.strip()
    return ret


def read_all_proteins(fasta):
    ret = {}
    with open(fasta) as f:
        seq = ""
        while True:
            line = f.readline()
            if line == "":
                if len(seq) > 0:
                    ret[ac] = Protein(ac, de, seq)
                break
            elif line.startswith(">"):
                if len(seq) > 0:
                    ret[ac] = Protein(ac, de, seq)
                line = line.strip()[1:].replace("\t", " ")
                ac, de = line.split(" ", 1)
                seq = ""
            else:
                seq += line.strip()
    return ret

def replace_Nterm(siteaa):
    if "N-term" in siteaa:
        return anyNterm
    else:
        return siteaa
    
def gen_mod_dict(modlist):
    mod_dict = {}

    def get_mod_AA(mod):
        return mod[mod.rfind('[') + 1:mod.rfind(']')]

    for mod in modlist:
        siteaa = get_mod_AA(mod)
        siteaa = replace_Nterm(siteaa)
        if siteaa not in mod_dict:
            mod_dict[siteaa] = [mod]
        else:
            mod_dict[siteaa].append(mod)
    return mod_dict

def generate_mod_dict(modstr):
    '''
    modstr: Oxidation[M],Deamidated[N] or Oxidation[M];Deamidated[N]
    '''
    if ',' in modstr: return gen_mod_dict(modstr.strip(',').split(','))
    else: return gen_mod_dict(modstr.strip(';').split(';'))

def modstr(mod, site):
    return str(site) + "," + mod + ";"


def get_fix_mod(seq, fixmod_dict):
    mod = ""
    s = 0
    if anyNterm in fixmod_dict:
        mod += modstr(fixmod_dict[anyNterm][0], 0)
        s = 1
    for i in range(s, len(seq)):
        if seq[i] in fixmod_dict:
            mod += modstr(fixmod_dict[seq[i]][0], i + 1)
    return mod, s


def add_modifications(peptide, varmod_dict, fixmod_dict, min_var_mod=0, max_var_mod=1):
    modseq_list = []
    mod, s = get_fix_mod(peptide, fixmod_dict)

    def add_mod_recur(modseq_list, peptide, i, mod, varmod_dict, min_var_mod, max_var_mod, n_var_mod):
        if len(modseq_list) >= max_peptidoforms_per_seq:
            return
        elif i > len(peptide):
            if n_var_mod >= min_var_mod: modseq_list.append((peptide, mod)) #last protein
        else:
            add_mod_recur(modseq_list, peptide, i + 1, mod, varmod_dict, min_var_mod, max_var_mod, n_var_mod)
            if peptide[i - 1] in varmod_dict and n_var_mod < max_var_mod:
                for modname in varmod_dict[peptide[i - 1]]:
                    add_mod_recur(modseq_list, peptide, i + 1, mod + modstr(modname, i), varmod_dict, min_var_mod,
                                  max_var_mod, n_var_mod + 1)

    if s == 0:
        add_mod_recur(modseq_list, peptide, 1, mod, varmod_dict, min_var_mod, max_var_mod, 0)
        if anyNterm in varmod_dict and max_var_mod > 0:
            for modname in varmod_dict[anyNterm]:
                add_mod_recur(modseq_list, peptide, 2, mod + modstr(modname, 0), varmod_dict, min_var_mod, max_var_mod, 1)
    else:
        add_mod_recur(modseq_list, peptide, 2, mod, varmod_dict, min_var_mod, max_var_mod, 0)
    return modseq_list


def get_peptidoforms(pep_set, varmod_dict, fixmod_dict, min_var_mod=0, max_var_mod=1):
    modseq_list = []
    for pep in pep_set:
        modseq_list.extend(add_modifications(pep, varmod_dict, fixmod_dict, min_var_mod, max_var_mod))
    return modseq_list
    
def get_peptidoforms_from_fasta(fasta, digest_config, varmods, fixmods, min_var_mod=0, max_var_mod=1, protein_AC_list = None):
    varmod_dict = generate_mod_dict(varmods)
    fixmod_dict = generate_mod_dict(fixmods)
    if not protein_AC_list: protein_dict = read_all_proteins(fasta)
    else: protein_dict = read_proteins_by_AC_list(fasta, protein_AC_list)
    peptide_set = set()
    for ac, protein in protein_dict.items():
        peptide_set = digest(protein, peptide_set, digest_config)
    return get_peptidoforms(peptide_set, varmod_dict, fixmod_dict, min_var_mod, max_var_mod), protein_dict

if __name__ == "__main__":
    fixmod_dict = gen_mod_dict(["fix[C]"])
    varmod_dict = gen_mod_dict(["a[A]", "nt[ProteinN-term]", "pS[S]", "nt[AnyN-term]"])
    l = get_peptidoforms(set(["ABCDEFSSABK", "ABCDEFSSABCK"]), varmod_dict, fixmod_dict)
    print(l)
