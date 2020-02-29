from ..utils.mass_calc import PeptideIonCalculator
from ..peptide.peptide import get_peptidoforms_from_fasta
from ..peptide.digest import DigestConfig
from ..peptide.protein_infer import *

class SeqLibrary(object):
    def __init__(self, min_charge = 2, max_charge = 4, 
                       min_precursor_mz = 400, max_precursor_mz = 1000, 
                       varmods = "Oxidation[M]", fixmods = "Carbamidomethyl[C]", 
                       min_varmod = 0, max_varmod = 1):
        self.ion_calc = PeptideIonCalculator()
        self.digest_config = DigestConfig()
        self.min_charge = min_charge
        self.max_charge = max_charge
        self.min_precursor_mz = min_precursor_mz
        self.max_precursor_mz = max_precursor_mz
        self.varmods = varmods
        self.fixmods = fixmods
        self.min_varmod = min_varmod
        self.max_varmod = max_varmod
        
        self.peptide_list = []
        self.peptide_to_protein_dict = {}
        
    def PeptideListFromFasta(self, fasta):
        self.peptide_list = []
        self.peptide_to_protein = {}
        
        modseq_list, protein_dict = get_peptidoforms_from_fasta(fasta, self.digest_config, self.varmods, self.fixmods, self.min_varmod, self.max_varmod)
        
        seq_for_proinfer = []
        for seq, modinfo in modseq_list:
            pepmass = self.ion_calc.calc_pepmass(seq, modinfo)
            for charge in range(self.min_charge, self.max_charge+1):
                mz = pepmass/charge + self.ion_calc.base_mass.mass_proton
                if mz >= self.min_precursor_mz and mz <= self.max_precursor_mz:
                    self.peptide_list.append((seq, modinfo, charge))
                    if seq_for_proinfer and seq_for_proinfer[-1] != seq:
                        seq_for_proinfer.append(seq)
                        
        pep_pro_dict = infer_protein(seq_for_proinfer, protein_dict)
        self.peptide_to_protein_dict = dict(zip([peptide, "/".join([pro_ac for pro_ac, site in pros.ites]) for peptide prosites in pep_pro_dict.items()]))
        return self.peptide_list, self.peptide_to_protein_dict
        
                