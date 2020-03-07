from ..utils.mass_calc import PeptideIonCalculator
from ..sequence.peptide import get_peptidoforms_from_fasta
from ..sequence.digest import DigestConfig
from ..sequence.protein_infer import *

class LibraryBase(object):
    def Open(self):
        pass
    def Close(self):
        pass
    def GetAllPeptides(self):
        pass
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}, peak_selection = "intensity", threshold = 0.01, mass_upper = 2000):
        pass
    def Submit(self):
        pass

class SequenceLibrary(object):
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
        self.protein_dict = {}
        
    def PeptideListFromFasta(self, fasta, protein_AC_list = None):
        self.peptide_list = []
        self.protein_dict = {}
        fmt = "Generated %d peptides (length: {} to {})".format(self.digest_config.min_len, self.digest_config.max_len)
        if not protein_AC_list:
            modseq_list, self.protein_dict = get_peptidoforms_from_fasta(fasta, self.digest_config, self.varmods, self.fixmods, self.min_varmod, self.max_varmod)
        else:
            modseq_list, self.protein_dict = get_peptidoforms_from_fasta(fasta, self.digest_config, self.varmods, self.fixmods, self.min_varmod, self.max_varmod, protein_AC_list)
        print(fmt%(len(modseq_list)))
        
        fmt = "Generated %d precursors (charge: {} to {}, m/z: {} to {})".format(self.min_charge, self.max_charge, self.min_precursor_mz, self.max_precursor_mz)
        seq_for_proinfer = []
        for seq, modinfo in modseq_list:
            pepmass = self.ion_calc.calc_pepmass(seq, modinfo)
            for charge in range(self.min_charge, self.max_charge+1):
                mz = pepmass/charge + self.ion_calc.base_mass.mass_proton
                if mz >= self.min_precursor_mz and mz <= self.max_precursor_mz:
                    self.peptide_list.append((seq, modinfo, charge))
                    if not seq_for_proinfer or seq_for_proinfer[-1] != seq:
                        seq_for_proinfer.append(seq)
                    if len(self.peptide_list) % 100000 == 0:
                        print(fmt%(len(self.peptide_list)))
        print(fmt%(len(self.peptide_list)))
        
        return self.peptide_list, self.protein_dict
        
if __name__ == "__main__":
    seqlib = SequenceLibrary(min_precursor_mz = 400, max_precursor_mz = 1000)
    seqlib.PeptideListFromFasta(r"C:\DataSets\fasta\SP-1503-Human_contaminant.fasta")
                