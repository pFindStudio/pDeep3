import numpy as np

from ..utils.mass_calc import PeptideIonCalculator
from ..sequence.peptide import get_peptidoforms_from_fasta, get_peptidoforms_from_pep2pro_dict
from ..sequence.digest import DigestConfig
from ..sequence.protein_infer import *

from ..config.modification import mod_dict

mod_mass_dict = {}
for modname, item in mod_dict.items():
    mod_mass_dict[modname] = float(item.split(" ")[2])
    
class LibraryBase(object):
    def __init__(self, pDeepParam = None):
        self._decoy = "reverse" # or pseudo_reverse or no_decoy
        self.decoy_tag = "DECOY_"
        self._ion_calc = PeptideIonCalculator()
        if pDeepParam:
            self.pDeepParam = pDeepParam
            self.ion_types = pDeepParam.library_ion_types
            # print("self.intensity_indices", self.intensity_indices)
            self.ion_terms = pDeepParam._ion_terms
            self.max_ion_charge = pDeepParam._max_ion_charge
        self.min_mz = 300
        self.max_mz = 2000
        self.min_intensity = 0.1
        self.least_n_peaks = 6
        self.inten_base = 10000
        
    @property
    def ion_types(self):
        return self._ion_types
    @ion_types.setter
    def ion_types(self, _ion_types):
        self._ion_types = _ion_types
        self.intensity_indices = self.pDeepParam.GetPredictedIonTypeIndices(self._ion_types)
        
    @property
    def decoy(self):
        return self._decoy
    @decoy.setter
    def decoy(self, _decoy):
        if _decoy == "no_decoy": self._decoy = None
        else: self._decoy = _decoy
        
    def Open(self):
        pass
    def CreateTables(self):
        pass
    def Close(self):
        pass
        
    def _get_decoy_peptide(self, seq, mod):
        if self.decoy == 'reverse':
            seq = seq[::-1]
            if mod:
                mods = mod.strip(";").split(";")
                modlist = []
                for onemod in mods:
                    site, modname = onemod.split(",")
                    site = int(site)
                    if site <= len(seq) and site != 0:
                        site = len(seq)-site+1
                    modlist.append((site, modname))
                modlist.sort(key=lambda x: x[0])
                mod = ";".join(['%d,%s'%(site, modname) for site, modname in modlist])
        else:
            seq = seq[:-1][::-1]+seq[-1]
            if mod:
                mods = mod.strip(";").split(";")
                modlist = []
                for onemod in mods:
                    site, modname = onemod.split(",")
                    site = int(site)
                    if site < len(seq) and site != 0:
                        site = len(seq)-site
                    modlist.append((site, modname))
                modlist.sort(key=lambda x: x[0])
                mod = ";".join(['%d,%s'%(site, modname) for site, modname in modlist])
        return seq, mod
    
    def _calc_ions(self, seq, mod, charge, intensities):
        pepmass, masses = self._ion_calc.calc_pepmass_and_ions_from_iontypes(seq, mod, self.ion_types, self.max_ion_charge)
        intens = intensities[:,self.intensity_indices]
            
        if self.decoy: 
            decoy_seq, decoy_mod = self._get_decoy_peptide(seq, mod)
            _, decoy_masses = self._ion_calc.calc_pepmass_and_ions_from_iontypes(decoy_seq, decoy_mod, self.ion_types, self.max_ion_charge)
        else:
            decoy_seq, decoy_mod, decoy_masses = None, None, None
            
        pepmass = pepmass / charge + self._ion_calc.base_mass.mass_proton
        # print(seq, mod, masses)
        
        _ion_sites = []
        _ion_types = []
        b_sites = np.tile(np.arange(1, len(seq)).reshape(-1,1), [1,self.max_ion_charge])
        for iontype in self.ion_types:
            iontype = iontype.format('')
            if self.ion_terms[iontype] == 'n':
                _ion_sites.append(b_sites)
            else:
                _ion_sites.append(len(seq)-b_sites)
            _ion_types.append(np.array([iontype]*((len(seq)-1)*self.max_ion_charge)).reshape(-1, self.max_ion_charge))
        sites = np.concatenate(_ion_sites, axis=1)
        types = np.concatenate(_ion_types, axis=1)
        
        charges = np.concatenate([np.full(len(seq)-1, i, dtype=int).reshape(-1,1) for i in range(1, self.max_ion_charge+1)], axis=1)
        charges = np.concatenate([charges]*len(self.ion_types), axis=1)
        
        masses = masses.reshape(-1)
        intens = intens.reshape(-1)
        sites = sites.reshape(-1)
        types = types.reshape(-1)
        charges = charges.reshape(-1)
        if self.decoy: decoy_masses = decoy_masses.reshape(-1)
        
        intens[np.abs(masses - pepmass) < 10] = 0 #delete ions around precursor m/z
        intens = intens/np.max(intens)
        
        masses = masses[intens > 0]
        sites = sites[intens > 0]
        types = types[intens > 0]
        charges = charges[intens > 0]
        if self.decoy: decoy_masses = decoy_masses[intens > 0]
        intens = intens[intens > 0]
        
            
        intens = intens[np.logical_and(masses<=self.max_mz, masses>=self.min_mz)]
        sites = sites[np.logical_and(masses<=self.max_mz, masses>=self.min_mz)]
        types = types[np.logical_and(masses<=self.max_mz, masses>=self.min_mz)]
        charges = charges[np.logical_and(masses<=self.max_mz, masses>=self.min_mz)]
        if self.decoy: decoy_masses = decoy_masses[np.logical_and(masses<=self.max_mz, masses>=self.min_mz)]
        masses = masses[np.logical_and(masses<=self.max_mz, masses>=self.min_mz)]
        
        if len(masses[intens > self.min_intensity]) >= self.least_n_peaks:
            masses = masses[intens > self.min_intensity]
            sites = sites[intens > self.min_intensity]
            types = types[intens > self.min_intensity]
            charges = charges[intens > self.min_intensity]
            if self.decoy: decoy_masses = decoy_masses[intens > self.min_intensity]
            intens = intens[intens > self.min_intensity] * self.inten_base
        else:
            indices = np.argsort(intens)[::-1]
            masses = masses[indices[:self.least_n_peaks]]
            sites = sites[indices[:self.least_n_peaks]]
            types = types[indices[:self.least_n_peaks]]
            charges = charges[indices[:self.least_n_peaks]]
            if self.decoy: decoy_masses = decoy_masses[indices[:self.least_n_peaks]]
            intens = intens[indices[:self.least_n_peaks]] * self.inten_base
        if self.decoy:
            u_masses = decoy_masses[decoy_masses < 10]
            decoy_masses[decoy_masses < 10] = np.random.random_sample(u_masses.shape)*(self.max_mz - self.min_mz) + self.min_mz
        return pepmass, masses, intens, sites, types, charges, decoy_seq, decoy_mod, decoy_masses
    
    def GetAllPeptides(self):
        pass
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}):
        pass

class SequenceLibrary(object):
    def __init__(self, min_charge = 2, max_charge = 4, 
                       min_peptide_len = 6, max_peptide_len = 60,
                       min_precursor_mz = 400, max_precursor_mz = 1200, 
                       varmod = "Oxidation[M]", fixmod = "Carbamidomethyl[C]", 
                       min_varmod = 0, max_varmod = 1):
        self.ion_calc = PeptideIonCalculator()
        self.digest_config = DigestConfig()
        self.min_charge = min_charge
        self.max_charge = max_charge
        self.min_peptide_len = min_peptide_len
        self.max_peptide_len = max_peptide_len
        self.min_precursor_mz = min_precursor_mz
        self.max_precursor_mz = max_precursor_mz
        self.varmod = varmod
        if fixmod: self.fixmod = fixmod
        else: self.fixmod = "Carbamidomethyl[C]"
        self.min_varmod = min_varmod
        self.max_varmod = max_varmod
        
        self.peptide_list = []
        self.protein_dict = {}
        
    def _add_charge(self, modseq_list):
        fmt = "Generated %d precursors (length: {} to {}, charge: {} to {}, m/z: {} to {})".format(
            self.min_peptide_len, self.max_peptide_len, self.min_charge, self.max_charge, self.min_precursor_mz, self.max_precursor_mz)
        for seq, modinfo in modseq_list:
            if len(seq) > self.max_peptide_len or len(seq) < self.min_peptide_len: continue
            pepmass = self.ion_calc.calc_pepmass(seq, modinfo)
            for charge in range(self.min_charge, self.max_charge+1):
                mz = pepmass/charge + self.ion_calc.base_mass.mass_proton
                if mz >= self.min_precursor_mz and mz <= self.max_precursor_mz:
                    self.peptide_list.append((seq, modinfo, charge))
                    if len(self.peptide_list) % 100000 == 0:
                        print(fmt%(len(self.peptide_list)))
        print(fmt%(len(self.peptide_list)))
        
    def PeptideListFromPeptideFile(self, pepfile):
        self.peptide_list = []
        self.pep2pro_dict = {}
        fmt = "Generated %d peptides"
        with open(pepfile) as f:
            head = f.readline().strip().split("\t")
            headidx = dict(zip(head, range(len(head))))
            if 'protein' in headidx: protein_idx = headidx['protein']
            else: protein_idx = None
            lines = f.readlines()
            for line in lines:
                items = line.strip().split("\t")
                if protein_idx is not None: self.pep2pro_dict[items[headidx['peptide']]] = items[protein_idx]
                else: self.pep2pro_dict[items[headidx['peptide']]] = "pDeep"
        modseq_list = get_peptidoforms_from_pep2pro_dict(self.pep2pro_dict, self.varmod, self.fixmod, self.min_varmod, self.max_varmod)
        print(fmt%(len(modseq_list)))
        
        self._add_charge(modseq_list)
        
        return self.peptide_list, self.pep2pro_dict
        
    def PeptideListFromFasta(self, fasta, protein_list = None):
        self.peptide_list = []
        self.protein_dict = {}
        fmt = "Generated %d peptides (length: {} to {})".format(self.digest_config.min_len, self.digest_config.max_len)
        modseq_list, self.protein_dict = get_peptidoforms_from_fasta(fasta, self.digest_config, self.varmod, self.fixmod, self.min_varmod, self.max_varmod, protein_list)
        print(fmt%(len(modseq_list)))
        
        self._add_charge(modseq_list)
        
        return self.peptide_list, self.protein_dict
        
if __name__ == "__main__":
    seqlib = SequenceLibrary(min_precursor_mz = 400, max_precursor_mz = 1000)
    seqlib.PeptideListFromFasta(r"C:\DataSets\fasta\SP-1503-Human_contaminant.fasta")
                