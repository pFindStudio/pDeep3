import numpy as np

from ..config.modification import mod_dict
from .modloss_priority import priority as PTM_NL_priority


ASCII = 128

class BaseMass:
    def __init__(self):
        self.mass_H = 1.0078250321
        self.mass_O = 15.9949146221
        self.mass_N = 14.0030740052
        self.mass_C = 12.00
        self.mass_isotope = 1.003
        self.mass_proton = 1.007276

        self.mass_H2O = self.mass_H * 2 + self.mass_O
        self.mass_CO = self.mass_C + self.mass_O
        self.mass_CO2 = self.mass_C + self.mass_O * 2
        self.mass_NH = self.mass_N + self.mass_H
        self.mass_NH3 = self.mass_N + self.mass_H * 3
        self.mass_HO = self.mass_H + self.mass_O

class PeptideIonCalculator:

    def __init__(self):
        self.base_mass = BaseMass()

        self.AAMass = np.zeros(ASCII, dtype=np.float64)
        self.AAMass[ord('A')] = 71.037114
        self.AAMass[ord('C')] = 103.009185
        self.AAMass[ord('D')] = 115.026943
        self.AAMass[ord('E')] = 129.042593
        self.AAMass[ord('F')] = 147.068414
        self.AAMass[ord('G')] = 57.021464
        self.AAMass[ord('H')] = 137.058912
        self.AAMass[ord('I')] = 113.084064
        self.AAMass[ord('J')] = 114.042927
        self.AAMass[ord('K')] = 128.094963
        self.AAMass[ord('L')] = 113.084064
        self.AAMass[ord('M')] = 131.040485
        self.AAMass[ord('N')] = 114.042927
        self.AAMass[ord('P')] = 97.052764
        self.AAMass[ord('Q')] = 128.058578
        self.AAMass[ord('R')] = 156.101111
        self.AAMass[ord('S')] = 87.032028
        self.AAMass[ord('T')] = 101.047679
        self.AAMass[ord('V')] = 99.068414
        self.AAMass[ord('W')] = 186.079313
        self.AAMass[ord('Y')] = 163.06332

        self.ModMass = {}
        for modname, modinfo in mod_dict.items():
            modinfo = modinfo.split(' ')
            modmass = float(modinfo[2])
            mod_neutral_loss = 0
            if modinfo[4] != '0':
                mod_neutral_loss = float(modinfo[5])
            self.ModMass[modname] = (modmass, mod_neutral_loss)

    def get_aamass(self, aa):
        return self.AAMass[ord(aa)]

    def set_aamass(self, aa, mass):
        self.AAMass[ord(aa)] = mass

    def calc_aamass_cumsum(self, peptide):
        return np.cumsum(self.AAMass[[ord(aa) for aa in peptide]])
        
    def calc_modification_mass(self, peptide, modinfo):
        modmass = np.zeros(len(peptide) + 2)
        if modinfo == "":
            return modmass, None, None
            
        lossmass = np.zeros(len(peptide) + 2)
        retmods = [""]*(len(peptide)+1)

        items = modinfo.strip(";").split(";")
        modlist = []

        for mod in items:
            strSite, modname = mod.split(",")
            site = int(strSite)
            modlist.append((site, modname))

        for site, modname in modlist:
            _mono, _loss = self.ModMass[modname]
            modmass[site] = _mono
            lossmass[site] = _loss
            if modname in PTM_NL_priority:
                retmods[site] = modname
            else:
                retmods[site] = ""  # 0 PTM_NL_priority
        return np.cumsum(modmass), lossmass, retmods

    def calc_b_ions_and_pepmass(self, peptide, mod_cumsum):
        aa_cumsum = self.calc_aamass_cumsum(peptide)
        b_ions = aa_cumsum[:-1] + mod_cumsum[1:-2]
        pepmass = aa_cumsum[-1] + mod_cumsum[-1] + self.base_mass.mass_H2O
        return b_ions, pepmass

    def calc_y_from_b(self, bions, pepmass):
        return pepmass - bions
        
    def calc_by_and_pepmass(self, peptide, modinfo, max_charge = 2):
        mod_cumsum, _, _ = self.calc_modification_mass(peptide, modinfo)
        bions, pepmass = self.calc_b_ions_and_pepmass(peptide, mod_cumsum)
        yions = self.calc_y_from_b(bions, pepmass)
        
        ions = [bions/charge+self.base_mass.mass_proton for charge in range(1, max_charge+1)]
        ions.extend([yions/charge+self.base_mass.mass_proton for charge in range(1, max_charge+1)])
        
        return np.array(ions).T, pepmass #.T, from 4xn to nx4, time step first

    def calc_a_from_b(self, bions):
        return bions - self.base_mass.mass_CO

    def calc_c_from_b(self, bions):
        return bions + self.base_mass.mass_NH3

    def calc_z_from_b(self, bions, pepmass):
        return (pepmass - self.base_mass.mass_NH3 + self.base_mass.mass_H) - bions

    def calc_H2O_loss(self, ions):
        return ions - self.base_mass.mass_H2O

    def calc_NH3_loss(self, ions):
        return ions - self.base_mass.mass_NH3

    def calc_Nterm_modloss(self, ions, modloss_list, modname_list):
        loss_nterm = modloss_list[0]
        prev_prior = PTM_NL_priority[modname_list[0]] if modname_list[0] in PTM_NL_priority else 0

        if not modloss_list or np.all(modloss_list == 0): return None

        ret = np.zeros(len(ions), dtype=float)
        for i in range(len(ions)):
            if modloss_list[i + 1] != 0:
                if modname_list[i + 1] in PTM_NL_priority:
                    if PTM_NL_priority[modname_list[i + 1]] >= prev_prior:
                        loss_nterm = modloss_list[i + 1]
                        prev_prior = PTM_NL_priority[modname_list[i + 1]]
                elif prev_prior == 0:
                    loss_nterm = modloss_list[i + 1]

            if loss_nterm != 0:
                ret[i] = ions[i] - loss_nterm
            else:
                ret[i] = 0.
        return ret

    def calc_Cterm_modloss(self, ions, modloss_list, modname_list):
        last_idx = len(modloss_list) - 1
        loss_cterm = modloss_list[last_idx]
        prev_prior = PTM_NL_priority[modname_list[last_idx]] if modname_list[last_idx] in PTM_NL_priority else 0

        if not modloss_list or np.all(modloss_list == 0): return None

        ret = np.zeros(len(ions), dtype=float)
        for i in range(len(ions) - 1, -1, -1):
            if modloss_list[i + 2] != 0:
                if modname_list[i + 2] in PTM_NL_priority:
                    if PTM_NL_priority[modname_list[i + 2]] >= prev_prior:
                        loss_cterm = modloss_list[i + 2]
                        prev_prior = PTM_NL_priority[modname_list[i + 2]]
                elif prev_prior == 0:
                    loss_cterm = modloss_list[i + 2]

            if loss_cterm != 0:
                ret[i] = ions[i] - loss_cterm
            else:
                ret[i] = 0.
        return ret 