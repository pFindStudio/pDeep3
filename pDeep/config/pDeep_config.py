from .modification import get_modification


class Common_Config(object):
    def __init__(self):
        self.ion_terms = {'a{}': 'n', 'x{}': 'c', 'b{}': 'n', 'y{}': 'c', 'c{}': 'n', 'z{}': 'c', 'b{}-ModLoss': 'n', 'y{}-ModLoss': 'c', 'b{}-H2O': 'n', 'y{}-H2O': 'c', 'b{}-NH3': 'n', 'y{}-NH3': 'c'}
        for iontype, term in list(self.ion_terms.items()):
            self.ion_terms[iontype.format("")] = term
            
        self.all_mod_dict = get_modification()

        self.time_step = 100  # at most time_step+1 length peptides
        self.max_ion_charge = 2

        self.fragmentation = 'HCD'
        self.instrument_list = ['QE', 'Velos', 'Elite', 'Fusion', 'Lumos']  # they will be converted into lower case
        self.max_instrument_num = 8  # left 3 pos for the some other instruments

        self.SetFixMod(['Carbamidomethyl[C]'])
        self.varmod = []
        self.min_var_mod_num = 0
        self.max_var_mod_num = 0

        self.element_list = ['C', 'H', 'N', 'O', 'S', 'P', 'Reserve']
        # self.element_list = [elem for elem in element_list]
        self.element_list.append("unknown")
        self.element_index = dict(zip(self.element_list, range(len(self.element_list))))

        self.SetIonTypes(['b{}', 'y{}'])

    def GetModFeatureSize(self):
        return len(self.element_list)

    def GetTFOutputSize(self):
        return len(self.ion_types) * self.max_ion_charge

    def SetVarMod(self, modlist):
        self.varmod = [mod for mod in modlist]

    def SetFixMod(self, modlist):
        self.fixmod = modlist
        self.fix_aa_mod = {}
        for item in self.fixmod:
            aa = item[item.find('[') + 1:item.find(']')]
            self.fix_aa_mod[aa] = item

    def CheckFixMod(self, peptide, modlist):
        return self.CheckFixMod_fixpass(peptide, modlist)

    def CheckFixMod_fixall(self, peptide, modlist):
        if len(self.fix_aa_mod) == 0: return True
        fixed = [0] * len(peptide)
        for i in range(len(peptide)):
            if peptide[i] in self.fix_aa_mod: fixed[i] = 1
        for idx, modname in modlist:
            if idx > 0 and idx <= len(peptide):
                if peptide[idx - 1] in self.fix_aa_mod:
                    if modname == self.fix_aa_mod[peptide[idx - 1]]: fixed[idx - 1] = 0
        return sum(fixed) == 0

    def CheckFixMod_fixpass(self, peptide, modlist):
        return True

    def SetIonTypes(self, ion_types):
        self.ion_types = ion_types

        self.SetPredictIonIndex()  # predicted intensity index of different ions in ndarray

    def GetIonTypeNames(self):
        return [ion_type.format("") for ion_type in self.ion_types]

    def GetIonNameBySite(self, peptide, site, ion_type):
        if self.ion_terms[ion_type] == 'c':
            return ion_type.format(len(peptide) - site)
        else:
            return ion_type.format(site)

    def GetIntenIdx(self, iontype):
        return self.pred_ion_idx[iontype]

    def SetPredictIonIndex(self):
        self.pred_ion_idx = dict(zip(self.ion_types, range(len(self.ion_types))))
        self.pred_ion_idx.update(dict(zip([ion_type.format("") for ion_type in self.ion_types], range(len(self.ion_types)))))

    def GetIonIndexByIonType(self, ion_type, ion_charge):
        if ion_charge > self.max_ion_charge: return None
        if not ion_type in self.pred_ion_idx: return None
        return self.pred_ion_idx[ion_type] * self.max_ion_charge + ion_charge - 1

    def GetIntenFromNDArrayByLossName(self, inten_ndarray, loss_name=None):
        # loss_name = None: noloss, "ModLoss", "H2O", "NH3"
        # shape of inten_ndarray: (peptides, cleavage_sites, ion_types)
        idxes = []
        if loss_name == None: loss_name = "{}"
        for ion_type in self.ion_types:
            if ion_type.endswith(loss_name):
                for ion_charge in range(1, self.max_ion_charge + 1):
                    idxes.append(self.GetIonIndexByIonType(ion_type, ion_charge))
        return inten_ndarray[:, :, idxes]


class HCD_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()


class ETD_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.SetIonTypes(['c{}', 'z{}'])
        self.fragmentation = 'ETD'


class EThcD_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.SetIonTypes(['b{}', 'y{}', 'c{}', 'z{}'])
        self.varmod.extend(['Oxidation[M]', 'Phospho[S]', 'Phospho[T]', 'Phospho[Y]'])
        self.fragmentation = 'EThcD'


class HCD_pho_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.varmod.extend(['Oxidation[M]', 'Phospho[S]', 'Phospho[T]', 'Phospho[Y]'])
        self.min_var_mod_num = 0
        self.max_var_mod_num = 2


class HCD_CommonMod_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.SetIonTypes(['b{}', 'y{}'])
        self.varmod = ['Oxidation[M]', 'Deamidated[N]', 'Deamidated[Q]', 'Acetyl[ProteinN-term]', 'Acetyl[AnyN-term]',
                       'Formyl[AnyN-term]', 'Gln->pyro-Glu[AnyN-termQ]']
        self.min_var_mod_num = 0
        self.max_var_mod_num = 2


class HCD_AllMod_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.SetIonTypes(['b{}', 'y{}'])
        self.varmod = [mod for mod in self.all_mod_dict]
        self.min_var_mod_num = 0
        self.max_var_mod_num = 2

