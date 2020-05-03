from .model_tf import use_tf2
from .config.pDeep_config import HCD_CommonMod_Config

class pDeepParameter:
    def __init__(self, cfg = None):
        self._max_ion_charge = 2 # ion charge = 1 or 2 for pDeep model
        self._ion_terms = {'a{}': 'n', 'x{}': 'c', 'b{}': 'n', 'y{}': 'c', 'c{}': 'n', 'z{}': 'c', 'b{}-ModLoss': 'n', 'y{}-ModLoss': 'c', 'b{}-H2O': 'n', 'y{}-H2O': 'c', 'b{}-NH3': 'n', 'y{}-NH3': 'c'}
        for iontype, term in list(self._ion_terms.items()):
            self._ion_terms[iontype.format("")] = term
        for iontype, term in list(self._ion_terms.items()):
            self._ion_terms[iontype.lower()] = term
            
        ######################################################################
        
        self.config = None
        
        self.model = "HCD"
        
        self.library_ion_types = self.ion_types
    
        self.RT_model = "tmp/model/RT-model.ckpt"
        
        self._pDeepModel = None
        self._pDeepRT = None
        
        self.threads = 4
        self.predict_batch = 4096
        
        self.fixmod = ["Carbamidomethyl[C]"]
        self.varmod = "Acetyl[ProteinN-term],Oxidation[M],Phospho[Y],Phospho[S],Phospho[T]".split(",")
        self.min_varmod = 0
        self.max_varmod = 3

        self.predict_input = "none"
        self.predict_instrument = "QE"
        self.predict_nce = 27
        self.predict_output = "none"

        # tune parameters:
        self.tune_psmlabels = []
        self.tune_instruments = []
        self.tune_nces = []
        self.tune_RT_psmlabel = ""
        self.epochs = 2
        self.n_tune_per_psmlabel = 1000
        self.tune_batch = 1024
        self.dropout = 0
        self.tune_save_as = None
        self.tune_RT_save_as = None
        
        self.test_psmlabels = []
        self.test_instruments = []
        self.test_nces = []
        self.test_RT_psmlabel = ""
        self.n_test_per_psmlabel = 100000000

        self.fasta = None

        if cfg: self._read_cfg(cfg)
        
    def InitConfig(self, config = None):
        if config is None:
            self.config = HCD_CommonMod_Config()
        else:
            self.config = config
        self.config.SetFixMod(self.fixmod)
        self.config.SetVarMod(self.varmod)
        self.config.SetIonTypes(self.ion_types)
        self.config.min_var_mod_num = self.min_varmod
        self.config.max_var_mod_num = self.max_varmod
        
    @property
    def tune_save_as(self):
        return self._tune_save
    @tune_save_as.setter
    def tune_save_as(self, save_as):
        if save_as is None: self._tune_save = None
        elif save_as.endswith(".ckpt"):
            self._tune_save = save_as
        else:
            self._tune_save = save_as + ".ckpt"
        
    @property
    def tune_RT_save_as(self):
        return self._tune_RT_save
    @tune_RT_save_as.setter
    def tune_RT_save_as(self, save_as):
        if save_as is None: self._tune_RT_save = None
        elif save_as.endswith(".ckpt"):
            self._tune_RT_save = save_as
        else:
            self._tune_RT_save = save_as + ".ckpt"
        
    @property
    def model(self):
        return self._model
    @model.setter
    def model(self, _model):
        if _model.upper() == "ETHCD":
            if use_tf2: self._model = "tmp/model/EThcD-tf2-bycz-Lumos-CE28.ckpt"
            else: self._model = "tmp/model/EThcD-bycz-Lumos-CE28.ckpt"
            self.ion_types = ['b{}', 'y{}', 'c{}', 'z{}']
        elif _model.upper() == "HCD":
            if use_tf2: self._model = "tmp/model/pretrain-tf2-200412.ckpt"
            else: self._model = "tmp/model/pretrain-180921-modloss-mod8D.ckpt"
            self.ion_types = ['b{}', 'y{}', 'b{}-ModLoss', 'y{}-ModLoss']
        elif _model.upper() == "PHOSPHO" or _model.upper() == "PHOS":
            if use_tf2: self._model = "tmp/model/pretrain-tf2-200412.ckpt"
            else: self._model = "tmp/model/pretrain-180921-modloss-mod8D.ckpt"
            self.ion_types = ['b{}', 'y{}', 'b{}-ModLoss', 'y{}-ModLoss']
        else:
            self._model = _model
            self.ion_types = ['b{}', 'y{}', 'b{}-ModLoss', 'y{}-ModLoss']
            
    @property
    def ion_types(self):
        return self._ion_types
    @ion_types.setter
    def ion_types(self, _ion_types):
        self._ion_types = _ion_types
        if self.config is not None: self.config.SetIonTypes(_ion_types)
        self._ion_type_idx = {}
        for i in range(len(self._ion_types)):
            self._ion_type_idx[self._ion_types[i]] = i
            self._ion_type_idx[self._ion_types[i].format('')] = i

    def GetPredictedIonTypeIndices(self, ion_types):
        return [self._ion_type_idx[iontype]*self._max_ion_charge+ch for iontype in ion_types for ch in range(self._max_ion_charge)]
    
    def _read_cfg(self, cfg):
        with open(cfg) as f:
            lines = f.readlines()

            # def parse_folder(line):
                # items = line.split("=")[1].split("|")
                # items[0] = items[0].strip()  # folder
                # items[1] = items[1].strip()  # instrument
                # items[2] = int(items[2].strip())  # NCE
                # if len(items) >= 4: items[3] = items[3].strip()
                # return items

            def get_int(line):
                return int(get_str(line))

            def get_str(line):
                return line[line.find("=")+1:].strip()

            for i in range(len(lines)):
                line = lines[i].strip()
                if line.startswith("model"):
                    self.model = get_str(line)
                elif line.startswith("RT_model"):
                    self.RT_model = get_str(line)
                elif line.startswith("threads"):
                    self.threads = get_int(line)
                elif line.startswith("predict_batch"):
                    self.predict_batch = get_int(line)
                    
                elif line.startswith("mod_no_check"):
                    self.fixmod = get_str(line).strip(",").split(",")
                elif line.startswith("mod_check"):
                    self.varmod = get_str(line).strip(",").split(",")
                elif line.startswith("min_mod_check"):
                    self.min_varmod = get_int(line)
                elif line.startswith("max_mod_check"):
                    self.max_varmod = get_int(line)

                # elif line.startswith("predict_input"):
                    # self.predict_input, self.predict_instrument, self.predict_nce = parse_folder(line)[:3]
                    # self.predict_nce = int(self.predict_nce)
                elif line.startswith("predict_instrument"):
                    self.predict_instrument = get_str(line)
                elif line.startswith("predict_instrument"):
                    self.predict_nce = float(get_str(line))
                elif line.startswith("predict_input"):
                    self.predict_input = get_str(line)
                elif line.startswith("predict_output"):
                    self.predict_output = get_str(line)
                elif line.startswith("fasta"):
                    self.fasta = get_str(line)

                # tune_parameters
                elif line.startswith("tune_psmlabels"):
                    line = get_str(line)
                    if line: self.tune_psmlabels = [s.strip() for s in line.split("|")]
                elif line.startswith("tune_instruments"):
                    line = get_str(line)
                    if line: self.tune_instruments = [s.strip() for s in line.split("|")]
                elif line.startswith("tune_nces"):
                    line = get_str(line)
                    if line: self.tune_nces = [float(s.strip()) for s in line.split("|")]
                elif line.startswith("tune_epochs"):
                    self.epochs = get_int(line)
                elif line.startswith("tune_batch"):
                    self.train_batch = get_int(line)
                elif line.startswith("n_tune_per_psmlabel"):
                    self.n_tune_per_psmlabel = get_int(line)
                elif line.startswith("tune_RT_psmlabel"):
                    self.tune_RT_psmlabel = get_str(line)
                elif line.startswith("tune_save_as"):
                    self.tune_save_as = get_str(line)
                elif line.startswith("tune_RT_save_as"):
                    self.tune_RT_save_as = get_str(line)
                    
                elif line.startswith("test_psmlabels"):
                    line = get_str(line)
                    if line: self.test_psmlabels = [s.strip() for s in line.split("|")]
                elif line.startswith("test_instruments"):
                    line = get_str(line)
                    if line: self.test_instruments = [s.strip() for s in line.split("|")]
                elif line.startswith("test_nces"):
                    line = get_str(line)
                    if line: self.test_nces = [float(s.strip()) for s in line.split("|")]
                elif line.startswith("test_RT_psmlabel"):
                    self.test_RT_psmlabel = get_str(line)
                elif line.startswith("n_test_per_psmlabel"):
                    self.n_test_per_psmlabel = get_int(line)
