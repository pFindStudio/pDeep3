class pDeepParameter:
    def __init__(self, cfg):
        self._ion_types = ['b{}', 'y{}', 'b{}-ModLoss', 'y{}-ModLoss']
    
        self.model = ""
        self.threads = 4
        self.predict_batch = 4096
        
        self.fixmod = []
        self.varmod = []
        self.min_varmod = 0
        self.max_varmod = 3

        self.predict_input = "peptide.txt"
        self.predict_instrument = "QE"
        self.predict_nce = 28
        self.predict_output = "predict.pdp"

        # tune parameters:
        self.tune_psmlabels = []
        self.epochs = 2
        self.n_tune_per_psmlabel = 100
        self.tune_batch = 1024
        
        self.test_psmlabels = []
        self.n_test_per_psmlabel = 100000000

        self.fasta = None

        self._read_cfg(cfg)

    def _read_cfg(self, cfg):
        with open(cfg) as f:
            lines = f.readlines()

            def parse_folder(line):
                items = line.split("=")[1].split("|")
                items[0] = items[0].strip()  # folder
                items[1] = items[1].strip()  # instrument
                items[2] = int(items[2].strip())  # NCE
                if len(items) >= 4: items[3] = items[3].strip()
                return items

            def get_int(line):
                return int(get_str(line))

            def get_str(line):
                return line.split("=")[1].strip()

            for i in range(len(lines)):
                line = lines[i].strip()
                if line.startswith("model"):
                    self.model = get_str(line)
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

                elif line.startswith("predict_input"):
                    self.predict_input, self.predict_instrument, self.predict_nce = parse_folder(line)[:3]
                    self.predict_nce = int(self.predict_nce)
                elif line.startswith("predict_output"):
                    self.predict_output = get_str(line)
                elif line.startswith("fasta"):
                    self.fasta = get_str(line)

                # tune_parameters
                elif line.startswith("tune_psmlabels"):
                    line = get_str(line)
                    if line: self.tune_psmlabels = [s.strip() for s in line.split("|")]
                elif line.startswith("tune_epochs"):
                    self.epochs = get_int(line)
                elif line.startswith("tune_batch"):
                    self.train_batch = get_int(line)
                elif line.startswith("n_tune_per_psmlabel"):
                    self.n_tune_per_psmlabel = get_int(line)
                    
                elif line.startswith("test_psmlabels"):
                    line = get_str(line)
                    if line: self.test_psmlabels = [s.strip() for s in line.split("|")]
                elif line.startswith("n_test_per_psmlabel"):
                    self.n_test_per_psmlabel = get_int(line)
