import time

from .library_base import LibraryBase

_mod_dict = {
    "Carbamidomethyl[C]": "CAM",
    "Label_13C(6)15N(2)[K]": "Label:13C(6)15N(2)", #SILAC
    "Label_13C(6)15N(4)[R]": "Label:13C(6)15N(4)", #SILAC
}

class MSP(LibraryBase):
    def __init__(self, pDeepParam):
        super(self.__class__, self).__init__(pDeepParam)
        self.peptide_dict = {}
        
        self.set_precision(10, 1)
        
    def set_precision(self, mass_precision, inten_precision):
        self._mass_precision = mass_precision
        self._inten_precision = inten_precision
        self._min_rel_inten = float("1e-%d"%inten_precision)
        self._mass_template = "{:.%df}"%mass_precision
        self._inten_template = "{:.%df}"%inten_precision
        self._peak_template = "{:.%df} {:.%df}"%(mass_precision, inten_precision)
        
    def _str_mass(self, mass):
        return self._mass_template.format(mass)
        
    def _str_inten(self, inten):
        return self._inten_template.format(inten)
        
    def Open(self, filename):
        self.file = filename
        
    def Close(self):
        self.file = None
        
    def _write_one_peptide(self, _file, seq, mod, charge, pepmass, masses, intens, charges, types, sites, RT, protein):
        def to_msp_mod(mod):
            if not mod: return "0"
            mods = mod.strip(";").split(";")
            modstr = str(len(mods))
            for mod in mods:
                site, modname = mod.split(",")
                site = int(site)
                if site > len(seq): site -= 2
                elif site >= 1: site -= 1
                if modname in _mod_dict: modstr += "(%d,%s,%s)"%(site,seq[site],_mod_dict[modname])
                else: modstr += "(%d,%s,%s)"%(site,seq[site],modname[:modname.rfind('[')])
            return modstr
        
        #Name: AAAAAAAAAAGAAGGRGSGPGR/3_2(0,A,Acetyl)(17,S,Phospho)_26eV
        # MW: 1833.8703
        # Comment: Single Pep=N-Semitryptic Mods=2(0,A,Acetyl)(17,S,Phospho) Fullname=M.AAAAAAAAAAGAAGGRGSGPGR.R Charge=3 Parent=611.2901 Mz_diff=-1.4ppm HCD=25.979772567749eV Scan=50232 Origfile="20120208_EXQ5_KiSh_SA_LabelFree_HeLa_Phospho_Nocodazole_rep2_FT1.raw.FT.hcd.ch.MGF" Sample="KusterSyn_MOD" Protein="sp|Q86U42|PABP2_HUMAN(pre=M,post=R)" Unassign_all=0.1561 Unassigned=0.0227 max_unassigned_ab=0.20 num_unassigned_peaks=1/20 FTResolution=17500 ms2IsolationWidth=1.60 ms1PrecursorAb=201081499.00 Precursor1MaxAb=191984705.53 PrecursorMonoisoMZ=611.2893 Filter="FTMS + p NSI d Full ms2 611.29@hcd25.00 [100.00-1890.00]"
        # Num peaks: 362
        # 462.7422	65190.5	"y11-H3PO4^2/1.4ppm"
        # 463.2442	28082.4	"y11-H3PO4^2+i/2.9ppm"
        # 466.2397	5450.3	"?"
        # 469.2402	671318	"b6/-0.7ppm"
        # 470.241	117243	"b6+i/-5.0ppm"
        # 471.2495	8158.7	"b6+2i/8.0ppm"
        # 471.7496	5915.2	"?"
        # 474.2256	4296.3	"?"
        # 476.2108	53809.1	"y10^2/-1.4ppm"
        # 476.7162	18087.4	"y10^2+i/7.3ppm"
        
        mod = to_msp_mod(mod)
        _file.write("Name: %s/%d_%s\n"%(seq, charge, mod))
        _file.write("MW: %s\n"%(self._str_mass((pepmass-1.007276)*charge)))
        _file.write('Comment: Spec=pDeep Mods={} Fullname=-.{}.- Charge={} Parent={} Protein="{}" RTInSeconds={}\n'.format(mod, seq, charge, pepmass, protein, RT))
        _file.write("Num peaks: %d\n"%len(masses))
        for mz, inten, charge, ion_type, site in zip(masses, intens, charges, types, sites):
            if len(ion_type) > 1:
                ion_type = ion_type[0] + str(site) + ion_type[1:]
            else:
                ion_type += str(site)
            if charge > 1: ion_type += "^%d"%charge
            ion_type += "/0"
            _file.write("%s\t%s\t%s\n"%(self._str_mass(mz), self._str_inten(inten), ion_type))
        _file.write("\n")
        
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}):
        f = open(self.file, "w")
        self.decoy = None
        print("[pDeep Info] updating msp ...")
        count = 0
        start = time.perf_counter()
        
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            seq, mod, charge = pepinfo.split("|")
            charge = int(charge)
            
            pepmass, masses, intens, sites, types, charges, decoy_seq, decoy_mod, decoy_masses = self._calc_ions(seq, mod, charge, intensities)
            
            RT = _prediction.GetRetentionTime(pepinfo)
            
            if seq in peptide_to_protein_dict:
                protein = peptide_to_protein_dict[seq]
            elif pepinfo in self.peptide_dict:
                protein = self.peptide_dict[pepinfo][-1]
            else:
                protein = "pDeep"
                
            count += 1
            self._write_one_peptide(f, seq, mod, charge, pepmass, masses, intens, charges, types, sites, RT, protein)
            if count%10000 == 0:
                print("[MSP UPDATE] {:.1f}%".format(100.0*count/len(_prediction.peptide_intensity_dict)), end="\r")
        print("[MSP UPDATE] 100%: {}".format(self.file))
        f.close()
        print("[pDeep Info] updating msp time = %.3fs"%(time.perf_counter()-start))
        if not self.decoy: print("[pDeep Info] only target peptides were considered!")
        