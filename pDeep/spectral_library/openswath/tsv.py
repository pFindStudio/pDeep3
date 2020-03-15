import sqlite3
import numpy as np
import time

from ...utils.mass_calc import PeptideIonCalculator as ion_calc
from ..library_base import LibraryBase, mod_mass_dict

from ...config.unimod import unimod_dict
    
_mod_dict = {
    "Carbamidomethyl[C]": "C(UniMod:4)",
    "Oxidation[M]": "M(UniMod:35)",
    "Phospho[S]": "S(UniMod:21)",
    "Phospho[T]": "T(UniMod:21)",
    "Phospho[Y]": "Y(UniMod:21)",
    "SILACnoLabel_13C(6)15N(2)[K]": "K(8.014199)",
    "SILACnoLabel_13C(6)15N(4)[R]": "R(10.008269)",
}
for modname, modmass in mod_mass_dict.items():
    if modname not in _mod_dict:
        _name = modname[:modname.rfind('[')]
        if modname.endswith("N-term]")or modname.endswith("C-term]"): _aa = '.'
        else: _aa = modname[-2]
        if _name in unimod_dict:
            _mod_dict[modname] = '%s(UniMod:%d)'%(_aa, unimod_dict[_name])
        else:
            print('[Warning] modification "%s" not found in unimod.xml'%modname)
            _mod_dict[modname] = '%s(%f)'%(_aa, modmass)

def pDeepFormat2PeptideModSeq(seq, modinfo):
    if not modinfo: return seq
    moditems = modinfo.strip(";").split(";")
    modlist = []
    for moditem in moditems:
        site, mod = moditem.split(",")
        modlist.append((int(site), mod))
    modlist.sort(reverse=True)
    for site, mod in modlist:
        if not mod in _mod_dict:
            print('[E] No {} in _mod_dict in openswath/pqp.py'.format(mod))
            return None
        if site == 0 or mod.endswith('N-term]'): seq = _mod_dict[mod]+seq
        else:
            seq = seq[:site-1] + _mod_dict[mod] + seq[site:]
    return seq

def PeptideModSeq2pDeepFormat(PeptideModSeq):
    site = PeptideModSeq.find('(')
    modlist = []
    while site != -1:
        if PeptideModSeq[site-1:].startswith('C(UniMod:4)'): modlist.append('%d,%s'%(site, 'Carbamidomethyl[C]'))
        elif PeptideModSeq[site-1:].startswith('M(UniMod:35)'): modlist.append('%d,%s'%(site, 'Oxidation[M]'))
        elif PeptideModSeq[site-1:].startswith('S(UniMod:21)'): modlist.append('%d,%s'%(site, 'Phospho[S]'))
        elif PeptideModSeq[site-1:].startswith('T(UniMod:21)'): modlist.append('%d,%s'%(site, 'Phospho[T]'))
        elif PeptideModSeq[site-1:].startswith('Y(UniMod:21)'): modlist.append('%d,%s'%(site, 'Phospho[Y]'))
        else: return None, None
        PeptideModSeq = PeptideModSeq[:site] + PeptideModSeq[PeptideModSeq.find(')')+1:]
        site = PeptideModSeq.find('(', site)
    return PeptideModSeq, ";".join(modlist)

class OSW_TSV(LibraryBase):
    def __init__(self):
        self._ion_calc = ion_calc()
        self.peptide_dict = {}
        self.head = "PrecursorMz	ProductMz	Tr_recalibrated	transition_name	CE	LibraryIntensity	transition_group_id	decoy	PeptideSequence	ProteinName	Annotation	FullUniModPeptideName	PrecursorCharge	PeptideGroupLabel	UniprotID	FragmentType	FragmentCharge	FragmentSeriesNumber	LabelType".split("\t")
        self.headidx = dict(zip(self.head, range(len(self.head))))
        
        self.set_precision(6, 1)
        
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
        
    def _init_row(self):
        items = ["0"] * len(self.head)
        self._set(items, 'CE', '-1')
        self._set(items, 'LabelType', 'light')
        self._set(items, 'decoy', '0')
        return items
        
    def _set(self, items, key, val):
        items[self.headidx[key]] = str(val)
    
    def _get(self, items, key):
        return items[self.headidx[key]]
        
    def Open(self, tsv_file):
        self.tsv_file = tsv_file
        
    def Close(self):
        self.tsv_file = None
    
    def GetAllPeptides(self):
        start = time.perf_counter()
        
        f = open(self.tsv_file)
        _head = f.readline().strip().split("\t")
        _headidx = dict(zip(_head, range(len(_head))))
        
        def _get(items, col):
            return items[_headidx[col]]
        
        self.peptide_dict = {}
        peptide_list = []
        while True:
            line = f.readline()
            if not line: break
            items = line.strip().split("\t")
            if _get(items, 'decoy') == "1": continue
            mod_seq = _get(items, "FullUniModPeptideName")
            seq, mod = PeptideModSeq2pDeepFormat(mod_seq)
            if seq is None: continue
            charge = int(_get(items, "PrecursorCharge"))
            if not peptide_list or peptide_list[-1] != (seq, mod, charge):
                RT = float(_get(items, "Tr_recalibrated"))*60
                protein = _get(items, "ProteinName")
                pro_list = protein.split("/")[1:]
                proteins = []
                for pro in pro_list:
                    if not pro.startswith("DECOY_"): proteins.append(pro)
                protein = str(len(proteins))+"/"+"/".join(proteins)
                peptide_list.append((seq, mod, charge))
                self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = [mod_seq, charge, RT, protein]
        print("reading tsv time = %.3fs"%(time.perf_counter() - start))
        f.close()
        
        return peptide_list
        
    def _write_one_peptide(self, _file, seq, mod, pre_charge, pepmass, masses, intens, charges, types, sites, RT, protein, pep_count, transition_count):
        labeled_seq = pDeepFormat2PeptideModSeq(seq, mod)
        #"PrecursorMz	ProductMz	Tr_recalibrated	transition_name	CE	LibraryIntensity	transition_group_id	decoy	PeptideSequence	ProteinName	Annotation	FullUniModPeptideName	PrecursorCharge	PeptideGroupLabel	UniprotID	FragmentType	FragmentCharge	FragmentSeriesNumber	LabelType"
        items = self._init_row()
        pep_count += 1
        self._set(items, "PrecursorMz", self._str_mass(pepmass))
        self._set(items, "PrecursorCharge", pre_charge)
        self._set(items, "FullUniModPeptideName", labeled_seq)
        self._set(items, "ProteinName", protein)
        self._set(items, "UniprotID", protein)
        self._set(items, "PeptideSequence", seq)
        pep_id = "%d_%s_%d"%(pep_count, labeled_seq, pre_charge)
        self._set(items, "transition_group_id", pep_id)
        self._set(items, "PeptideGroupLabel", pep_id)
        self._set(items, "Tr_recalibrated", RT/60)
        count = 0
        for mz, inten, charge, ion_type, site in zip(masses, intens, charges, types, sites):
            # if abs(mz - pepmass) < 10: continue
            transition_count += 1
            self._set(items, "FragmentType", ion_type)
            self._set(items, "FragmentCharge", charge)
            self._set(items, "ProductMz", self._str_mass(mz))
            self._set(items, "LibraryIntensity", self._str_inten(inten))
            self._set(items, 'FragmentSeriesNumber', site)
            self._set(items, 'Annotation', '%s%d^%d/0'%(ion_type, site, charge) if charge > 1 else '%s%d/0'%(ion_type, site))
            self._set(items, "transition_name", "%d_%s_%d"%(transition_count, labeled_seq, pre_charge))
            
            _file.write("\t".join(items)+"\n")
        return pep_count, transition_count
        
    # peak_selection = "rank" or "absolute"
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}, peak_selection = "absolute", threshold = 0.05, mass_upper = 2000):
        f = open(self.tsv_file, "w")
        f.write("\t".join(self.head)+"\n")
        print("updating tsv ...")
        transition_count = 0
        count = 0
        start = time.perf_counter()
        
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            seq, mod, charge = pepinfo.split("|")
            charge = int(charge)
            max_ion_charge = 2
            masses, pepmass = self._ion_calc.calc_by_and_pepmass(seq, mod, max_ion_charge)
            # print(seq, mod, masses)
            
            b_sites = np.tile(np.arange(1, len(seq)).reshape(-1,1), [1,max_ion_charge])
            y_sites = len(seq)-b_sites
            sites = np.concatenate((b_sites, y_sites), axis=1)
            
            b_types = np.array(['b']*((len(seq)-1)*max_ion_charge)).reshape(-1, max_ion_charge)
            y_types = np.array(['y']*((len(seq)-1)*max_ion_charge)).reshape(-1, max_ion_charge)
            types = np.concatenate((b_types, y_types), axis=1)
            
            charges = np.concatenate([np.full(len(seq)-1, i, dtype=int).reshape(-1,1) for i in range(1, max_ion_charge+1)], axis=1)
            charges = np.concatenate((charges, charges), axis=1)
            
            intens = intensities[:,:masses.shape[1]]
            intens[0, 0:2] = 0 #b1+/b1++ = 0
            intens[-1, 2:] = 0 #y1+/y1++ = 0, do not consider y1/b1 in the library
            masses = masses.reshape(-1)
            intens = intens.reshape(-1)
            sites = sites.reshape(-1)
            types = types.reshape(-1)
            charges = charges.reshape(-1)
            intens = intens/np.max(intens)
            
            if peak_selection == "absolute": 
                masses = masses[intens > threshold]
                sites = sites[intens > threshold]
                types = types[intens > threshold]
                charges = charges[intens > threshold]
                intens = intens[intens > threshold]*10000
            else:
                indices = numpy.argsort(intens)[::-1]
                masses = masses[indices[:threshold]]
                sites = sites[indices[:threshold]]
                types = types[indices[:threshold]]
                charges = charges[indices[:threshold]]
                intens = intens[indices[:threshold]]*10000
                
            intens = intens[masses < mass_upper]
            sites = sites[masses < mass_upper]
            types = types[masses < mass_upper]
            charges = charges[masses < mass_upper]
            
            masses = masses[masses < mass_upper]
            
            RT = _prediction.GetRetentionTime(pepinfo)
            
            pepmass = pepmass / charge + self._ion_calc.base_mass.mass_proton
            if seq in peptide_to_protein_dict:
                protein = peptide_to_protein_dict[seq]
                protein = str(protein.count("/")+1)+"/"+protein
            elif pepinfo in self.peptide_dict:
                protein = self.peptide_dict[pepinfo][-1]
            else:
                protein = "1/pDeep"
                
            count, transition_count = self._write_one_peptide(f, seq, mod, charge, pepmass, masses, intens, charges, types, sites, RT, protein, count, transition_count)
            if count%10000 == 0:
                print("[TSV UPDATE] {:.1f}%".format(100.0*count/len(_prediction.peptide_intensity_dict)), end="\r")
        print("[TSV UPDATE] 100%: {}".format(self.tsv_file))
        f.close()
        print("updating tsv time = %.3fs"%(time.perf_counter()-start))
        
        
        