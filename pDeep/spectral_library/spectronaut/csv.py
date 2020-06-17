import sqlite3
import zlib
import struct
import numpy as np
import time

from ...utils.mass_calc import PeptideIonCalculator as ion_calc
from ..library_base import LibraryBase, mod_mass_dict
from ..openswath.tsv import OSW_TSV

_mod_dict = {
    "Carbamidomethyl[C]": "[Carbamidomethyl (C)]",
    "Oxidation[M]": "[Oxidation (M)]",
    "Phospho[S]": "[Phospho (STY)]",
    "Phospho[T]": "[Phospho (STY)]",
    "Phospho[Y]": "[Phospho (STY)]",
}

def NCtermLoss(peptide, modinfo):
    if not modinfo: return [""]*len(peptide), [""]*len(peptide)
    moditems = modinfo.strip(";").split(";")
    modlist = []
    for moditem in moditems:
        site, mod = moditem.split(",")
        modlist.append((int(site), mod))
    modlist.sort()
    
    def XTermLoss(peptide, modlst, Nterm = True):
        XtermLoss = [""]*len(peptide)
        for site, mod in modlist:
            if site == 0:
                XtermLoss[0] = mod
            elif site == len(peptide)+1:
                XtermLoss[len(peptide)-1] = mod
            else:
                XtermLoss[site-1] = mod
        if not Nterm: XtermLoss = XtermLoss[::-1]
        for i in range(1, len(XtermLoss)):
            if XtermLoss[i].startswith("Phospho") or XtermLoss[i-1].startswith("Phospho"):
                XtermLoss[i] = "Phospho"
            elif XtermLoss[i].startswith("Oxidation") or XtermLoss[i-1].startswith("Oxidation"):
                XtermLoss[i] = "Oxidation"
        for i in range(len(XtermLoss)):
            if XtermLoss[i] == "Phospho": XtermLoss[i] = "H3PO4"
            elif XtermLoss[i] == "Oxidation": XtermLoss[i] = "H4COS"
        return XtermLoss
        
    NtermLoss = XTermLoss(peptide, modlist, True)
    CtermLoss = XTermLoss(peptide, modlist, False)
    return NtermLoss, CtermLoss

def pDeepFormat2PeptideModSeq(seq, modinfo):
    if not modinfo: return "_" + seq + "_"
    moditems = modinfo.strip(";").split(";")
    modlist = []
    for moditem in moditems:
        site, mod = moditem.split(",")
        modlist.append((int(site), mod))
    modlist.sort(reverse=True)
    for site, mod in modlist:
        if not mod in _mod_dict: return None
        seq = seq[:site] + _mod_dict[mod] + seq[site:]
    return "_" + seq + "_"

def PeptideModSeq2pDeepFormat(PeptideModSeq):
    PeptideModSeq = PeptideModSeq.strip("_")
    site = PeptideModSeq.find('(')
    modlist = []
    while site != -1:
        if PeptideModSeq[site-1] == 'C': modlist.append('%d,%s;'%(site, 'Carbamidomethyl[C]'))
        elif PeptideModSeq[site-1] == 'M': modlist.append('%d,%s;'%(site, 'Oxidation[M]'))
        elif PeptideModSeq[site-1] == 'S': modlist.append('%d,%s;'%(site, 'Phospho[S]'))
        elif PeptideModSeq[site-1] == 'T': modlist.append('%d,%s;'%(site, 'Phospho[T]'))
        elif PeptideModSeq[site-1] == 'Y': modlist.append('%d,%s;'%(site, 'Phospho[Y]'))
        PeptideModSeq = PeptideModSeq[:site] + PeptideModSeq[PeptideModSeq.find(')')+1:]
        site = PeptideModSeq.find('(', site)
    return PeptideModSeq, "".join(modlist)
    
class SPN_CSV(OSW_TSV):
    def __init__(self, pDeepParam = None):
        super(OSW_TSV, self).__init__(pDeepParam)
        
        self.peptide_dict = {}
        self.head = "PrecursorMz	FragmentMz	iRT	RelativeFragmentIntensity	StrippedPeptide	ModifiedPeptide	LabeledPeptide	PrecursorCharge	ProteinName	ProteinId	FragmentType	FragmentCharge	FragmentNumber	FragmentLossType".split("\t")
        self.headidx = dict(zip(self.head, range(len(self.head))))
        
        self.set_precision(10, 4)
        self.decoy = None
        self.col_sep = ","
        
    def _init_row(self):
        return ["0"] * len(self.head)
        
    def _str_inten(self, inten):
        return self._inten_template.format(inten/self.inten_base)
    
    def GetAllPeptides(self):
        start = time.perf_counter()
        
        f = open(self.tsv_file)
        _head = f.readline().strip().split(self.col_sep)
        _headidx = dict(zip(_head, range(len(_head))))
        
        def _get(items, col):
            return items[_headidx[col]]
        
        self.peptide_dict = {}
        peptide_list = []
        while True:
            line = f.readline()
            if not line: break
            items = line.strip().split(self.col_sep)
            mod_seq = _get(items, "ModifiedPeptide")
            seq, mod = PeptideModSeq2pDeepFormat(mod_seq)
            if seq is None: continue
            charge = int(_get(items, "PrecursorCharge"))
            if not peptide_list or peptide_list[-1] != (seq, mod, charge):
                if "Tr_recalibrated" in _headidx: RT = float(_get(items, "Tr_recalibrated"))*60
                elif "RetentionTime" in _headidx: RT = float(_get(items, "RetentionTime"))*60
                elif "iRT" in _headidx: RT = float(_get(items, "iRT"))*60
                else: RT = 0
                protein = _get(items, "ProteinName")
                pro_list = protein.split("/")[1:]
                proteins = []
                for pro in pro_list:
                    if not pro.startswith("DECOY_"): proteins.append(pro)
                protein = "/".join(proteins)
                peptide_list.append((seq, mod, charge))
                self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = (mod_seq, charge, RT, '', -1, protein) # scan=-1, unknown
        print("reading tsv time = %.3fs"%(time.perf_counter() - start))
        f.close()
        
        return peptide_list
        
    def _write_one_peptide(self, _file, seq, mod, pre_charge, pepmass, masses, intens, charges, types, sites, RT, protein, pep_count, transition_count):
        labeled_seq = pDeepFormat2PeptideModSeq(seq, mod)
        if not labeled_seq: return pep_count, transition_count
        
        NtermLoss,CtermLoss = NCtermLoss(seq, mod)
        
        #"PrecursorMz	FragmentMz	iRT	RelativeFragmentIntensity	StrippedPeptide	ProteinName	ModifiedPeptide	PrecursorCharge	ProteinId	FragmentType	FragmentCharge	FragmentNumber	FragmentLossType"
        items = self._init_row()
        pep_count += 1
        self._set(items, "PrecursorMz", self._str_mass(pepmass))
        self._set(items, "PrecursorCharge", pre_charge)
        self._set(items, "ModifiedPeptide", labeled_seq)
        self._set(items, "LabeledPeptide", labeled_seq)
        self._set(items, "ProteinName", protein)
        self._set(items, "ProteinId", protein)
        self._set(items, "StrippedPeptide", seq)
        self._set(items, "iRT", RT/60)
        for mz, inten, charge, ion_type, site in zip(masses, intens, charges, types, sites):
            transition_count += 1
            if '-' in ion_type:
                # loss_type = ion_type[ion_type.find('-')+1:].lower()
                # loss_type = "H3PO4"
                if ion_type[0] == 'b': loss_type = NtermLoss[site-1]
                elif ion_type[0] == 'y': loss_type = CtermLoss[site-1]
                if not loss_type: loss_type = ion_type[ion_type.find('-')+1:].lower()
            else:
                loss_type = "noloss"
                
            self._set(items, "FragmentType", ion_type[0])
            self._set(items, "FragmentCharge", charge)
            self._set(items, "FragmentMz", self._str_mass(mz))
            self._set(items, "RelativeFragmentIntensity", self._str_inten(inten))
            self._set(items, 'FragmentNumber', site)
            self._set(items, 'FragmentLossType', loss_type)
            
            _file.write(self.col_sep.join(items)+"\n")
        return pep_count, transition_count
        
        
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}):
        self.decoy = False
        f = open(self.tsv_file, "w")
        f.write(self.col_sep.join(self.head)+"\n")
        print("[pDeep Info] updating csv ...")
        transition_count = 0
        count = 0
        start = time.perf_counter()
        
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            seq, mod, charge = pepinfo.split("|")
            charge = int(charge)
            
            pepmass, masses, intens, sites, types, charges, decoy_seq, decoy_mod, decoy_masses = self._calc_ions(seq, mod, charge, intensities)
            
            RT = _prediction.GetRetentionTime(pepinfo)
            
            if seq in peptide_to_protein_dict:
                protein = peptide_to_protein_dict[seq]
                protein = ";".join(protein.split("/"))
            else:
                protein = "pDeep"
                
            count, transition_count = self._write_one_peptide(f, seq, mod, charge, pepmass, masses, intens, charges, types, sites, RT, protein, count, transition_count)
            if count%10000 == 0:
                print("[CSV UPDATE] {:.1f}%".format(100.0*count/len(_prediction.peptide_intensity_dict)), end="\r")
        print("[CSV UPDATE] 100%: {}".format(self.tsv_file))
        f.close()
        print("[pDeep Info] updating csv time = %.3fs"%(time.perf_counter()-start))
        if not self.decoy: print("[pDeep Info] only target transitions were generated!")

