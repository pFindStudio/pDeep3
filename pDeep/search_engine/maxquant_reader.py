from .reader_base import *

MQ_mod_dict = {}
MQ_mod_dict['(Acetyl (Protein N-term))'] = 'Acetyl[ProteinN-term]'
MQ_mod_dict['C(Carbamidomethyl (C))'] = 'Carbamidomethyl[C]'
MQ_mod_dict['M(Oxidation (M))'] = 'Oxidation[M]'
MQ_mod_dict['S(Phospho (S))'] = 'Phospho[S]'
MQ_mod_dict['T(Phospho (T))'] = 'Phospho[T]'
MQ_mod_dict['Y(Phospho (Y))'] = 'Phospho[Y]'
MQ_mod_dict['K(GlyGly (K))'] = 'GlyGly[K]'
MQ_mod_dict['(ac)'] = 'Acetyl[ProteinN-term]'
MQ_mod_dict['M(ox)'] = 'Oxidation[M]'
MQ_mod_dict['S(ph)'] = 'Phospho[S]'
MQ_mod_dict['T(ph)'] = 'Phospho[T]'
MQ_mod_dict['Y(ph)'] = 'Phospho[Y]'
MQ_mod_dict['K(gl)'] = 'GlyGly[K]'

def PeptideModSeq2pDeepFormat(PeptideModSeq, fixed_C = True):
    PeptideModSeq = PeptideModSeq.strip("_")
    modlist = []
    if PeptideModSeq.startswith('(Acetyl (Protein N-term))'): 
        modlist.append((0, 'Acetyl[ProteinN-term]'))
        PeptideModSeq = PeptideModSeq[len('(Acetyl (Protein N-term))'):]
    site = PeptideModSeq.find('(')
    while site != -1:
        site_end = PeptideModSeq.find(')')+1
        if site_end < len(PeptideModSeq) and PeptideModSeq[site_end] == ')': site_end += 1
        if PeptideModSeq[site-1:site_end] in MQ_mod_dict: modlist.append((site, MQ_mod_dict[PeptideModSeq[site-1:site_end]]))
        else: return None, None
        PeptideModSeq = PeptideModSeq[:site] + PeptideModSeq[site_end:]
        site = PeptideModSeq.find('(', site)
        site_end = PeptideModSeq.find(')')+1
        if site_end < len(PeptideModSeq) and PeptideModSeq[site_end] == ')': site_end += 1
    if fixed_C:
        site = PeptideModSeq.find('C')
        while site != -1:
            modlist.append((site+1,'Carbamidomethyl[C]'))
            site = PeptideModSeq.find('C',site+1)
    modlist.sort(key=lambda x:x[0])
    return PeptideModSeq, ";".join(["%d,%s"%(site,mod) for site, mod in modlist])

class MaxQuantEvidenceReader(ReaderBase):
    def __init__(self):
        super(self.__class__, self).__init__()
    def GetAllPeptides(self):
        self.peptide_dict = {}
        head = self._file.readline().strip().split("\t")
        headidx = dict(zip([name.lower() for name in head], range(len(head))))
        
        peptide_list = []
        lines = self._file.readlines()
        for line in lines:
            items = line.strip().split("\t")
            modseq = items[headidx['modified sequence']]
            seq, mod = PeptideModSeq2pDeepFormat(modseq)
            if seq is None: continue
            charge = int(items[headidx['charge']])
            RT = float(items[headidx['retention time']])*60
            scan = int(items[headidx['ms/ms scan number']])
            raw = items[headidx['raw file']]
            protein = "/".join([pro for pro in items[headidx['proteins']].split(";")])
            self.peptide_dict['%s|%s|%d'%(seq,mod,charge)] = (modseq, charge, RT, raw, scan, protein)
            peptide_list.append((seq, mod, charge))
        return peptide_list
            
if __name__ == "__main__":
    modseq = "_(ac)SGGVYGCM(ox)EVGALK(gl)VFDIGSYTVR_"
    print(PeptideModSeq2pDeepFormat(modseq))
