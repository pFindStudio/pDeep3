from .reader_base import *

def PeptideModSeq2pDeepFormat(PeptideModSeq, fixed_C = True):
    PeptideModSeq = PeptideModSeq.strip("_")
    modlist = []
    if PeptideModSeq.startswith('(Acetyl (Protein N-term))'): 
        modlist.append((0, 'Acetyl[ProteinN-term]'))
        PeptideModSeq = PeptideModSeq[len('(Acetyl (Protein N-term))'):]
    site = PeptideModSeq.find('(')
    while site != -1:
        if not fixed_C and PeptideModSeq[site-1:].startswith('C(Carbamidomethyl (C))'): modlist.append((site, 'Carbamidomethyl[C]'))
        elif PeptideModSeq[site-1:].startswith('M(Oxidation (M))'): modlist.append((site, 'Oxidation[M]'))
        elif PeptideModSeq[site-1:].startswith('S(Phospho (S))'): modlist.append((site, 'Phospho[S]'))
        elif PeptideModSeq[site-1:].startswith('T(Phospho (T))'): modlist.append((site, 'Phospho[T]'))
        elif PeptideModSeq[site-1:].startswith('Y(Phospho (Y))'): modlist.append((site, 'Phospho[Y]'))
        else: return None, None
        PeptideModSeq = PeptideModSeq[:site] + PeptideModSeq[PeptideModSeq.find('))')+2:]
        site = PeptideModSeq.find('(', site)
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
    modseq = "_(Acetyl (Protein N-term))SGGVYGCM(Oxidation (M))EVGALVFDIGSYTVR_"
    print(PeptideModSeq2pDeepFormat(modseq))
