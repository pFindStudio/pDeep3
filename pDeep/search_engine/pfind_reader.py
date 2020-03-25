from .reader_base import *

class pFindSpectraReader(ReaderBase):
    def __init__(self):
        super(self.__class__, self).__init__()
    def GetAllPeptides(self):
        self.peptide_dict = {}
        head = self._file.readline().strip().split("\t")
        headidx = dict(zip(head, range(len(head))))
        
        peptide_list = []
        lines = self._file.readlines()
        for line in lines:
            items = line.strip().split("\t")
            if len(items) < 10: break
            if float(items[headidx['Q-value']]) > self.FDR: continue
            seq = items[headidx['Sequence']]
            mod = items[headidx['Modification']].strip(";")
            charge = int(items[headidx['Charge']])
            RT = -1
            scan = int(items[headidx['Scan_No']])
            raw = items[headidx['File_Name']]
            raw = raw[:raw.find('.%d.%d.'%(scan, scan))]
            protein = items[headidx['Proteins']]
            self.peptide_dict['%s|%s|%d'%(seq,mod,charge)] = (modseq, charge, RT, raw, scan, protein)
            peptide_list.append((seq, mod, charge))
        return peptide_list
