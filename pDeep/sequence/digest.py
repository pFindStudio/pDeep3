import re

class DigestConfig:
    def __init__(self):
        self._digest_pattern = "[KR]" # Trypsin/P = "([KR](?=[^P]))", AspN = "\w(?=D)"
        self.digest_pattern = self._digest_pattern
        self.min_len = 7
        self.max_len = 30
        self.max_miss_cleave = 2
        self.cleave_type = "full"

    @property
    def digest_pattern(self):
        return self._digest_pattern

    @digest_pattern.setter
    def digest_pattern(self, _pattern):
        self._digest_pattern = _pattern
        self.regex = re.compile(_pattern)

def digest(protein, pepset, digest_config):
    if digest_config.cleave_type == "full":
        return digest_full(protein, pepset, digest_config.regex, digest_config.min_len, digest_config.max_len,
                           digest_config.max_miss_cleave)
    else:
        return pepset


def digest_full(protein, pepset, regex, min_len=7, max_len=30, max_miss_cleave=2):
    seq = protein.seq
    digest_sites = [m.start()+1 for m in regex.finditer(seq)]
    digest_sites.insert(0,0)
    if digest_sites[-1] != len(seq) - 1:
        # no matter what, cleavage at the last aa
        digest_sites.append(len(seq))

    return cleave_full(seq, pepset, digest_sites, min_len, max_len, max_miss_cleave)


def cleave_full(seq, seq_set, sites, min_len=7, max_len=30, max_miss_cleave=2):
    def add_Nterm_M_loss(seq_set, sub_seq):
        if sub_seq[0] == "M" and len(sub_seq) - 1 >= min_len and len(sub_seq) - 1 <= max_len:
            seq_set.add(sub_seq[1:])

    for i in range(len(sites)):
        for msclv in range(max_miss_cleave + 1):
            if i + msclv + 1 >= len(sites): break
            sub_seq = seq[sites[i]:sites[i + msclv + 1]]
            if len(sub_seq) > max_len or len(sub_seq) < min_len: continue
            if i == 0:
                add_Nterm_M_loss(seq_set, sub_seq)
            seq_set.add(sub_seq)
    return seq_set
