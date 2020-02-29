
class DigestConfig:
    def __init__(self):
        self.digest_AAs = "KR"
        self.Nterm = False
        self.min_len = 9
        self.max_len = 30
        self.max_miss_cleave = 2
        self.cleave_type = "full"


def digest(protein, pepset, conf):
    if conf.cleave_type == "full":
        return digest_full(protein, pepset, conf.digest_AAs, conf.Nterm, conf.min_len, conf.max_len,
                           conf.max_miss_cleave)
    else:
        return pepset


def digest_full(protein, pepset, digest_AAs="KR", Nterm=False, min_len=9, max_len=30, max_miss_cleave=2):
    seq = protein.seq
    if Nterm: seq = seq[::-1]
    digest_sites = [0]  # no matter what, cleavage at pos before the first aa
    for i in range(len(seq)):
        if seq[i] in digest_AAs:
            digest_sites.append(i + 1)  # i+1? seq=KAAK, sites=0 and 3, result seq is seq[0+1:3+1]
    if digest_sites[-1] != len(seq) - 1:
        # no matter what, cleavage at the last aa
        digest_sites.append(len(seq))

    return cleave_full(seq, pepset, digest_sites, Nterm, min_len, max_len, max_miss_cleave)


def cleave_full(seq, seq_set, sites, Nterm=False, min_len=7, max_len=30, max_miss_cleave=2):
    def add_Nterm_M_loss(seq_set, sub_seq):
        if sub_seq[0] == "M" and len(sub_seq) - 1 >= min_len and len(sub_seq) - 1 <= max_len:
            seq_set.add(sub_seq[1:])

    for i in range(len(sites)):
        for msclv in range(max_miss_cleave + 1):
            if i + msclv + 1 >= len(sites): break
            sub_seq = seq[sites[i]:sites[i + msclv + 1]]
            if len(sub_seq) > max_len or len(sub_seq) < min_len: continue
            if Nterm:
                sub_seq = sub_seq[::-1]
                if i == len(sites) - 1:
                    add_Nterm_M_loss(seq_set, sub_seq)
            else:
                if i == 0:
                    add_Nterm_M_loss(seq_set, sub_seq)
            seq_set.add(sub_seq)
    return seq_set
