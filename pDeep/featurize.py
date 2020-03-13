import warnings

import numpy as np

from .bucket import bucket_item_dict

# warnings.filterwarnings('error')
warnings.simplefilter("default")

base_dtype = np.int8


def mod_feature_idx(idx, peplen):
    if idx == 0 or idx == 1:
        idx == 0
    elif idx >= peplen:
        idx = peplen - 1
    else:
        idx -= 1
    return idx


def to_ndarray(bucket_value):
    ret = []
    for i in range(len(bucket_value[0])):
        if isinstance(bucket_value[0][i], np.ndarray):
            x = np.zeros((len(bucket_value),) + bucket_value[0][i].shape, dtype=bucket_value[0][i].dtype)
            for j in range(len(bucket_value)):
                x[j] = bucket_value[j][i]
            ret.append(x)
        else:
            dtype = base_dtype
            if isinstance(bucket_value[0][i], str):
                dtype = None
            elif isinstance(bucket_value[0][i], float):
                dtype = np.float32
            x = []
            for j in range(len(bucket_value)):
                x.append(bucket_value[j][i])
            ret.append(np.array(x, dtype=dtype))
    return ret


def to_numpy(buckets):
    ret_buckets = {}
    for peplen, bucket_value in buckets.items():
        ret_buckets[peplen] = to_ndarray(bucket_value)
    return ret_buckets


class Seq2Tensor:
    def __init__(self, conf, prev=1, next=1):
        self.config = conf
        self.prev = prev
        self.next = next
        self.__parse_mod__()

        self.max_samples = 100000000
        self.__parse_instrument__()

    def __parse_instrument__(self):
        self.instrument_feature = {}
        for i in range(self.config.max_instrument_num):
            feature = [0] * self.config.max_instrument_num
            feature[i] = 1
            if i < len(self.config.instrument_list):
                self.instrument_feature[self.config.instrument_list[i].lower()] = np.array(feature, dtype=base_dtype)
            if i == self.config.max_instrument_num - 1:
                self.instrument_feature['unknown'] = np.array(feature, dtype=base_dtype)

    def __parse_mod__(self):
        self.mod_feature_size = self.config.GetModFeatureSize()

        def parse_element(modname, elem_str):
            feature = [0] * self.mod_feature_size
            elems = elem_str.split(')')[:-1]
            for elem in elems:
                chem, num = elem.split('(')
                num = int(num)
                if chem in self.config.element_index:
                    idx = self.config.element_index[chem]
                elif chem in self.config.all_mod_dict:
                    idx = -2
                else:
                    if modname in self.config.varmod or modname in self.config.fixmod:
                        warnings.warn("[W] Unsupported element '%s' in modification formula '%s'" % (chem, elem_str))
                    idx = -1
                feature[idx] = num
            return np.array(feature, dtype=base_dtype)

        self.mod_feature = {}
        for modname, elem_str in self.config.all_mod_dict.items():
            self.mod_feature[modname] = parse_element(modname, elem_str.split(' ')[-1])

    def get_mod_x(self, peptide, modinfo):
        mod_x = np.zeros((len(peptide) - 1, self.mod_feature_size * 2), dtype=base_dtype)

        if modinfo:
            moditems = modinfo.strip(";").split(';')

            for mod in moditems:
                modtmp = mod.split(',')
                idx = int(modtmp[0])
                modname = modtmp[1]
                idx = mod_feature_idx(idx, len(peptide))
                if idx < 0: mod_x[idx - 1, self.mod_feature_size:] = self.mod_feature[modname]
                if idx < len(peptide) - 1: mod_x[idx, :self.mod_feature_size] = self.mod_feature[modname]
        return mod_x

    def CountVarMod(self, peptide, modlist):
        var_mod_count = 0
        for idx, modname in modlist:
            if modname in self.config.fixmod:
                continue
            elif modname in self.config.varmod:
                var_mod_count += 1
            else:
                return -1  # unexpected
        return var_mod_count

    def FeaturizeOnePeptide(self, peptide, modinfo):
        if len(peptide) > self.config.time_step + 1: return None

        mod_idx_feature = [np.array([0] * self.mod_feature_size, dtype=base_dtype) for i in range(len(peptide))]
        if modinfo:
            moditems = modinfo.split(';')
            unexpected_mod = False
            modlist = []
            var_mod_count = 0

            for mod in moditems:
                modtmp = mod.split(',')
                idx = int(modtmp[0])
                modname = modtmp[1]
                modlist.append((idx, modname))
                if modname in self.config.fixmod:
                    idx = mod_feature_idx(idx, len(peptide))
                    mod_idx_feature[idx] = self.mod_feature[modname]
                elif modname in self.config.varmod:
                    idx = mod_feature_idx(idx, len(peptide))
                    mod_idx_feature[idx] = self.mod_feature[modname]
                    var_mod_count += 1
                else:
                    unexpected_mod = True
                    break
            if var_mod_count < self.config.min_var_mod_num or var_mod_count > self.config.max_var_mod_num: return None
            if unexpected_mod: return None
            if not self.config.CheckFixMod(peptide, modlist): return None
        if not CheckPeptide(peptide): return None

        x = _seq2vector(peptide, self.prev, self.next)
        mod_x = []
        for site in range(1, len(peptide)):
            mod_x.append(np.append(mod_idx_feature[site - 1], mod_idx_feature[site]))
        return np.array(x), np.array(mod_x)

    def FeaturizeOnePeptide_buckets(self, peptide, modinfo, pre_charge, nce, instrument):
        nce /= 100.0
        peptide_info = "{}|{}|{}".format(peptide, modinfo, pre_charge)
        x = self.FeaturizeOnePeptide(peptide, modinfo)
        if x is None:
            print("[Error] Illegal peptide: {}".format(peptide_info))
            return None
        instrument = instrument.lower()
        if instrument in self.instrument_feature:
            inst_feature = self.instrument_feature[instrument]
        else:
            warnings.warn("[W] Unknown instrument: %s, pDeep uses 'unknown' instrument" % instrument, UserWarning)
            inst_feature = self.instrument_feature['unknown']
        buckets = {len(peptide): [(x[0], x[1], pre_charge, nce, inst_feature, peptide_info)]}
        return to_numpy(buckets)

    def Featurize_buckets_predict(self, peptide_list, nce, instrument):
        nce /= 100.0
        instrument = instrument.lower()
        if instrument in self.instrument_feature:
            inst_feature = self.instrument_feature[instrument]
        else:
            warnings.warn("[W] Unknown instrument: %s, pDeep uses 'unknown' instrument" % instrument, UserWarning)
            inst_feature = self.instrument_feature['unknown']

        buckets = {}

        for peptide, modinfo, pre_charge in peptide_list:
            pre_charge = int(pre_charge)

            # print(items[0])
            x = self.FeaturizeOnePeptide(peptide, modinfo)
            if x is None: continue

            peplen = len(peptide)

            peptide_info = "{}|{}|{}".format(peptide, modinfo, pre_charge)
            if peplen in buckets:
                buckets[peplen].append((x[0], x[1], pre_charge, nce, inst_feature, peptide_info))
            else:
                buckets[peplen] = [(x[0], x[1], pre_charge, nce, inst_feature, peptide_info)]

        return to_numpy(buckets)

    def get_instrument_x(self, instrument):
        instrument = instrument.lower()
        if instrument in self.instrument_feature:
            return self.instrument_feature[instrument]
        else:
            warnings.warn("[W] Unknown instrument: %s, pDeep uses 'unknown' instrument" % instrument, UserWarning)
            return self.instrument_feature['unknown']
        
    def Featurize_RT_buckets(self, ion_file, nce = None, instrument = None):
        f = open(ion_file)

        buckets = {}

        items = f.readline().strip().split('\t')
        headeridx = dict(zip(items, range(len(items))))
        
        if 'RT' in headeridx: RTidx = headeridx['RT']
        elif 'RTInSeconds' in headeridx: RTidx = headeridx['RTInSeconds']
        else: 
            print("[W] no column 'RT' or 'RTInSeconds' in '{}'".format(ion_file))
            print("[W] all 'RT's will be assigned as zero")
            RTidx = None

        sample_count = 0
        while True:
            line = f.readline()
            if line == "": break
            items = line.strip().split("\t")
            peptide = items[headeridx["peptide"]]
            modinfo = items[headeridx["modinfo"]].strip(";")

            x = self.FeaturizeOnePeptide(peptide, modinfo)
            if x is None: continue
            pre_charge = int(items[headeridx["charge"]])
            if RTidx: RT = float(items[RTidx])
            else: RT = 0

            peplen = len(peptide)

            pepinfo = "{}|{}|{}".format(peptide, modinfo, pre_charge)
            if peplen in buckets:
                buckets[peplen].append((x[0], x[1], pre_charge, 0, 0, RT, pepinfo))
            else:
                buckets[peplen] = [(x[0], x[1], pre_charge, 0, 0, RT, pepinfo)]
            sample_count += 1
            if sample_count >= self.max_samples: break
        f.close()

        return to_numpy(buckets)

    def Featurize_buckets(self, ion_file, nce, instrument):
        nce /= 100.0
        f = open(ion_file)
        inst_feature = self.get_instrument_x(instrument)

        buckets = {}

        items = f.readline().strip().split('\t')
        headeridx = dict(zip(items, range(len(items))))

        charge_in_spec = True
        if "charge" in headeridx: charge_in_spec = False

        sample_count = 0
        while True:
            line = f.readline()
            if line == "": break
            type2inten = {}
            allinten = []
            items = line.rstrip("\r\n").split("\t")
            peptide = items[headeridx["peptide"]]
            modinfo = items[headeridx["modinfo"]].strip(";")
            if charge_in_spec:
                pre_charge = int(items[0].split(".")[-3])
            else:
                pre_charge = int(items[headeridx["charge"]])

            x = self.FeaturizeOnePeptide(peptide, modinfo)
            if x is None: continue

            type2inten = {}
            for ion_type in self.config.GetIonTypeNames():
                peaks = items[headeridx[ion_type]]
                if len(peaks) < 2: continue
                peaks = [peak.split(",") for peak in peaks.strip().strip(";").split(";")]
                type2inten.update(dict([(peak[0], float(peak[1])) for peak in peaks]))
            if len(type2inten) < len(peptide): continue

            intenvec = []
            for site in range(1, len(peptide)):
                v = []
                for ion_type in self.config.ion_types:
                    ion_name = self.config.GetIonNameBySite(peptide, site, ion_type)
                    for charge in range(1, self.config.max_ion_charge + 1):
                        ion_name_charge = ion_name + "+{}".format(charge)
                        if ion_name_charge in type2inten:
                            v.append(type2inten[ion_name_charge])
                        else:
                            v.append(0)
                intenvec.append(np.array(v, dtype=np.float32))

            intenvec = np.array(intenvec)
            intenvec /= np.max(intenvec)
            peplen = len(peptide)

            pepinfo = "{}|{}|{}".format(peptide, modinfo, pre_charge)
            if peplen in buckets:
                buckets[peplen].append((x[0], x[1], pre_charge, nce, inst_feature, intenvec, pepinfo))
            else:
                buckets[peplen] = [(x[0], x[1], pre_charge, nce, inst_feature, intenvec, pepinfo)]
            sample_count += 1
            if sample_count >= self.max_samples: break
        f.close()

        return to_numpy(buckets)


class Seq2Tensor_noCheck(Seq2Tensor):
    def __init__(self, conf, prev=1, next=1):
        super(self.__class__, self).__init__(conf, prev, next)

    def FeaturizeOnePeptide(self, peptide, modinfo):
        if not CheckPeptide(peptide):
            print("[W] invalid aa in sequence '%s', ignore this peptide!" % peptide)
            return None

        x = _seq2vector(peptide, self.prev, self.next)
        mod_x = self.get_mod_x(peptide, modinfo)
        return np.array(x), np.array(mod_x)


ASCII = 128

def _AAVec():
    AAVec = [None for _ in range(ASCII)]
    AAVec[ord('A')] = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('C')] = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('D')] = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('E')] = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('F')] = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('G')] = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('H')] = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('I')] = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('K')] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('L')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('M')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('N')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('P')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    AAVec[ord('Q')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    AAVec[ord('R')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    AAVec[ord('S')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
    AAVec[ord('T')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    AAVec[ord('V')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
    AAVec[ord('W')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    AAVec[ord('Y')] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    return AAVec


def _AAIdx():
    AAIdx = [-1 for _ in range(ASCII)]
    AAIdx[ord('A')] = 0
    AAIdx[ord('C')] = 1
    AAIdx[ord('D')] = 2
    AAIdx[ord('E')] = 3
    AAIdx[ord('F')] = 4
    AAIdx[ord('G')] = 5
    AAIdx[ord('H')] = 6
    AAIdx[ord('I')] = 7
    AAIdx[ord('K')] = 8
    AAIdx[ord('L')] = 9
    AAIdx[ord('M')] = 10
    AAIdx[ord('N')] = 11
    AAIdx[ord('P')] = 12
    AAIdx[ord('Q')] = 13
    AAIdx[ord('R')] = 14
    AAIdx[ord('S')] = 15
    AAIdx[ord('T')] = 16
    AAIdx[ord('V')] = 17
    AAIdx[ord('W')] = 18
    AAIdx[ord('Y')] = 19
    return AAIdx


nAAs = 20
AAVec = _AAVec()
AAIdx = _AAIdx()


def CheckPeptide(peptide):
    for aa in peptide:
        if AAVec[ord(aa)] is None:
            return 0
    return 1


def _seq2vector(peptide, prev, next):
    x = []
    base_NtermAA_Count = [0] * nAAs
    base_CtermAA_Count = [0] * nAAs
    for aa in peptide[prev + next:]:
        base_CtermAA_Count[AAIdx[ord(aa)]] += 1
    for seqidx in range(1, len(peptide)):
        v = []
        for i in range(seqidx - prev, seqidx):
            if i < 0:
                v.extend([0] * nAAs)
            else:
                v.extend(AAVec[ord(peptide[i])])
        for i in range(seqidx, seqidx + next):
            if i >= len(peptide):
                v.extend([0] * nAAs)
            else:
                v.extend(AAVec[ord(peptide[i])])
        prev_idx = seqidx - 1 - prev
        next_idx = seqidx + next
        if prev_idx >= 0:
            base_NtermAA_Count[AAIdx[ord(peptide[prev_idx])]] += 1
        if next_idx < len(peptide):
            base_CtermAA_Count[AAIdx[ord(peptide[next_idx])]] -= 1
        v.extend(base_NtermAA_Count)
        v.extend(base_CtermAA_Count)

        if seqidx == 1:
            v.append(1)
        else:
            v.append(0)
        if seqidx == len(peptide) - 1:
            v.append(1)
        else:
            v.append(0)

        x.append(np.array(v, dtype=base_dtype))
    return x


def _mod2vector(peptide, mod_idx_feature):
    mod_x = []
    for site in range(1, len(peptide)):
        mod_x.append(np.append(mod_idx_feature[site - 1], mod_idx_feature[site]))
    return mod_x
