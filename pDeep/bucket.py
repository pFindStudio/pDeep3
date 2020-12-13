import sys

import numpy as np
from .config.pDeep_config import Common_Config

instrument_list = Common_Config().instrument_list
instrument_num = Common_Config().max_instrument_num
instrument_dict = dict(zip(instrument_list, range(len(instrument_list))))

_feature_name_list = ['x', 'mod_x', 'charge', 'nce', 'instrument', 'y', 'pepinfo']
_feature_name_dict = dict(zip(_feature_name_list, range(len(_feature_name_list))))
bucket_item_dict = _feature_name_dict

def change_instrument_nce(buckets, instrument, nce):
    nce_idx = _feature_name_dict['nce']
    ins_idx = _feature_name_dict['instrument']
    ins_feature = np.zeros((1, instrument_num), dtype=np.int8)
    ins_feature[0, instrument_dict[instrument]] = 1
    for key, val in buckets.items():
        buckets[key][nce_idx] = np.ones(val[nce_idx].shape[0])*nce/100.0
        buckets[key][ins_idx] = np.repeat(ins_feature, val[ins_idx].shape[0], axis=0)
    return buckets

def get_data(one_bucket, name):
    return one_bucket[_feature_name_dict[name]]

def set_mod_zero_buckets(buckets):
    mod_idx = _feature_name_dict['mod_x']
    for key, value in buckets.items():
        buckets[key][mod_idx] = np.zeros(value[mod_idx].shape)
    return buckets

def peptide_as_key(pep_buckets, predict_buckets):
    peptide_prediction_dict = {}
    for key, value in pep_buckets.items():
        predictions = predict_buckets[key][-1]
        for i in range(predictions.shape[0]):
            pepinfo = value[-1][i] # peptide|modinfo|charge, for example 'ABCDEF|3,Carbamidomethyl[C]|2'
            peptide_prediction_dict[pepinfo] = predictions[i]
    return peptide_prediction_dict

# def write_buckets_mgf(outfile, buckets, predict_buckets, fconfig, ioncalc, iontypes=['b{}', 'y{}']):
    # def write_one(f, pepinfo, pred):
        # peptide, mod, charge = pepinfo.split("|")
        # f.write('BEGIN IONS\n')
        # f.write('TITLE=' + pepinfo + '\n')
        # f.write('CHARGE=' + charge + '+\n')
        # pre_charge = int(charge)
        # f.write('pepinfo=' + pepinfo + '\n')

        # ions = {}
        # modmass, lossmass, modname = ioncalc.calc_mod_mass_list(peptide, mod)
        # bions = ioncalc.calc_b_ions(peptide, modmass)
        # pepmass = ioncalc.calc_pepmass_from_b(peptide, modmass, bions)
        # yions = ioncalc.calc_y_from_b(bions, pepmass)

        # if 'b{}' in fconfig.ion_types and 'b{}' in iontypes: ions['b{}'] = bions
        # if 'y{}' in fconfig.ion_types and 'y{}' in iontypes: ions['y{}'] = yions
        # if 'c{}' in fconfig.ion_types and 'c{}' in iontypes: ions['c{}'] = ioncalc.calc_c_from_b(bions)
        # if 'z{}' in fconfig.ion_types and 'z{}' in iontypes: ions['z{}'] = ioncalc.calc_z_from_b(bions, pepmass)
        # if 'b{}-ModLoss' in fconfig.ion_types and 'b{}-ModLoss' in iontypes: ions[
            # 'b{}-ModLoss'] = ioncalc.calc_Nterm_modloss(bions, lossmass, modname)
        # if 'y{}-ModLoss' in fconfig.ion_types and 'y{}-ModLoss' in iontypes: ions[
            # 'y{}-ModLoss'] = ioncalc.calc_Cterm_modloss(yions, lossmass, modname)

        # max_charge = fconfig.max_ion_charge if pre_charge >= fconfig.max_ion_charge else pre_charge

        # peak_list = []

        # for ion_type in ions.keys():
            # x_ions = np.array(ions[ion_type])
            # for charge in range(1, max_charge + 1):
                # intens = pred[:, fconfig.GetIonIndexByIonType(ion_type, charge)]
                # f.write('{}={}\n'.format(ion_type.format("+" + str(charge)),
                                         # ','.join(['%.5f' % inten for inten in intens])))
                # peak_list.extend(zip(x_ions / charge + ioncalc.base_mass.mass_proton, intens))

        # pepmass = pepmass / pre_charge + ioncalc.base_mass.mass_proton
        # f.write("PEPMASS=%.5f\n" % pepmass)

        # peak_list.sort()
        # for mz, inten in peak_list:
            # if inten > 1e-8: f.write("%f %.8f\n" % (mz, inten))

        # f.write('END IONS\n')

    # with open(outfile, 'w') as f:
        # for key, value in buckets.items():
            # preds = predict_buckets[key][-1]
            # for i in range(value[-1].shape[0]):
                # write_one(f, value[-1][i], preds[i])


# def write_buckets(outfile, buckets, predict_buckets, iontypes=['b+1', 'b+2', 'y+1', 'y+2']):
    # def write_one(f, pepinfo, pred):
        # f.write('BEGIN IONS\n')
        # f.write('pepinfo=' + pepinfo + '\n')
        # for i in range(len(iontypes)):
            # f.write('{}={}\n'.format(iontypes[i], ','.join(['%.5f' % inten for inten in pred[:, i]])))
        # f.write('END IONS\n')

    # with open(outfile, 'w') as f:
        # for key, value in buckets.items():
            # preds = predict_buckets[key][-1]
            # for i in range(value[-1].shape[0]):
                # write_one(f, value[-1][i], preds[i])


# write_predict = write_buckets

def print_buckets(buckets, print_peplen=True, print_file=sys.stdout):
    total_size = 0
    for key, value in buckets.items():
        if print_peplen:
            str = '[I] '
            str += 'peplen = %d' % key
            for i in range(len(value)):
                str += ', x{}.shape = {}'.format(i, value[i].shape)
            print(str, file=print_file)
        total_size += value[0].shape[0]
    print('[I] total data size = {}'.format(total_size), file=print_file)


def count_buckets(buckets):
    ret = {}
    ret["total"] = 0
    for key, value in buckets.items():
        ret[str(key)] = value[0].shape[0]
        ret["total"] += value[0].shape[0]
    return ret


def merge_buckets(buckets, buckets_new):
    def merge_buckets_tuples(t1, t2):
        ret = []
        for i in range(len(t1)):
            ret.append(np.append(t1[i], t2[i], axis=0))
        return ret

    for key, value in buckets_new.items():
        if key in buckets:
            buckets[key] = merge_buckets_tuples(buckets[key], value)
        else:
            buckets[key] = value
    return buckets


class Bucket_Batch(object):
    def __init__(self, buckets, batch_size=1024, shuffle=True):
        self._buckets = buckets
        self._bucket_keys = np.sort(np.array(list(buckets.keys()), dtype=np.int32))
        self._shuffle = shuffle
        self._batch_size = batch_size
        self._feature_name_list = _feature_name_list
        self._tuple_idx = dict(zip(self._feature_name_list, range(len(self._feature_name_list))))
        self._tuple_idx['peplen'] = -1

    def get_data_from_batch(self, batch, name):
        return batch[self._tuple_idx[name]]

    def generate_batch(self):
        if self._shuffle:
            self._bucket_keys = np.random.permutation(self._bucket_keys)
        for key in self._bucket_keys:
            _cur_bucket = self._buckets[key]
            _cur_idxes = np.arange(len(_cur_bucket[0]))
            if self._shuffle:
                _cur_idxes = np.random.permutation(_cur_idxes)
            for i in range(0, len(_cur_idxes), self._batch_size):
                if i + self._batch_size > len(_cur_idxes):
                    _end = len(_cur_idxes)
                else:
                    _end = i + self._batch_size

                def get_one_batch(bucket_value):
                    ret = []
                    for j in range(len(bucket_value)):
                        ret.append(bucket_value[j][_cur_idxes[i:_end]])
                    return ret

                ret = get_one_batch(_cur_bucket)
                ret.append(np.array([key] * (_end - i), dtype=np.int32))

                yield ret
