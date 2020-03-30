import os

from .bucket import merge_buckets
from .featurize import Seq2Tensor, Seq2Tensor_noCheck, to_numpy


def load_plabel_as_buckets(filenames, config, nce, instrument, max_n_samples=10000000000):
    seq2vec = Seq2Tensor(conf=config, prev=1, next=1)
    seq2vec.max_samples = max_n_samples
    buckets = {}
    count = 0
    for filename in filenames:
        count += 1
        print("%dth psmlabel" % count, end="\r")
        _buckets = seq2vec.Featurize_buckets(filename, nce, instrument)
        buckets = merge_buckets(buckets, _buckets)
    return buckets
    
    
def load_RT_file_as_buckets(filename, config, nce = None, instrument = None, max_n_samples=10000000000):
    seq2vec = Seq2Tensor(conf=config, prev=1, next=1)
    seq2vec.max_samples = max_n_samples
    return seq2vec.Featurize_RT_buckets(filename, nce, instrument)


def load_plabel_as_feature_list(filenames, config, nce, instrument, max_n_samples=10000000000):
    seq2vec = Seq2Tensor(conf=config, prev=1, next=1)
    seq2vec.max_samples = max_n_samples
    feature_list = []
    count = 0
    for filename in filenames:
        count += 1
        print("%dth psmlabel" % count, end="\r")
        feature_list.extend(seq2vec.Featurize_list_with_end_scan(filename, nce, instrument))
    return feature_list


def feature_list_to_buckets(feature_list, old_buckets={}):
    buckets = {}
    for feature in feature_list:
        if len(feature[-2]) in buckets:
            buckets[len(feature[-2])].append(feature[:-1])
        else:
            buckets[len(feature[-2])] = [feature[:-1]]
    return merge_buckets(to_numpy(buckets), old_buckets)


def load_folder_as_buckets(dataset_folder, config, nce, instrument='QE', max_n_samples=10000000000):
    print("Loading %s .." % dataset_folder)
    filenames = []
    for input_file in os.listdir(dataset_folder):
        if input_file.endswith(".psmlabel") or input_file.endswith(".plabel"):
            filenames.append(os.path.join(dataset_folder, input_file))
    return load_plabel_as_buckets(filenames, config, nce, instrument, max_n_samples)


def load_files_as_buckets(filenames, config, nce, instrument='QE', max_n_samples=10000000000):
    print("Loading data from files...")
    return load_plabel_as_buckets(filenames, config, nce, instrument, max_n_samples)


# format 'peptide	modinfo	charge'
def load_peptide_file_as_buckets(filename, config, nce, instrument='QE'):
    peptide_list = []
    with open(filename) as f:
        head = f.readline().strip().split("\t")
        headidx = dict(zip([name.lower() for name in head], range(len(head))))
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not line: continue
        items = line.split("\t")
        peptide_list.append((items['peptide'], items['modinfo'], items['charge']))
    return load_peptides_as_buckets(peptide_list, config, nce, instrument)


# format list of (peptide,modification,charge)
def load_peptides_as_buckets(peptide_list, config, nce, instrument='QE'):
    seq2vec = Seq2Tensor_noCheck(conf=config, prev=1, next=1)
    buckets = seq2vec.Featurize_buckets_predict(peptide_list, nce, instrument)
    return buckets
