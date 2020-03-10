import os
import sys
import time

import pDeep.config.pDeep_config as fconfig
import pDeep.rt_model as model
from pDeep.bucket import merge_buckets, print_buckets, count_buckets
from pDeep.load_data import load_RT_file_as_buckets

out_model = 'tmp/model/RT_model.ckpt'
epochs = 2
n = 100

mod_config = fconfig.HCD_CommonMod_Config()
mod_config.time_step = 100
mod_config.min_var_mod_num = 0
mod_config.max_var_mod_num = 3

pdeep_rt = model.pDeepRTModel(mod_config)

pdeep_rt.learning_rate = 0.001
pdeep_rt.layer_size = 256
pdeep_rt.batch_size = 1024
pdeep_rt.dropout = 0.0
pdeep_rt.BuildModel(aa_size=82, mod_size=mod_config.GetModFeatureSize() * 2, output_size=mod_config.GetTFOutputSize(),
                 nlayers=1)

pdeep_rt.epochs = epochs

RTfile = sys.argv[1]

start_time = time.perf_counter()

buckets = load_RT_file_as_buckets(RTfile, mod_config, max_n_samples = n)

print('[I] train data:')
print_buckets(buckets, print_peplen=False)
buckets_count = count_buckets(buckets)
print(buckets_count)
print(buckets_count["total"])

load_time = time.perf_counter()

pdeep_rt.TrainModel(buckets, save_as=out_model)

predict_buckets = pdeep_rt.Predict(buckets)
import numpy as np
pred_list = np.array([])
real_list = np.array([])
for key, val in predict_buckets.items():
    pred_list = np.append(pred_list, val)
    real_list = np.append(real_list, buckets[key][-2])
# print(pred_arr, real_arr)
# for pepinfo, pred, real in zip(pep_list, pred_list, real_list):
    # print(pepinfo ,pred, real)
from scipy.stats import pearsonr
print(pearsonr(pred_list ,real_list)[0])

train_time = time.perf_counter()

print("load = {:.3f}s, train = {:.3f}s".format(load_time - start_time, train_time - load_time))
