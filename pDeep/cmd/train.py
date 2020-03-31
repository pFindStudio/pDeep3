import os
import sys
import time

import pDeep.config.pDeep_config as fconfig
import pDeep.evaluate as evaluate
import pDeep.similarity_calc as sim_calc
from pDeep.bucket import merge_buckets, print_buckets, count_buckets
from pDeep.load_data import load_folder_as_buckets as load_folder

out_folder = 'tmp/model'
# model_name = 'pretrain-180921-modloss-mod8D.ckpt'
model_name = sys.argv[1]

argd = {}
epochs = 1
n = 100
if len(sys.argv) > 2:
    for i in range(2, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i + 1]
    epochs = int(argd['-epoch'])
    n = int(argd['-n'])

ion_types = ['b{}', 'y{}', 'b{}-ModLoss', 'y{}-ModLoss']
mod_config = fconfig.HCD_CommonMod_Config()
mod_config.SetIonTypes(ion_types)
mod_config.time_step = 100
mod_config.min_var_mod_num = 0
mod_config.max_var_mod_num = 2

phos_config = fconfig.HCD_pho_Config()
phos_config.SetIonTypes(ion_types)
phos_config.time_step = 100
phos_config.min_var_mod_num = 0
phos_config.max_var_mod_num = 1

if model_name.endswith(".ckpt"):
    import pDeep.model_tf as model
else:
    print("Unknown pDeep model!!")
    sys.exit(-1)

pdeep = model.pDeepModel(mod_config)

pdeep.learning_rate = 0.001
pdeep.layer_size = 256
pdeep.batch_size = 1024
pdeep.dropout = 0.4
pdeep.BuildModel(aa_size=82, mod_size=mod_config.GetModFeatureSize() * 2, output_size=mod_config.GetTFOutputSize(), nlayers=2)

pdeep.epochs = epochs
max_n_test = 100 if n > 100 else n

instrument_list = ["QE", "Velos", "Elite", "Lumos", "Fusion"]

strQE = "QE"
strVelos = "Velos"
strElite = "Elite"
strLumos = "Lumos"
strFusion = "Fusion"

# strUnknown = "XXXX"
# strQE = strUnknown
# strVelos = strUnknown
# strElite = strUnknown
# strLumos = strUnknown
# strFusion = strUnknown

# plot_folder = os.path.join(out_folder, 'log/plots/%s' % model_name)

# try:
    # os.makedirs(plot_folder)
# except:
    # pass

tr_test = r"e:\DIAData\Specter\HEK_SpikeP100\108ng"

start_time = time.perf_counter()

buckets = {}
if "QE" in instrument_list: buckets = merge_buckets(buckets, load_folder(tr_test, mod_config, nce=27, instrument=strQE, max_n_samples=n))

print('[I] train data:')
print_buckets(buckets, print_peplen=False)
buckets_count = count_buckets(buckets)
print(buckets_count)
print(buckets_count["total"])

load_time = time.perf_counter()

pdeep.TrainModel(buckets, save_as=os.path.join(out_folder, model_name))

train_time = time.perf_counter()

print('[I] test data:')
test_buckets = load_folder(tr_test, mod_config, nce=25, instrument='Lumos', max_n_samples=max_n_test)

_start_test_time = time.perf_counter()

def test(pcc, cos, spc, kdt, SA, saveplot):
    sim_names = ['PCC', 'COS', 'SPC', 'KDT', 'SA']
    print(evaluate.cum_plot([pcc, cos, spc, kdt, SA], sim_names, evaluate.thres_list, saveplot=saveplot))

output_buckets = pdeep.Predict(test_buckets)
pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, test_buckets)
test(pcc, cos, spc, kdt, SA, saveplot=None)

test_time = time.perf_counter()

print("load = {:.3f}s, train = {:.3f}s, test = {:.3f}s".format(load_time - start_time, train_time - load_time, test_time - _start_test_time))
