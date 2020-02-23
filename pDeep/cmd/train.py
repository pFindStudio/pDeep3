import os
import sys
import time

import pDeep.config.pDeep_config as fconfig
import pDeep.evaluate as evaluate
import pDeep.similarity_calc as sim_calc
from pDeep.bucket import merge_buckets, print_buckets, count_buckets
from pDeep.load_data import load_folder_as_buckets as load_folder

out_folder = './pDeep-models/model-180921-modloss'
# model_name = 'pretrain-180921-modloss-mod8D.ckpt'
model_name = sys.argv[1]

argd = {}
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
pdeep.BuildModel(aa_size=82, mod_size=mod_config.GetModFeatureSize() * 2, output_size=mod_config.GetTFOutputSize(),
                 nlayers=2)

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

plot_folder = os.path.join(out_folder, 'log/plots/%s' % model_name)

try:
    os.makedirs(plot_folder)
except:
    pass

QEMM = "../datasets/Mann-MouseBrain-NNeu-2015-QEHF-27/new_plabel"
QEYM = "../datasets/Mann-FissionYeast-NM-2014-QE-25/new_plabel"
QEHG = "../datasets/Gygi-HEK293-NBT-2015-QE-25/plabel"
QEHpO = "../datasets/Olsen-CellSys-2017/HelaPhos-QE-28/plabel"  # Hela Phospho
EliteHP_CD4 = "../datasets/Pandey-Human-Nature-2014/Adult_CD4Tcells_bRP_Elite_28_CE32/plabel"
VelosHP_CD4 = "../datasets/Pandey-Human-Nature-2014/Adult_CD4Tcells_bRP_Velos_29_CE41/plabel"
VelosHP_Gel_CD4 = "../datasets/Pandey-Human-Nature-2014/Adult_CD4Tcells_Gel_Velos_30_CE41/plabel"
VelosHP_Lung = "../datasets/Pandey-Human-Nature-2014/Adult_Lung_bRP_Velos_12_CE39/plabel"
VelosHK_Rectum = "../datasets/Kuster-Human-Nature-2014/rectum-Velos-30/plabel"
PT25_tr = "../datasets/zengwenfeng-ProteomeTools/plabel/HCD25/train"
PT30_tr = "../datasets/zengwenfeng-ProteomeTools/plabel/HCD30/train"
PT35_tr = "../datasets/zengwenfeng-ProteomeTools/plabel/HCD35/train"

PT25_ts = "../datasets/zengwenfeng-ProteomeTools/plabel/HCD25/test"
PT30_ts = "../datasets/zengwenfeng-ProteomeTools/plabel/HCD30/test"
PT35_ts = "../datasets/zengwenfeng-ProteomeTools/plabel/HCD35/test"

start_time = time.perf_counter()

buckets = {}
if "QE" in instrument_list: buckets = merge_buckets(buckets, load_folder(QEMM, mod_config, nce=27, instrument=strQE,
                                                                         max_n_samples=n))
if "QE" in instrument_list: buckets = merge_buckets(buckets, load_folder(QEHpO, phos_config, 28, strQE, n))
if "QE" in instrument_list: buckets = merge_buckets(buckets, load_folder(QEYM, mod_config, 25, strQE, n))
if "QE" in instrument_list: buckets = merge_buckets(buckets, load_folder(QEHG, mod_config, 25, strQE, n))
if "Lumos" in instrument_list: buckets = merge_buckets(buckets, load_folder(PT25_tr, mod_config, 25, strLumos, n))
if "Lumos" in instrument_list: buckets = merge_buckets(buckets, load_folder(PT30_tr, mod_config, 30, strLumos, n))
if "Lumos" in instrument_list: buckets = merge_buckets(buckets, load_folder(PT35_tr, mod_config, 35, strLumos, n))
if "Elite" in instrument_list: buckets = merge_buckets(buckets, load_folder(EliteHP_CD4, mod_config, 32, strElite, n))
if "Velos" in instrument_list: buckets = merge_buckets(buckets, load_folder(VelosHP_CD4, mod_config, 41, strVelos, n))
if "Velos" in instrument_list: buckets = merge_buckets(buckets,
                                                       load_folder(VelosHP_Gel_CD4, mod_config, 41, strVelos, n))
if "Velos" in instrument_list: buckets = merge_buckets(buckets, load_folder(VelosHP_Lung, mod_config, 39, strVelos, n))
if "Velos" in instrument_list: buckets = merge_buckets(buckets,
                                                       load_folder(VelosHK_Rectum, mod_config, 30, strVelos, n))

print('[I] train data:')
print_buckets(buckets, print_peplen=False)
buckets_count = count_buckets(buckets)
print(buckets_count)
print(buckets_count["total"])

load_time = time.perf_counter()

pdeep.TrainModel(buckets, save_as=os.path.join(out_folder, model_name))

train_time = time.perf_counter()

print('[I] test data:')
test_buckets_PT25 = load_folder(PT25_ts, mod_config, nce=25, instrument='Lumos', max_n_samples=max_n_test)
print_buckets(test_buckets_PT25, print_peplen=False)
test_buckets_PT30 = load_folder(PT30_ts, mod_config, nce=30, instrument='Lumos', max_n_samples=max_n_test)
print_buckets(test_buckets_PT30, print_peplen=False)
test_buckets_PT35 = load_folder(PT35_ts, mod_config, nce=35, instrument='Lumos', max_n_samples=max_n_test)
print_buckets(test_buckets_PT35, print_peplen=False)

_start_test_time = time.perf_counter()

with open('{}/log/train_{}.txt'.format(out_folder, model_name), 'w') as log_out:
    def test(pcc, cos, spc, kdt, SA, saveplot):
        sim_names = ['PCC', 'COS', 'SPC', 'KDT', 'SA']
        print(evaluate.cum_plot([pcc, cos, spc, kdt, SA], sim_names, evaluate.thres_list, saveplot=saveplot,
                                print_file=log_out), file=log_out)


    output_buckets = pdeep.Predict(buckets)
    pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, buckets)
    test(pcc, cos, spc, kdt, SA, saveplot=os.path.join(plot_folder, 'train_trainset.png'))

    output_buckets = pdeep.Predict(test_buckets_PT25)
    pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, test_buckets_PT25)
    test(pcc, cos, spc, kdt, SA, saveplot=os.path.join(plot_folder, 'train_PT25.png'))

    output_buckets = pdeep.Predict(test_buckets_PT30)
    pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, test_buckets_PT30)
    test(pcc, cos, spc, kdt, SA, saveplot=os.path.join(plot_folder, 'train_PT30.png'))

    output_buckets = pdeep.Predict(test_buckets_PT35)
    pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, test_buckets_PT35)
    test(pcc, cos, spc, kdt, SA, saveplot=os.path.join(plot_folder, 'train_PT35.png'))

test_time = time.perf_counter()

print("load = {:.3f}s, train = {:.3f}s, test = {:.3f}s".format(load_time - start_time, train_time - load_time,
                                                               test_time - _start_test_time))
