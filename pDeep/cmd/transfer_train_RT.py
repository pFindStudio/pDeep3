import sys
import os

from ..parameter import pDeepParameter
from .. import evaluate as evaluate

from . import tune_and_predict

# train_folder = sys.argv[1]
# test_folder = sys.argv[2]

train_GG = "/Users/zengwenfeng/Workspace/Data/pDeep3/Fynn-GlyGly/psmlabel/train", 30, "QE"
test_GG = "/Users/zengwenfeng/Workspace/Data/pDeep3/Fynn-GlyGly/psmlabel/test", 30, "QE"

param = pDeepParameter()
param.model = "pDeep/tmp/model/pretrain-180921-modloss-mod8D.ckpt"
# param.model = "tmp/model/pretrain-phos.ckpt"
param.RT_model = "pDeep/tmp/model/RT-model.ckpt"

param.fixmod = "Carbamidomethyl[C],Oxidation[M]".split(",")
param.varmod = "GlyGly[K]".split(",")
param.predict_instrument = "QE"
param.predict_nce = 27
param.min_varmod = 1
param.max_varmod = 3
param.epochs = 4
param.dropout = 0.2
param.n_tune_per_psmlabel = 1000000
param.n_test_per_psmlabel = param.n_tune_per_psmlabel
param.tune_RT_psmlabel = []
param.test_RT_psmlabel = []

def get_psmlabels(data_folder, nce = 27, instrument = "QE"):
    ret = []
    instruments = []
    nces = []
    for psmlabel in os.listdir(data_folder):
        if psmlabel.endswith(".psmlabel") or psmlabel.endswith(".plabel"):
            ret.append(os.path.join(data_folder, psmlabel))
            instruments.append(instrument)
            nces.append(nce)
    return ret, nces, instruments
    
def add_to_train(param, psmlabels, nces, instruments):
    param.tune_RT_psmlabel.extend(psmlabels)
    param.tune_nces.extend(nces)
    param.tune_instruments.extend(instruments)
    
def add_to_test(param, psmlabels, nces, instruments):
    param.test_RT_psmlabel.extend(psmlabels)
    param.test_nces.extend(nces)
    param.test_instruments.extend(instruments)
    
add_to_train(param, *get_psmlabels(*train_GG))
add_to_test(param, *get_psmlabels(*test_GG))

tune_and_predict.init_config(param)

_, pdeep_RT = tune_and_predict.tune(param)

pdeep_RT.SaveModel("pDeep/tmp/model/RT-GlyGly.ckpt")


