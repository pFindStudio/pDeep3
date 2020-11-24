import sys
import os

from ..parameter import pDeepParameter
from .. import evaluate as evaluate

from . import tune_and_predict

# train_folder = sys.argv[1]
# test_folder = sys.argv[2]

train_Vls = ("/home/pfind/pDeepDev/pDeepData/phospho/PhosSynVelosptmRS/train", 40, "Velos")
train_QE = ("/home/pfind/pDeepDev/pDeepData/phospho/Olsen-NC-QEHFX-28/train", 28, "QE")
train_Lumos = ("/home/pfind/pDeepDev/pDeepData/phospho/Olsen-NC-QEHFX-28/train", 32, "Lumos")
test_Vls = ("/home/pfind/pDeepDev/pDeepData/phospho/PhosSynVelosptmRS/test", 40, "Velos")
test_QE = ("/home/pfind/pDeepDev/pDeepData/phospho/Olsen-NC-QEHFX-28/test", 28, "QE")

train_WB = "/home/pfind/pDeepDev/pDeepData/phospho/WenBoPhos/train", 30, "QE"
test_WB = "/home/pfind/pDeepDev/pDeepData/phospho/WenBoPhos/test", 30, "QE"

param = pDeepParameter()
param.model = "tmp/model/pretrain-180921-modloss-mod8D.ckpt"
# param.model = "tmp/model/pretrain-phos.ckpt"
param.RT_model = "RT-model.ckpt"

param.fixmod = "Carbamidomethyl[C],Oxidation[M]".split(",")
param.varmod = "Phospho[Y],Phospho[S],Phospho[T]".split(",")
param.predict_instrument = "QE"
param.predict_nce = 30
param.min_varmod = 1
param.max_varmod = 3
param.epochs = 50
param.dropout = 0.2
param.n_tune_per_psmlabel = 1000000
param.n_test_per_psmlabel = param.n_tune_per_psmlabel

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
    param.tune_psmlabels.extend(psmlabels)
    param.tune_nces.extend(nces)
    param.tune_instruments.extend(instruments)
    
def add_to_test(param, psmlabels, nces, instruments):
    param.test_psmlabels.extend(psmlabels)
    param.test_nces.extend(nces)
    param.test_instruments.extend(instruments)
    
add_to_train(param, *get_psmlabels(*train_Vls))
add_to_train(param, *get_psmlabels(*train_QE))
add_to_train(param, *get_psmlabels(*train_Lumos))
add_to_test(param, *get_psmlabels(*test_Vls))
add_to_test(param, *get_psmlabels(*test_QE))
# add_to_train(param, *get_psmlabels(*train_WB))
# add_to_test(param, *get_psmlabels(*test_WB))

tune_and_predict.init_config(param)

pdeep, _ = tune_and_predict.tune(param)

pdeep.SaveModel("tmp/model/pretrain-phos.ckpt")


