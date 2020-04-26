import sys
import os

from ..parameter import pDeepParameter
from .. import evaluate as evaluate

from . import tune_and_predict

train_folder = sys.argv[1]
test_folder = sys.argv[2]

param = pDeepParameter()
param.model = "tmp/model/pretrain-180921-modloss-mod8D.ckpt"
param.RT_model = ""

param.fixmod = "Carbamidomethyl[C],Oxidation[M]".split(",")
param.varmod = "Phospho[Y],Phospho[S],Phospho[T]".split(",")
param.predict_instrument = "QE"
param.predict_nce = 28
param.min_varmod = 1
param.max_varmod = 3
param.epochs = 20
param.n_tune_per_psmlabel = 1000000
param.n_test_per_psmlabel = param.n_tune_per_psmlabel

def get_psmlabels(data_folder):
    ret = []
    for psmlabel in os.listdir(data_folder):
        if psmlabel.endswith(".psmlabel") or psmlabel.endswith(".plabel"):
            ret.append(os.path.join(data_folder, psmlabel))
    return ret
    
param.tune_psmlabels = get_psmlabels(train_folder)
param.test_psmlabels = get_psmlabels(test_folder)

tune_and_predict.init_config(param)

pdeep, _ = tune_and_predict.tune(param)

pdeep.SaveModel("tmp/model/pretrain_phospho.ckpt")


