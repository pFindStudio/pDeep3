import sys
import time
import numpy as np
import tensorflow as tf

from ..parameter import pDeepParameter
from ..bucket import peptide_as_key
from ..config import pDeep_config as pDconfig
from .. import load_data as load_data
from .. import evaluate as evaluate
from .. import similarity_calc as sim_calc
from .. import model_tf as model
from ..rt_model import pDeepRTModel
from ..prediction import pDeepPrediction
        
def load_param(pDeep_cfg):
    return pDeepParameter(pDeep_cfg)
    
def init_config(param):
    mod_config = pDconfig.HCD_CommonMod_Config()
    mod_config.SetFixMod(param.fixmod)
    mod_config.SetVarMod(param.varmod)
    mod_config.SetIonTypes(param._ion_types)
    mod_config.min_var_mod_num = param.min_varmod
    mod_config.max_var_mod_num = param.max_varmod
    param.config = mod_config
    
def init_pdeep(param):
    pdeep = model.pDeepModel(param.config)
    pdeep.epochs = param.epochs
    pdeep.num_threads = param.threads
    pdeep.dropout = 0
    if param.RT_model:
        pdeep_RT = pDeepRTModel(param.config)
    else:
        pdeep_RT = None
        
    return pdeep, pdeep_RT

def tune(param):
    nce = param.predict_nce
    instrument = param.predict_instrument
    pdeep, pdeep_RT = init_pdeep(param)
    
    if param.tune_psmlabels:
        pdeep.BuildTransferModel(param.model) # load the trainable pre-trained model
    else:
        pdeep.LoadModel(param.model) # no transfer learning
    if pdeep_RT:
        pdeep_RT.LoadModel(param.RT_model)
    
    if param.test_psmlabels:
        print("\n############################\ntest of pre-trained model")
        start_time = time.perf_counter()
        def eval_model(pdeep, buckets):
            output_buckets = pdeep.Predict(buckets)
            pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, buckets)
            sim_names = ['PCC']
            return evaluate.cum_plot([pcc], sim_names, evaluate.thres_list)
        test_buckets = load_data.load_plabel_as_buckets(param.test_psmlabels, param.config, nce, instrument, max_n_samples=param.n_test_per_psmlabel)
        eval_model(pdeep, test_buckets)
        print("testing time = %.3fs"%(time.perf_counter() - start_time))
        print("\n")

    if param.tune_psmlabels:
        pdeep.batch_size = param.tune_batch
        start_time = time.perf_counter()
        train_buckets = load_data.load_plabel_as_buckets(param.tune_psmlabels, param.config, nce, instrument, max_n_samples=param.n_tune_per_psmlabel)
        pdeep.TrainModel(train_buckets, save_as=None)
        print("tuning time = %.3fs"%(time.perf_counter() - start_time))
        train_buckets = None # release the memory
    
    pdeep.batch_size = param.predict_batch
    
    if param.test_psmlabels:
        print("\n############################\ntest of fine-tuned model")
        start_time = time.perf_counter()
        def eval_model(pdeep, buckets):
            output_buckets = pdeep.Predict(buckets)
            pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, buckets)
            sim_names = ['PCC']
            return evaluate.cum_plot([pcc], sim_names, evaluate.thres_list)
        eval_model(pdeep, test_buckets)
        print("testing time = %.3fs"%(time.perf_counter() - start_time))
        print("\n")
        test_buckets = None # release the memory
        
    return pdeep, pdeep_RT
     
def predict(pdeep, param, peptide_list = None):
    if peptide_list is not None:
        pep_buckets = load_data.load_peptides_as_buckets(peptide_list, param.config, nce=param.predict_nce, instrument=param.predict_instrument)
    else:
        pep_buckets = load_data.load_peptide_file_as_buckets(param.predict_input, param.config, nce=param.predict_nce, instrument=param.predict_instrument)
    start_time = time.perf_counter()
    print("predicting ...")
    predict_buckets = pdeep.Predict(pep_buckets)
    print('predicting time = {:.3f}s'.format(time.perf_counter() - start_time))
    return pep_buckets,predict_buckets
    
def predict_RT(pdeep_RT, param, pep_buckets):
    start_time = time.perf_counter()
    print("RT predicting ...")
    predict_buckets = pdeep_RT.Predict(pep_buckets)
    # print(predict_buckets)
    print('predicting time = {:.3f}s'.format(time.perf_counter() - start_time))
    return predict_buckets
    
def run(pDeep_cfg, peptide_list = None):
    '''
    @param pDeep_cfg. pDeep_cfg is the parameter file of pDeep, see "/pDeep-tune.cfg" for details.
    
    @return dict with peptide as key, and prediction as value. Example of a peptide: "ACDEFG|2,Carbamidomethyl[C];|2", where the peptide sequence is "ACDEFG", modification is "2,Carbamidomethyl[C];", precursor charge is "2". Prediction is a np.ndarray with shape (n-1, 8), where n is len(sequence), 8 is len(ion_types)*max_ion_charge ( ion_types=[b, y, b-ModLoss, y-ModLoss], max_ion_charge=2, 1+ and 2+ prodoct ions). The order for 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++]. You can use pDeepPrediction.GetPredictionByIonType to get the intensities for a given ion_type and ion_charge. All the prediction order is from N-term to C-term.
    '''
    param = load_param(pDeep_cfg)
    init_config(param)
    pdeep, pdeep_RT = tune(param)
    pep_buckets, predict_buckets = predict(pdeep, param, peptide_list)
    if pdeep_RT:
        RT_buckets = predict_RT(pdeep_RT, param, pep_buckets)
    else:
        RT_buckets = None
    return pDeepPrediction(param.config, pep_buckets, predict_buckets, RT_buckets)
    
if __name__ == "__main__":
    pdeep_prediction = run(sys.argv[1])
    for peptide, intensities in pdeep_prediction.peptide_intensity_dict.items():
        print("b+ ions of %s ="%peptide, pdeep_prediction.GetIntensitiesByIonType(intensities, "b", 1)) # get b+ ions
    for peptide, intensities in pdeep_prediction.peptide_intensity_dict.items():
        print("y+ ions of %s ="%peptide, pdeep_prediction.GetIntensitiesByIonType(intensities, "y", 1)) # get y+ ions
