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
from .generate_predicted_speclib import _from_fasta,_from_tsv,_from_dlib,_from_pqp
from ..data_generator import *
        
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
        pdeep_RT.epochs = param.epochs*2
        pdeep_RT.num_threads = param.threads
        pdeep_RT.dropout = 0
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
        if param.tune_RT_psmlabel:
            pdeep_RT.BuildTransferModel(param.RT_model)
        else:
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
        print("[pDeep Info] testing time = %.3fs"%(time.perf_counter() - start_time))
        print("\n")

    if param.tune_psmlabels:
        pdeep.batch_size = param.tune_batch
        start_time = time.perf_counter()
        train_buckets = load_data.load_plabel_as_buckets(param.tune_psmlabels, param.config, nce, instrument, max_n_samples=param.n_tune_per_psmlabel)
        pdeep.TrainModel(train_buckets, save_as=None)
        print("[pDeep Info] tuning time = %.3fs"%(time.perf_counter() - start_time))
        print("\n\n")
        
    if param.tune_RT_psmlabel and pdeep_RT:
        try:
            start_time = time.perf_counter()
            train_buckets = load_data.load_RT_file_as_buckets(param.tune_RT_psmlabel, param.config, nce, instrument, max_n_samples=param.n_tune_per_psmlabel)
            pdeep_RT.TrainModel(train_buckets, save_as=None)
            print("[pDeep Info] RT tuning time = %.3fs"%(time.perf_counter() - start_time))
        except:
            print("[Error] exception in tuning RT: '{}'".format(param.tune_RT_psmlabel))
        print("\n\n")
    
    pdeep.batch_size = param.predict_batch
    if pdeep_RT: pdeep_RT.batch_size = param.predict_batch
    
    if param.test_psmlabels:
        print("\n############################\ntest of fine-tuned model")
        start_time = time.perf_counter()
        def eval_model(pdeep, buckets):
            output_buckets = pdeep.Predict(buckets)
            pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, buckets)
            sim_names = ['PCC']
            return evaluate.cum_plot([pcc], sim_names, evaluate.thres_list)
        eval_model(pdeep, test_buckets)
        print("[pDeep Info] testing time = %.3fs"%(time.perf_counter() - start_time))
        print("\n")
        test_buckets = None # release the memory
    
    param._pDeepModel = pdeep
    param._pDeepRT = pdeep_RT
    return pdeep, pdeep_RT
     
def predict(pdeep, param, peptide_list = None):
    print("[pDeep Info] loading peptides ...")
    if peptide_list is not None:
        pep_buckets = load_data.load_peptides_as_buckets(peptide_list, param.config, nce=param.predict_nce, instrument=param.predict_instrument)
    else:
        pep_buckets = load_data.load_peptide_file_as_buckets(param.predict_input, param.config, nce=param.predict_nce, instrument=param.predict_instrument)
    start_time = time.perf_counter()
    print("[pDeep Info] predicting ...")
    predict_buckets = pdeep.Predict(pep_buckets)
    print('[pDeep Info] predicting time = {:.3f}s'.format(time.perf_counter() - start_time))
    return pep_buckets,predict_buckets
    
def predict_RT(pdeep_RT, param, pep_buckets):
    start_time = time.perf_counter()
    print("[pDeep Info] predicting RT ...")
    predict_buckets = pdeep_RT.Predict(pep_buckets)
    # print(predict_buckets)
    print('[pDeep Info] predicting time = {:.3f}s'.format(time.perf_counter() - start_time))
    return predict_buckets
    
def run_predict(param, peptide_list):
    pep_buckets, predict_buckets = predict(param._pDeepModel, param, peptide_list)
    if param._pDeepRT:
        RT_buckets = predict_RT(param._pDeepRT, param, pep_buckets)
    else:
        RT_buckets = None
    return pDeepPrediction(param.config, pep_buckets, predict_buckets, RT_buckets)
    
def run(pDeep_cfg, peptide_list = None):
    '''
    @param pDeep_cfg. pDeep_cfg is the parameter file of pDeep, see "/pDeep-tune.cfg" for details.
    
    @return dict with peptide as key, and prediction as value. Example of a peptide: "ACDEFG|2,Carbamidomethyl[C];|2", where the peptide sequence is "ACDEFG", modification is "2,Carbamidomethyl[C];", precursor charge is "2". Prediction is a np.ndarray with shape (n-1, 8), where n is len(sequence), 8 is len(ion_types)*max_ion_charge ( ion_types=[b, y, b-ModLoss, y-ModLoss], max_ion_charge=2, 1+ and 2+ prodoct ions). The order for 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++]. You can use pDeepPrediction.GetPredictionByIonType to get the intensities for a given ion_type and ion_charge. All the prediction order is from N-term to C-term.
    '''
    if type(pDeep_cfg) is not pDeepParameter: param = load_param(pDeep_cfg)
    else: param = pDeep_cfg
    init_config(param)
    pdeep, pdeep_RT = tune(param)
    pep_buckets, predict_buckets = predict(pdeep, param, peptide_list)
    if pdeep_RT:
        RT_buckets = predict_RT(pdeep_RT, param, pep_buckets)
    else:
        RT_buckets = None
    return pDeepPrediction(param.config, pep_buckets, predict_buckets, RT_buckets)
    
def get_prediction(input_file, tune_psm = None, raw = None, instrument = 'QE', ce = 27):
    '''
    @param input_file, a peptide list file containing title "peptide\tmodinfo\tcharge\tprotein".
    @param tune_psm .osw (OpenSWATH), .elib (EncyclopDIA), evidence.txt (MaxQuant) or .spectra (pFind) file for tuning pDeep and pDeepRT. If it is None, the model will not be tuned (default None).
    @param raw, raw file for tuning pDeep and pDeepRT (default None).
    @param instrument, instrument type for prediction (default "QE").
    @param ce, collision energy for prediction (default 27).
    '''
    
    mod_set = set()
    def _add_mod(peptide_list):
        for seq,mod,charge in peptide_list:
            if mod:
                mods = mod.strip(";").split(";")
                for mod in mods:
                    modname = mod[mod.find(',')+1:]
                    if modname not in mod_set: 
                        mod_set.add(modname)
    
    peptide_list, pep_pro_dict = ReadModSeq(input_file)
    _add_mod(peptide_list)
    
    if tune_psm and raw and RawFileReader:
        psmRT = GeneratePSMFile(tune_psm, raw)
        psmLabel = Run_psmLabel(psmRT, raw)
    else:
        psmLabel = None
        psmRT = None
    # if psmLabel:
        # Sort_psmLabel(psmLabel)
    
    
    param = pDeepParameter()
        
    param = Set_pDeepParam(param, instrument=instrument, ce=ce, psmLabel=psmLabel, psmRT=psmRT, fixmod=",".join(mod_set), varmod=None)
    
    prediction = tune_and_predict.run(param, peptide_list) #return a pDeep.prediction.pDeepPrediction object
    
    # example: 
    # ion_types = ['b','y','b-ModLoss','y-ModLoss']
    # indices = prediction.GetIonTypeIndices(ion_types)
    # GetIntensitiesByIndices('ACDMNLK', '2,Carbamidomethyl[C];4,Oxidation[M]', 3, indices)
    
    
    
    