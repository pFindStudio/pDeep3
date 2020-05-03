import sys
import time
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import numpy as np
import tensorflow as tf
from scipy.stats import pearsonr

from ..parameter import pDeepParameter
from ..bucket import peptide_as_key, merge_buckets
from ..config import pDeep_config as pDconfig
from .. import load_data as load_data
from .. import evaluate as evaluate
from .. import similarity_calc as sim_calc
from .. import model_tf as model
from ..rt_model import pDeepRTModel
from ..prediction import pDeepPrediction
from ..data_generator import *
        
def load_param(pDeep_cfg):
    return pDeepParameter(pDeep_cfg)
    
def init_config(param):
    mod_config = pDconfig.HCD_CommonMod_Config()
    mod_config.SetFixMod(param.fixmod)
    mod_config.SetVarMod(param.varmod)
    mod_config.SetIonTypes(param.ion_types)
    mod_config.min_var_mod_num = param.min_varmod
    mod_config.max_var_mod_num = param.max_varmod
    param.config = mod_config
    
def init_pdeep(param):
    pdeep = model.pDeepModel(param.config)
    pdeep.epochs = param.epochs
    pdeep.num_threads = param.threads
    pdeep.dropout = param.dropout
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
            sim_names = ['PCC', 'SPC']
            return evaluate.cum_plot([pcc, spc], sim_names, evaluate.thres_list)
        if not param.test_instruments or not param.test_nces:
            test_buckets = load_data.load_plabel_as_buckets(param.test_psmlabels, param.config, nce, instrument, max_n_samples=param.n_test_per_psmlabel)
        else:
            test_buckets = {}
            for psmlabel, instrument, nce in zip(param.test_psmlabels, param.test_instruments, param.test_nces):
                print(psmlabel, instrument, nce)
                test_buckets = merge_buckets(test_buckets, load_data.load_plabel_as_buckets([psmlabel], param.config, nce, instrument, max_n_samples=param.n_test_per_psmlabel))
        eval_model(pdeep, test_buckets)
        print("[pDeep Info] testing time = %.3fs"%(time.perf_counter() - start_time))
        print("\n")
    
    if param.test_RT_psmlabel:
        print("\n############################\ntest of pre-trained RT model")
        start_time = time.perf_counter()
        def eval_model(pdeep_RT, buckets):
            output_buckets = pdeep_RT.Predict(buckets)
            pred_list = np.array([])
            real_list = np.array([])
            for key, val in output_buckets.items():
                pred_list = np.append(pred_list, val)
                real_list = np.append(real_list, buckets[key][-2])
            pcc = pearsonr(pred_list ,real_list)[0]
            print("[pDeep Info] PCC of test RT = %.3f"%pcc)
            return pcc
        test_RT_buckets = load_data.load_RT_file_as_buckets(param.tune_RT_psmlabel, param.config, nce, instrument, max_n_samples=param.n_test_per_psmlabel)
        eval_model(pdeep_RT, test_RT_buckets)
        print("[pDeep Info] testing RT time = %.3fs"%(time.perf_counter() - start_time))
        print("\n")
        

    if param.tune_psmlabels:
        pdeep.batch_size = param.tune_batch
        start_time = time.perf_counter()
        if not param.tune_instruments or not param.tune_nces:
            train_buckets = load_data.load_plabel_as_buckets(param.tune_psmlabels, param.config, nce, instrument, max_n_samples=param.n_tune_per_psmlabel)
        else:
            train_buckets = {}
            for psmlabel, instrument, nce in zip(param.tune_psmlabels, param.tune_instruments, param.tune_nces):
                print(psmlabel, instrument, nce)
                train_buckets = merge_buckets(train_buckets, load_data.load_plabel_as_buckets([psmlabel], param.config, nce, instrument, max_n_samples=param.n_tune_per_psmlabel))
        # train_buckets = load_data.load_plabel_as_buckets(param.tune_psmlabels, param.config, nce, instrument, max_n_samples=param.n_tune_per_psmlabel)
        pdeep.TrainModel(train_buckets, save_as=param.tune_save_as)
        print("[pDeep Info] tuning time = %.3fs"%(time.perf_counter() - start_time))
        print("\n\n")
        
    if param.tune_RT_psmlabel and pdeep_RT:
        try:
            start_time = time.perf_counter()
            train_buckets = load_data.load_RT_file_as_buckets(param.tune_RT_psmlabel, param.config, nce, instrument, max_n_samples=param.n_tune_per_psmlabel)
            pdeep_RT.TrainModel(train_buckets, save_as=param.tune_RT_save_as)
            print("[pDeep Info] RT tuning time = %.3fs"%(time.perf_counter() - start_time))
        except:
            print("[Error] exception in tuning RT: '{}'".format(param.tune_RT_psmlabel))
        print("\n\n")
    
    pdeep.batch_size = param.predict_batch
    if pdeep_RT: pdeep_RT.batch_size = param.predict_batch
    
    if param.tune_psmlabels and param.test_psmlabels:
        print("\n############################\ntest of fine-tuned model")
        start_time = time.perf_counter()
        def eval_model(pdeep, buckets):
            output_buckets = pdeep.Predict(buckets)
            pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, buckets)
            sim_names = ['PCC', 'SPC']
            return evaluate.cum_plot([pcc, spc], sim_names, evaluate.thres_list)
        eval_model(pdeep, test_buckets)
        print("[pDeep Info] testing time = %.3fs"%(time.perf_counter() - start_time))
        print("\n")
        test_buckets = None # release the memory
        
    if param.tune_RT_psmlabel and pdeep_RT and param.test_RT_psmlabel:
        print("\n############################\ntest of fine-tuned RT model")
        start_time = time.perf_counter()
        def eval_model(pdeep_RT, buckets):
            output_buckets = pdeep_RT.Predict(buckets)
            pred_list = np.array([])
            real_list = np.array([])
            for key, val in output_buckets.items():
                pred_list = np.append(pred_list, val)
                real_list = np.append(real_list, buckets[key][-2])
            pcc = pearsonr(pred_list ,real_list)[0]
            print("[pDeep Info] PCC of test RT = %.3f"%pcc)
            return pcc
        eval_model(pdeep_RT, test_RT_buckets)
        print("[pDeep Info] testing RT time = %.3fs"%(time.perf_counter() - start_time))
        print("\n")
    
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
    
def get_prediction(input_peptides, tune_psm = None, raw = None, n_psm_to_tune = 1000, instrument = 'QE', ce = 27, model = "HCD"):
    '''
    @param input_peptides, could be a peptide list [(sequence1, mod1, charge1), (seq2, mod2, charge2), ...] to be predicted, or a file containing tab seperated head "peptide, modinfo, charge, protein".
    @param tune_psm evidence.txt (MaxQuant), .spectra (pFind) or *.psm.txt/*.txt (with tab seperated head "raw_name, scan, peptide, modinfo, charge, RTInSeconds") file for tuning pDeep and pDeepRT. If it is None, the model will not be tuned (default None).
    @param raw, raw file for tuning pDeep and pDeepRT (default None).
    @param instrument, instrument type for prediction (default "QE").
    @param ce, collision energy for prediction (default 27).
    @return prediction, pDeep.prediction.pDeepPrediction object, 
    # example: 
    # ion_types = ['b','y','b-ModLoss','y-ModLoss'] or ['b','y']
    # ion_indices, used_ion_types = prediction.GetIonTypeIndices(ion_types)
    # print(used_ion_types)
    # intensities = GetIntensitiesByIndices('ACDMNLK', '2,Carbamidomethyl[C];4,Oxidation[M]', 3, ion_indices)
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
    
    if type(input_peptides) is list: peptide_list = input_peptides
    else: peptide_list, pep_pro_dict = ReadModSeq(input_peptides)
    _add_mod(peptide_list)
    
    if tune_psm and raw and RawFileReader:
        if tune_psm.endswith('.psm.txt'):
            psmRT = tune_psm
            psmLabel = Run_psmLabel(psmRT, raw)
        elif tune_psm.endswith('evidence.txt') or tune_psm.endswith(".spectra"):
            psmRT = GeneratePSMFile(tune_psm, raw)
            psmLabel = Run_psmLabel(psmRT, raw)
        else:
            psmRT = tune_psm
            psmLabel = Run_psmLabel(psmRT, raw)
    else:
        psmLabel = None
        psmRT = None
    if psmLabel:
        Sort_psmLabel(psmLabel)
    
    param = pDeepParameter()
        
    param = Set_pDeepParam(param, model, instrument=instrument, ce=ce, psmLabel=psmLabel, psmRT=psmRT, fixmod=",".join(mod_set), varmod=None, n_tune=n_psm_to_tune)
    
    return run(param, peptide_list) #return a pDeep.prediction.pDeepPrediction object
    
if __name__ == "__main__":
    input_peptides = [('ACDMNLK', '2,Carbamidomethyl[C];4,Oxidation[M]', 3)]
    ion_types = ['b','y', 'c', 'z']
    # prediction = get_prediction(input_peptides, tune_psm=r"e:\DDATools\MaxQuant_1.6.12.0\test_data\combined\txt\evidence.txt", raw=r"e:\DDATools\MaxQuant_1.6.12.0\test_data\20141010_DIA_20x5mz_700to800.raw")
    prediction = get_prediction(input_peptides, model="EThcD")
    ion_indices, used_ion_types = prediction.GetIonTypeIndices(ion_types)
    print(used_ion_types)
    print(prediction.GetIntensitiesByIndices(*input_peptides[0], ion_indices))
    