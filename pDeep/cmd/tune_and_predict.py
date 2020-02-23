import sys
import time
import numpy as np

from ..parameter import pDeepParameter
from ..bucket import peptide_as_key
from ..config import pDeep_config as pDconfig
from .. import load_data as load_data
from .. import evaluate as evaluate
from ..import similarity_calc as sim_calc
from ..import model_tf as model
    
class pDeepPrediction:
    def __init__(self, param, peptide_prediction_dict):
        '''
        @param param. Parameters from pDeep.cfg.
        @param peptide_prediction_dict. pDeep predicted intensites for all peptides, stored in a dict where the key is a peptide and the value is predicted intensities (np.ndarray, shape=[n-1,8]) of peptide's fragments.
        '''
        self.param = param
        self.peptide_prediction_dict = peptide_prediction_dict
        
    def GetIntensitiesByIonType(self, intensities, ion_type, ion_charge):
        '''
        @param intensities. intensities is predicted intensities of a peptide, which is a np.ndarray with shape [len(sequence)-1, 8]. The order for 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++].
        @param ion_type. ion_type can be "b","y","b-ModLoss","y-ModLoss".
        @param ion_charge. ion_charge can be 1 and 2.
        
        @return a predicted 1-D intensity np.ndarray for the given ion_type and ion_charge. Note that the order of the intensities is from N-term to C-term. For example, for 2+ b ions, the intensities are [b1++, b2++, b3++, ..., y(n-1)++]; for 2+ y ions, the intensities are [y(n-1)++, y(n-2)++, ..., y2++, y1++]. If ion_type or ion_charge is not in the given list, return np.zeros.
        '''
        idx = self.param.config.GetIonIndexByIonType(ion_type, ion_charge)
        if idx is None: return np.zeros((intensities.shape[0], 1), dtype=np.float32)
        else: return intensities[:,idx] 
    
    def GetIntensities(self, sequence, modinfo, precursor_charge):
        '''
        Get the predicted intensities (np.ndarray with shape=[n-1, 8]) for the given sequence, modinfo, precursor_charge. The order for 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++]. If the peptide is not in predictions, return np.zeros.
        '''
        pepinfo = "%s|%s|%d"%(sequence, modinfo, precursor_charge)
        if pepinfo not in self.peptide_prediction_dict: return np.zeros((len(sequence)-1, len(self.param._ion_types)*self.param.config.max_ion_charge), dtype=np.float32)
        return self.peptide_prediction_dict[pepinfo]
        
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
    return pdeep

def tune(param):
    nce = param.predict_nce
    instrument = param.predict_instrument
    pdeep = init_pdeep(param)
    
    if param.tune_psmlabels:
        pdeep.BuildTransferModel(param.model) # load the pre-trained model
    else:
        pdeep.LoadModel(param.model) # no transfer learning
    
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
        
    return pdeep
     
def predict(pdeep, param):
    pep_buckets = load_data.load_peptide_file_as_buckets(param.predict_input, param.config, nce=param.predict_nce, instrument=param.predict_instrument)
    start_time = time.perf_counter()
    predict_buckets = pdeep.Predict(pep_buckets)
    print('predicting time = {:.3f}s'.format(time.perf_counter() - start_time))
    peptide_prediction_dict = peptide_as_key(pep_buckets, predict_buckets)
    return peptide_prediction_dict
    
def run(pDeep_cfg):
    '''
    @param pDeep_cfg. pDeep_cfg is the parameter file of pDeep, see "/pDeep-tune.cfg" for details.
    
    @return dict with peptide as key, and prediction as value. Example of a peptide: "ACDEFG|2,Carbamidomethyl[C];|2", where the peptide sequence is "ACDEFG", modification is "2,Carbamidomethyl[C];", precursor charge is "2". Prediction is a np.ndarray with shape (n-1, 8), where n is len(sequence), 8 is len(ion_types)*max_ion_charge ( ion_types=[b, y, b-ModLoss, y-ModLoss], max_ion_charge=2, 1+ and 2+ prodoct ions). The order for 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++]. You can use pDeepPrediction.GetPredictionByIonType to get the intensities for a given ion_type and ion_charge. All the prediction order is from N-term to C-term.
    '''
    param = load_param(pDeep_cfg)
    init_config(param)
    pdeep = tune(param)
    return pDeepPrediction(param, predict(pdeep, param))
    
if __name__ == "__main__":
    pdeep_prediction = run(sys.argv[1])
    for peptide, intensities in pdeep_prediction.peptide_prediction_dict.items():
        print("b+ ions of %s ="%peptide, pdeep_prediction.GetIntensitiesByIonType(intensities, "b", 1)) # get b+ ions
    for peptide, intensities in pdeep_prediction.peptide_prediction_dict.items():
        print("y+ ions of %s ="%peptide, pdeep_prediction.GetIntensitiesByIonType(intensities, "y", 1)) # get y+ ions
