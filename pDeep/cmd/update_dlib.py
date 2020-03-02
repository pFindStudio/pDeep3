import sys
import time
import numpy as np

from ..parameter import pDeepParameter
from ..config import pDeep_config as pDconfig
from ..prediction import pDeepPrediction
from ..spectral_library.encyclopedia.dlib import DLIB

from .. import load_data as load_data
from .. import model_tf as model

from . import tune_and_predict

argd = {}
for i in range(1, len(sys.argv), 2):
    argd[sys.argv[i]] = sys.argv[i+1]

dlib = DLIB()
dlib.Open(argd['-dlib'])
prediction = tune_and_predict.run(argd['-cfg'], dlib.GetAllPeptides())
dlib.UpdateByPrediction(prediction)
dlib.Close()
