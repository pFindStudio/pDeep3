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

dlib = DLIB()
dlib.Open(sys.argv[1])
prediction = tune_and_predict.run(sys.argv[2], dlib.GetAllPeptides())
dlib.UpdateByPrediction(prediction)
dlib.Close()
