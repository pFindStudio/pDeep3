import numba
from numba import float64, int32
import numpy as np
import math

# @numba.jit((float64[:], float64[:], float64))
@numba.jit
def DoCentroid(mz_array, inten_array, merge_tol = 0.01):
    ret_mzs = []
    ret_intens = []
    
    start_idx = 0
    end_idx = 0
    while start_idx < mz_array.shape[0]:
        end_idx = __find_sister_peaks(mz_array, inten_array, start_idx, merge_tol)
        if np.sum(inten_array[start_idx:end_idx+1]) <= 3:
            start_idx = end_idx + 1
            continue
        center_mz = __merge_sister_peaks(mz_array, inten_array, start_idx, end_idx)
        center_inten = np.max(inten_array[start_idx:end_idx+1])
        ret_mzs.append(center_mz)
        ret_intens.append(center_inten)
        start_idx = end_idx + 1
    return np.array(ret_mzs), np.array(ret_intens)
    
@numba.jit(int32(float64[:], float64[:], int32, float64), nopython=True)
def __find_sister_peaks(mz_array, inten_array, start_idx, merge_tol):
    end_idx = start_idx
    for i in range(start_idx + 1, mz_array.shape[0]):
        if inten_array[i] == 0:
            end_idx += 1
            break
        if mz_array[i] - mz_array[end_idx] <= merge_tol:
            end_idx += 1
        else:
            break
    return end_idx
    
@numba.jit(float64(float64[:], float64[:], int32, int32), nopython=True)
def __merge_sister_peaks(mz_array, inten_array, start_idx, end_idx):
    # center_idx = np.argmax(inten_array[start_idx:end_idx+1])+start_idx
    center_mass = 0
    weight = 0
    for i in range(start_idx, end_idx+1):
        if inten_array[i] != 0:
            _w = math.log(inten_array[i])
            center_mass += mz_array[i] * _w
            weight += _w
    return center_mass / weight
    
    