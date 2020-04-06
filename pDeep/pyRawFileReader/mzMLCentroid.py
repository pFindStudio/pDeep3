import zlib
import struct
import numpy as np
from base64 import b64decode, b64encode

from ..pyRawFileReader.spectrum_centroid import DoCentroid

centroid_tag = '1000127"'
profile_tag = '1000128"'

def CheckCentroid(line):
    item = line[line.find(':')+1:]
    if item.startswith(centroid_tag): return True
    elif item.startswith(profile_tag):
        line = line.replace(':1000128"', ':1000127"')
        line = line.replace('profile', 'centroid')
        return line
    else: return None

def DecodeDoubleArray(code_str):
    code = b64decode(code_str)
    # print(len(code))
    bytes = zlib.decompress(code)
    return np.array(struct.unpack("%dd"%(len(bytes)/8), bytes), dtype=np.float64)
    
def DecodeFloatArray(code_str):
    code = b64decode(code_str)
    bytes = zlib.decompress(code)
    return np.array(struct.unpack("%df"%(len(bytes)/4), bytes), dtype=np.float32)
    
def EncodeDoubleArray(lst):
    bytes = struct.pack("%dd"%len(lst), *lst)
    code = zlib.compress(bytes)
    return b64encode(code)
    
def EncodeFloatArray(lst):
    bytes = struct.pack("%df"%len(lst), *lst)
    code = zlib.compress(bytes)
    return b64encode(code)
    
def CentroidOneSpectrum(one_spectrum_lines, merge_tol = 0.015):
    i = 0
    while i < len(one_spectrum_lines):
        line = one_spectrum_lines[i]
        centroid_line = CheckCentroid(line)
        if centroid_line is True: return
        elif type(centroid_line) is str: one_spectrum_lines[i] = centroid_line
        item = line.strip()
        if item.startswith('<binaryDataArrayList'):
            binary_array_start_end = []
            while True:
                i += 1
                line = one_spectrum_lines[i]
                item = line.strip()
                if item.startswith('</binaryDataArrayList>'):
                    binary_end = i
                    break
                if item.startswith('<binaryDataArray '):
                    array_start = i
                    while True:
                        i += 1
                        line = one_spectrum_lines[i]
                        item = line.strip()
                        if item.startswith('</binaryDataArray>'):
                            array_end = i
                            break
                    binary_array_start_end.append((array_start, array_end))
        i += 1
        
    def _extract_binary_xml(start, end):
        for i in range(start, end+1):
            item = one_spectrum_lines[i].strip()
            if item.startswith('<binary>'):
                return DecodeDoubleArray(item[item.find('>')+1:item.rfind('<')])
    mz_array = _extract_binary_xml(*binary_array_start_end[0])
    inten_array = _extract_binary_xml(*binary_array_start_end[1])
    
    mz_array, inten_array = DoCentroid(mz_array, inten_array, merge_tol)
    mz_code = EncodeDoubleArray(mz_array)
    inten_code = EncodeDoubleArray(inten_array)
    
    
    # rewrite defaultArrayLength
    num_peak = mz_array.shape[0]
    first_line = one_spectrum_lines[0]
    len_idx = first_line.rfind('defaultArrayLength')
    idx1 = first_line.find('"', len_idx+17)+1
    idx2 = first_line.find('"', idx1+1)
    # print(mz_array.shape[0], inten_array.shape[0], first_line, first_line[idx1:idx2])
    one_spectrum_lines[0] = first_line[:idx1] + str(num_peak) + first_line[idx2:]
    
    
    def _reset_binary_xml(start, end, code):
        code = code.decode('ascii')
        one_spectrum_lines[start] = '<binaryDataArray encodedLength="%d">\n'%len(code)
        one_spectrum_lines[end-1] = '<binary>%s</binary>\n'%code
    
    _reset_binary_xml(*binary_array_start_end[0], mz_code)
    _reset_binary_xml(*binary_array_start_end[1], inten_code)
    
def Centroid_mzML(input_mzml, output_mzml):
    with open(input_mzml) as f, open(output_mzml, 'w') as out:
        print(output_mzml)
        count = 0
        while True:
            line = f.readline()
            if not line: break
            out.write(line)
            item = line.strip()
            if item.startswith('<spectrumList '):
                total = item[item.find('"')+1:]
                total = int(total[:total.find('"')])
                while True:
                    line = f.readline()
                    item = line.strip()
                    if item.startswith('</spectrumList>'): 
                        out.write(line)
                        break
                    if item.startswith('<spectrum '):
                        one_spectrum = [line]
                    elif item.startswith('</spectrum>'):
                        one_spectrum.append(line)
                        CentroidOneSpectrum(one_spectrum)
                        out.writelines(one_spectrum)
                        count += 1
                        print('%d / %d'%(count, total), end='\r')
                    else:
                        one_spectrum.append(line)
    print(total)


if __name__ == "__main__":
    import sys
    import time
    start_time = time.perf_counter()
    Centroid_mzML(sys.argv[1], sys.argv[2])
    print("time = %.1fs"%(time.perf_counter()-start_time))
    