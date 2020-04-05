import zlib
import struct
import numpy as np
from base64 import b64decode, b64encode

from ..pyRawFileReader.spectrum_processing import DoCentroid

centroid_tag = '1000127"'
profile_tag = '1000128"'

def CheckCentroid(line):
    item = line[line.find(':')+1:]
    if item.startswith(centroid_tag): return True
    elif item.startswith(profile_tag):
        line = line.replace(':1000128"', ':1000127"')
        line = line.replace('profile', 'centroid')
        return line
    else: return True

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
    
def CentroidOneSpectrum(mz_binary_lines, inten_binary_lines, merge_tol = 0.015):
    def _extract_binary_xml(binary_lines):
        for i, line in enumerate(binary_lines):
            item = line.strip()
            if item.startswith('<binary>'):
                return DecodeDoubleArray(item[item.find('>')+1:item.rfind('<')])
    mz_array = _extract_binary_xml(mz_binary_lines)
    inten_array = _extract_binary_xml(inten_binary_lines)
    mz_array, inten_array = DoCentroid(mz_array, inten_array, merge_tol)
    mz_code = EncodeDoubleArray(mz_array)
    inten_code = EncodeDoubleArray(inten_array)
    
    def _reset_binary_xml(binary_lines, code):
        binary_lines[0] = '<binaryDataArray encodedLength="%d">\n'%len(code)
        binary_lines[-2] = '<binary>%s</binary>'%code
    
    _reset_binary_xml(mz_binary_lines, mz_code)
    _reset_binary_xml(inten_binary_lines, inten_code)
    
def Centroid_mzML(input_mzml, output_mzml):
    with open(input_mzml) as f, open(output_mzml, 'w') as out:
        print(output_mzml)
        count = 0
        while True:
            line = f.readline()
            if not line: break
            item = line.strip()
            centroid_line = CheckCentroid(line)
            if centroid_line is True:
                out.write(line)
            else:
                out.write(centroid_line)
            if item.startswith('<binaryDataArrayList'):
                binary_line_list = []
                while True:
                    line = f.readline()
                    item = line.strip()
                    if item.startswith('</binaryDataArrayList>'):
                        out.write(line)
                        break
                    if item.startswith('<binaryDataArray '):
                        binary_lines = [line]
                        while True:
                            line = f.readline()
                            binary_lines.append(line)
                            item = line.strip()
                            if item.startswith('</binaryDataArray>'):
                                break
                        binary_line_list.append(binary_lines)
                CentroidOneSpectrum(binary_line_list[0], binary_line_list[1])
                out.writelines(binary_line_list[0])
                out.writelines(binary_line_list[1])
                count += 1
                print(count, end='\r')
    print(count)


if __name__ == "__main__":
    import sys
    import time
    start_time = time.perf_counter()
    Centroid_mzML(sys.argv[1], sys.argv[2])
    print("time = %.1fs"%(time.perf_counter()-start_time))
    