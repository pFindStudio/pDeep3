class pFindMGFReader(object):
    def __init__(self, filename, **kwargs):
        self.scan_to_RT_dict = {}
        with open(filename) as f:
            scan = -1
            while True:
                line = f.readline()
                if not line: break
                if line.startswith("TITLE"):
                    if scan == -1: scan = int(line.split(".")[-4])
                elif line.startswith('SCAN'):
                    scan = int(line.split('='))
                elif line.startswith("RTINSECONDS"):
                    rt = float(line.split("=")[-1])
                elif line.startswith('END IONS'):
                    self.scan_to_RT_dict[scan] = rt
                    scan = -1
            
    def Close(self):
        pass
        
    def RTInSecondsFromScanNum(self, scanNum):
        return self.scan_to_RT_dict[scanNum]