class pFindMGFReader(object):
    def __init__(self, filename, **kwargs):
        self.scan_to_RT_dict = {}
        with open(filename) as f:
            while True:
                line = f.readline()
                if not line: break
                if line.startswith("TITLE"):
                    scan = int(line.split(".")[-4])
                elif line.startswith("RTINSECONDS"):
                    rt = float(line.split("=")[-1])
                    self.scan_to_RT_dict[scan] = rt
            
    def Close(self):
        pass
        
    def RTInSecondsFromScanNum(self, scanNum):
        return self.scan_to_RT_dict[scanNum]