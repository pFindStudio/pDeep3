class ReaderBase(object):
    def __init__(self):
        self.peptide_dict = {}
        self._file = None
        self._filename = None
    def Open(self, file):
        self._file = open(file)
        self._filename = file
    def Close(self):
        if self._file is not None:
            self._file.close()
            self._file = None
            self._filename = None
    def GetAllPeptides(self):
        pass