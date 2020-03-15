import sqlite3
import numpy as np
import time

from ...utils.mass_calc import PeptideIonCalculator as ion_calc
from ..library_base import LibraryBase
from .tsv import pDeepFormat2PeptideModSeq, PeptideModSeq2pDeepFormat

# not support insertion, only create from empty.pqp
class PQP(LibraryBase):
    def __init__(self):
        self.sql_conn = None
        self._ion_calc = ion_calc()
        self.peptide_dict = {}
        self.peptide_list = []
        self.precursor_peptide_dict = {}
        self.peptide_protein_dict = {}
        
    def Open(self, pqp_file):
        self.sql_conn = sqlite3.connect(pqp_file)
        self.cursor = self.sql_conn.cursor()
        self.pqp_file = pqp_file
        
    def Close(self):
        self.sql_conn.close()
        self.pqp_file = None
    
    def GetAllPeptides(self):
        start = time.perf_counter()
        self.peptide_dict = {}
        peptide_list = []
        cursor = self.cursor.execute("SELECT c.MODIFIED_SEQUENCE, a.CHARGE, a.LIBRARY_RT, b.precursor_id, b.peptide_id FROM precursor a INNER JOIN precursor_peptide_mapping b ON a.ID=b.precursor_id INNER JOIN peptide c on b.peptide_id=c.ID")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            charge = int(row[1])
            RT = float(row[2])*60
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = [row[0], charge, RT, row[3], row[4]] #items = [PeptideModSeq, CHARGE, RT, precursor_id, peptide_id]
        print("reading pqp time = %.3fs"%(time.perf_counter() - start))
        return peptide_list
        
    def _clear_transition(self):
        self.cursor.excute("DELETE FROM transition")
        self.cursor.excute("DELETE FROM transition_precursor_mapping")
        self.sql_conn.commit()
        
    # CREATE TABLE TRANSITION(ID INT PRIMARY KEY NOT NULL,TRAML_ID TEXT NULL,PRODUCT_MZ REAL NOT NULL,CHARGE INT NULL,TYPE CHAR(1) NULL,ANNOTATION TEXT NULL,ORDINAL INT NULL,DETECTING INT NOT NULL,IDENTIFYING INT NOT NULL,QUANTIFYING INT NOT NULL,LIBRARY_INTENSITY REAL NULL,DECOY INT NOT NULL)
    def _insert_transition(self):
        pass
        
    # peak_selection = "topK"? or "intensity"
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}, peak_selection = "intensity", threshold = 0.05, mass_upper = 2000):
        pass
        
        
        