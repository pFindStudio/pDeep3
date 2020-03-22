import sqlite3
import numpy as np
import time

from ...utils.mass_calc import PeptideIonCalculator as ion_calc
from ..library_base import LibraryBase
from .tsv import pDeepFormat2PeptideModSeq, PeptideModSeq2pDeepFormat

class OSW(LibraryBase):
    def __init__(self, pDeepParam=None):
        self.sql_conn = None
        
    def Open(self, osw_file):
        self.sql_conn = sqlite3.connect(osw_file)
        self.cursor = self.sql_conn.cursor()
        self.osw_file = osw_file
        
    def Close(self):
        self.sql_conn.close()
        self.osw_file = None
    
    def GetAllPeptides(self):
        start = time.perf_counter()
        self.peptide_dict = {}
        peptide_list = []
        cursor = self.cursor.execute("SELECT pep.MODIFIED_SEQUENCE, pre.CHARGE, fea.EXP_RT FROM score_ms2 ms2 INNER JOIN feature fea ON ms2.FEATURE_ID=fea.ID INNER JOIN precursor pre ON pre.ID=fea.precursor_id INNER JOIN precursor_peptide_mapping prepep ON prepep.precursor_id=fea.precursor_id INNER JOIN peptide pep ON pep.ID=prepep.peptide_id where ms2.rank=1 and ms2.qvalue<0.01")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            charge = int(row[1])
            RT = float(row[2])
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = (row[0], charge, RT, '', -1, '') #items = [PeptideModSeq, CHARGE, RT, raw, scan, protein]
        print("[pDeep Info] reading osw time = %.3fs"%(time.perf_counter() - start))
        return peptide_list

# not support insertion, only create from empty.pqp
class PQP(LibraryBase):
    def __init__(self, pDeepParam = None):
        super(self.__class__, self).__init__(pDeepParam)
        self.sql_conn = None
        self.peptide_dict = {}
        self.peptide_list = []
        self.precursor_peptide_dict = {}
        self.peptide_protein_dict = {}
        
        #CREATE TABLE TRANSITION(ID INT PRIMARY KEY NOT NULL,TRAML_ID TEXT NULL,PRODUCT_MZ REAL NOT NULL,CHARGE INT NULL,TYPE CHAR(1) NULL,ANNOTATION TEXT NULL,ORDINAL INT NULL,DETECTING INT NOT NULL,IDENTIFYING INT NOT NULL,QUANTIFYING INT NOT NULL,LIBRARY_INTENSITY REAL NULL,DECOY INT NOT NULL)
        self.transition_insert_sql = "INSERT INTO TRANSITION(ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ANNOTATION, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES(?, ?, ?, ?, ?, ?, ?, 1, 0, 1, ?, ?)"
        
        # CREATE TABLE PRECURSOR(ID INT PRIMARY KEY NOT NULL,TRAML_ID TEXT NULL,GROUP_LABEL TEXT NULL,PRECURSOR_MZ REAL NOT NULL,CHARGE INT NULL,LIBRARY_INTENSITY REAL NULL,LIBRARY_RT REAL NULL,LIBRARY_DRIFT_TIME REAL NULL,DECOY INT NOT NULL)
        self.precursor_insert_sql = "INSERT INTO PRECURSOR(ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, LIBRARY_RT, LIBRARY_DRIFT_TIME, DECOY) VALUES(?, ?, ?, ?, ?, '', ?, -1.0, ?)"
        
        # CREATE TABLE PEPTIDE(ID INT PRIMARY KEY NOT NULL,UNMODIFIED_SEQUENCE TEXT NOT NULL,MODIFIED_SEQUENCE TEXT NOT NULL,DECOY INT NOT NULL)
        self.peptide_insert_sql = "INSERT INTO PEPTIDE(ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES(?, ?, ?, ?)"
        
        # CREATE TABLE PROTEIN(ID INT PRIMARY KEY NOT NULL,PROTEIN_ACCESSION TEXT NOT NULL,DECOY INT NOT NULL)
        self.protein_insert_sql = "INSERT INTO protein(ID, PROTEIN_ACCESSION, DECOY) VALUES(?, ?, ?)"
        
        # CREATE TABLE PEPTIDE_GENE_MAPPING(PEPTIDE_ID INT NOT NULL,GENE_ID INT NOT NULL)
        self.peptide2gene_sql = "INSERT INTO PEPTIDE_GENE_MAPPING(PEPTIDE_ID, GENE_ID) VALUES(?, 0)"
        
        # CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID INT NOT NULL,PROTEIN_ID INT NOT NULL)
        self.peptide2protein_sql = "INSERT INTO PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID, PROTEIN_ID) VALUES(?, ?)"
        
        # CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INT NOT NULL,PEPTIDE_ID INT NOT NULL)
        self.precursor2peptide_sql = "INSERT INTO PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID, PEPTIDE_ID) VALUES(?, ?)"
        
        # CREATE TABLE TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID INT NOT NULL,PRECURSOR_ID INT NOT NULL)
        self.transition2precursor_sql = "INSERT INTO TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID, PRECURSOR_ID) VALUES(?, ?)"
        
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
        cursor = self.cursor.execute("SELECT pep.MODIFIED_SEQUENCE, pre.CHARGE, pre.LIBRARY_RT FROM precursor pre INNER JOIN precursor_peptide_mapping prepep ON pre.ID=prepep.precursor_id INNER JOIN peptide pep on prepep.peptide_id=pep.ID")
        for row in cursor:
            seq, mod = PeptideModSeq2pDeepFormat(row[0])
            charge = int(row[1])
            RT = float(row[2])*60
            peptide_list.append((seq, mod, charge))
            self.peptide_dict["%s|%s|%d"%(seq,mod,charge)] = [row[0], charge, RT, -1, ""] #items = [PeptideModSeq, CHARGE, RT, scan, protein]
        print("[pDeep Info] reading pqp time = %.3fs"%(time.perf_counter() - start))
        return peptide_list
        
    def _write_one_peptide(self, seq, mod, pre_charge, pepmass, masses, intens, charges, types, sites, RT, protein, transition_count, transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list, decoy=0):
    
        labeled_seq = pDeepFormat2PeptideModSeq(seq, mod)
        
        if decoy: protein = self.decoy_tag+protein
        
        precursor_row = []
        # "INSERT INTO PRECURSOR(ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, LIBRARY_RT, LIBRARY_DRIFT_TIME, DECOY) VALUES(?, ?, ?, ?, ?, NULL, ?, -1.0, ?)"
        precursor_id = len(precursor2peptide)
        label = "%s%d_%s_%d"%("" if not decoy else self.decoy_tag, precursor_id, labeled_seq, pre_charge)
        precursor_row.append(precursor_id)
        precursor_row.append(precursor_id)
        precursor_row.append(label)
        precursor_row.append(pepmass)
        precursor_row.append(pre_charge)
        precursor_row.append(RT/60)
        precursor_row.append(decoy)
        precursor_list.append(precursor_row)
        
        if labeled_seq in peptide_dict:
            precursor2peptide[precursor_id] = peptide_dict[labeled_seq]
            precursor2peptide_list.append((precursor_id, peptide_dict[labeled_seq]))
        else:
            precursor2peptide[precursor_id] = len(peptide_dict)
            precursor2peptide_list.append((precursor_id, len(peptide_dict)))
            peptide_list.append((len(peptide_dict), seq, labeled_seq, decoy))
            
            peptide_dict[labeled_seq] = len(peptide_dict)
            
        if protein in protein_dict:
            peptide2protein[peptide_dict[labeled_seq]] = protein_dict[protein]
            peptide2protein_list.append((peptide_dict[labeled_seq],protein_dict[protein]))
        else:
            peptide2protein[peptide_dict[labeled_seq]] = len(protein_dict)
            peptide2protein_list.append((peptide_dict[labeled_seq],len(protein_dict)))
            protein_list.append((len(protein_dict), protein, decoy))
            
            protein_dict[protein] = len(protein_dict)
            
        for mz, inten, charge, ion_type, site in zip(masses, intens, charges, types, sites):
            # "INSERT INTO TRANSITION(ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ANNOTATION, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES(?, ?, ?, ?, ?, ?, ?, 1, 0, 1, ?, ?)"
            transition_row = []
            transition_row.append(transition_count)
            transition_row.append("%s%d_%s_%d"%("" if not decoy else self.decoy_tag, transition_count, labeled_seq, pre_charge))
            transition_row.append(mz)
            transition_row.append(int(charge))
            transition_row.append(ion_type)
            transition_row.append('%s%d^%d/0'%(ion_type, site, charge) if charge > 1 else '%s%d/0'%(ion_type, site))
            transition_row.append(int(site))
            transition_row.append(float(inten))
            transition_row.append(decoy)
            transition_list.append(transition_row)
            transition2precursor_list.append((transition_count, precursor_id))
            transition_count += 1
            
        return transition_count
        
    def UpdateByPrediction(self, _prediction, peptide_to_protein_dict = {}):
        print("[pDeep Info] updating pqp ...")
        transition_count = 0
        count = 0
        transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list = [], [], [], {}, [], {}, [], {}, [], {}, []
        start = time.perf_counter()
        
        if self.decoy: total = len(_prediction.peptide_intensity_dict)*2
        else: total = len(_prediction.peptide_intensity_dict)
        
        for pepinfo, intensities in _prediction.peptide_intensity_dict.items():
            seq, mod, charge = pepinfo.split("|")
            
            pre_charge = int(charge)
            pepmass, masses, intens, sites, types, charges, decoy_seq, decoy_mod, decoy_masses = self._calc_ions(seq, mod, pre_charge, intensities)
            
            RT = _prediction.GetRetentionTime(pepinfo)
            
            if seq in peptide_to_protein_dict:
                protein = peptide_to_protein_dict[seq]
                protein = str(protein.count("/")+1)+"/"+protein
            elif pepinfo in self.peptide_dict:
                protein = self.peptide_dict[pepinfo][-1]
            else:
                protein = "1/pDeep"
                
            transition_count = self._write_one_peptide(seq, mod, pre_charge, pepmass, masses, intens, charges, types, sites, RT, protein, transition_count, transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list)
            count += 1
            if self.decoy:
                transition_count = self._write_one_peptide(decoy_seq, decoy_mod, pre_charge, pepmass, decoy_masses, intens, charges, types, sites, RT, protein, transition_count, transition_list, transition2precursor_list, precursor_list, precursor2peptide, precursor2peptide_list, peptide_dict, peptide_list, peptide2protein, peptide2protein_list, protein_dict, protein_list, 1)
                count += 1
            if count%10000 == 0:
                print("[PQP UPDATE] {:.1f}%".format(100.0*count/total), end="\r")
                if count%100000 == 0:
                    self.cursor.executemany(self.transition_insert_sql, transition_list)
                    self.cursor.executemany(self.precursor_insert_sql, precursor_list)
                    self.cursor.executemany(self.peptide_insert_sql, peptide_list)
                    self.cursor.executemany(self.protein_insert_sql, protein_list)
                    self.cursor.executemany(self.transition2precursor_sql, transition2precursor_list)
                    self.cursor.executemany(self.precursor2peptide_sql, precursor2peptide_list)
                    self.cursor.executemany(self.peptide2protein_sql, peptide2protein_list)
                    self.sql_conn.commit()
                    transition_list, transition2precursor_list, precursor_list, precursor2peptide_list, peptide_list, peptide2protein_list, protein_list = [], [], [], [], [], [], []
                    
        self.cursor.executemany(self.transition_insert_sql, transition_list)
        self.cursor.executemany(self.precursor_insert_sql, precursor_list)
        self.cursor.executemany(self.peptide_insert_sql, peptide_list)
        self.cursor.executemany(self.protein_insert_sql, protein_list)
        self.cursor.executemany(self.transition2precursor_sql, transition2precursor_list)
        self.cursor.executemany(self.precursor2peptide_sql, precursor2peptide_list)
        self.cursor.executemany(self.peptide2protein_sql, peptide2protein_list)
        self.cursor.executemany(self.peptide2gene_sql, [(peptide_id, ) for peptide_id in peptide_dict.values()])
        self.sql_conn.commit()
        print("[PQP UPDATE] 100%: {}".format(self.pqp_file))
        print("[pDeep Info] updating pqp time = %.3fs"%(time.perf_counter()-start))
        if not self.decoy: print("[pDeep Info] only target transitions were generated!")
        
        
        