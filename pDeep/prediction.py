from .bucket import peptide_as_key

class pDeepPrediction:
    def __init__(self, config, peptide_buckets, predict_buckets, predict_RT_buckets = None):
        '''
        @param config. from pDeep.config.pDeep_config
        @param peptide_buckets. pDeep input buckets
        @param predict_buckets. pDeep predict buckets
        '''
        self.config = config
        self.peptide_intensity_dict = peptide_as_key(peptide_buckets, predict_buckets)
        if predict_RT_buckets:
            self.peptide_RT_dict = peptide_as_key(peptide_buckets, predict_RT_buckets)
        else: 
            self.peptide_RT_dict = {}
        
    def GetIntensitiesByIonType(self, intensities, ion_type, ion_charge):
        '''
        @param intensities. self.GetIntensities, which is a np.ndarray with shape [len(sequence)-1, 8]. The order for the default 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++].
        @param ion_type. ion_type can be "b","y","b-ModLoss","y-ModLoss".
        @param ion_charge. ion_charge can be 1 and 2.
        
        @return a predicted 1-D intensity np.ndarray for the given ion_type and ion_charge. Note that the order of the intensities is from N-term to C-term. For example, for 2+ b ions, the intensities are [b1++, b2++, b3++, ..., y(n-1)++]; for 2+ y ions, the intensities are [y(n-1)++, y(n-2)++, ..., y2++, y1++]. If ion_type or ion_charge is not in the given list, return np.zeros.
        '''
        idx = self.config.GetIonIndexByIonType(ion_type, ion_charge)
        if idx is None: return None
        else: return intensities[:,idx] 
    
    def GetRetentionTime(self, pepinfo, modinfo = None, precursor_charge = None):
        if precursor_charge is not None:
            pepinfo = "%s|%s|%d"%(pepinfo, modinfo, precursor_charge)
        if pepinfo not in self.peptide_RT_dict: return None
        else: return float(self.peptide_RT_dict[pepinfo])
        
    def GetIntensities(self, pepinfo, modinfo = None, precursor_charge = None):
        '''
        Get the predicted intensities (np.ndarray with shape=[n-1, 8]) for the given sequence, modinfo, precursor_charge. The order for the default 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++]. If the peptide is not in predictions, return np.zeros.
        '''
        if precursor_charge is not None:
            pepinfo = "%s|%s|%d"%(pepinfo, modinfo, precursor_charge)
        if pepinfo not in self.peptide_intensity_dict: return None
        else: return self.peptide_intensity_dict[pepinfo]
        