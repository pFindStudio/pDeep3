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
            
    def GetIonTypeIndices(self, ion_types):
        indices = []
        for ion_type in ion_types:
            for ion_charge in range(1,self.config.max_ion_charge+1):
                idx = self.config.GetIonIndexByIonType(ion_type, ion_charge)
                if idx: indices.append(idx)
        return indices
        
    def GetIntensitiesByIndices(self, pepinfo, indices):
        return self.GetIntensities(pepinfo)[:indices]
        
    def GetIntensitiesByIndices(self, sequence, modification, charge, indices):
        return self.GetIntensities(sequence, modification, charge)[:indices]
            
    def GetIntensitiesByIonTypes(self, pepinfo, ion_types):
        '''
        @param pepinfo. "sequence|modification|precursor_charge".
        @param ion_types. ion_types can be a list containing "b","y","b-ModLoss","y-ModLoss".
        '''
        return self.IntensitiesByIonType(self.GetIntensities(pepinfo), ion_types)
            
    def GetIntensitiesByIonTypes(self, sequence, modification, precursor_charge, ion_types):
        '''
        @param pepinfo. "sequence|modification|precursor_charge".
        @param ion_types. ion_types can be a list containing "b","y","b-ModLoss","y-ModLoss".
        '''
        return self.IntensitiesByIonType(self.GetIntensities(sequence, modification, precursor_charge), ion_types)
        
    def IntensitiesByIonType(self, intensities, ion_types):
        '''
        @param intensities. self.GetIntensities, which is a np.ndarray with shape [len(sequence)-1, 8]. The order for the default 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++].
        @param ion_types. ion_types can be a list containing "b","y","b-ModLoss","y-ModLoss".
        
        @return a predicted n-D intensity np.ndarray for the given n ion_types. Note that the order of the intensities is from N-term to C-term. For example, for 2+ b ions, the intensities are [b1++, b2++, b3++, ..., y(n-1)++]; for 2+ y ions, the intensities are [y(n-1)++, y(n-2)++, ..., y2++, y1++]. If ion_type or ion_charge is not in the given list, return np.zeros.
        '''
        self.GetIonTypeIndices(ion_types)
        if indices: return intensities[:,indices]
        else: return None
    
    def GetRetentionTime(self, pepinfo, modinfo = None, precursor_charge = None):
        if modinfo is not None:
            modinfo = modinfo.strip(";")
            pepinfo = "%s|%s|%d"%(pepinfo, modinfo, precursor_charge)
        if pepinfo not in self.peptide_RT_dict: return None
        else: return float(self.peptide_RT_dict[pepinfo])
        
    def GetIntensities(self, pepinfo, modinfo = None, precursor_charge = None):
        '''
        Get the predicted intensities (np.ndarray with shape=[n-1, 8]) for the given sequence, modinfo, precursor_charge. The order for the default 8 ion_types is [b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+, y-ModLoss++]. If the peptide is not in predictions, return np.zeros.
        '''
        if modinfo is not None:
            modinfo = modinfo.strip(";")
            pepinfo = "%s|%s|%d"%(pepinfo, modinfo, precursor_charge)
        if pepinfo not in self.peptide_intensity_dict: return None
        else: return self.peptide_intensity_dict[pepinfo]
        