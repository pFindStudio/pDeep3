#coding=utf-8
'''
Created on 2013.12.13

@author: dell
'''

import numpy as np
from scipy.stats.stats import pearsonr

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure as Figure
import os
import struct

from ..utils.mass_calc import PeptideIonCalculator
from ..pyRawFileReader.RawFileReader import RawFileReader

ioncalc = PeptideIonCalculator()

def get_ax(fig):
    ax = fig.add_subplot(111)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    return ax

class b2bplot(object):
    '''
    classdocs
    '''
    #
    def __init__(self, pdeep_prediction, plot_ion_types, tol = 20, tol_type = "ppm"):
        self.mass_proton = ioncalc.base_mass.mass_proton
        
        self.tol = tol
        self.tol_type = tol_type
        self.mz_bin_size = 0.02
        if self.tol_type.upper() == "DA":
            self.mz_bin_size = self.tol
        self.max_plot_mz = 2100
        self.min_plot_inten = 1e-4
        
        self.ion_color = {'b{}':'#1E1EE5','y{}':'#E51E1E','c{}':'darkgreen','z{}':'#9370db','b{}-ModLoss':'#1E1EE57F','y{}-ModLoss':'#E51E1E7F','b{}-H2O':'#1E1EE5','y{}-H2O':'#E51E1E','b{}-NH3':'#1E1EE5','y{}-NH3':'#E51E1E'}
        self.config = pdeep_prediction.config
        self.ion_types = plot_ion_types
        self.pdeep_prediction = pdeep_prediction
        self.ion_indices, self.used_ion_types = pdeep_prediction.GetIonTypeIndices(self.ion_types)
        
        self.show_plot = True
        self.raw_reader = None
        
    
    def calc_tol(self, mz):
        if self.tol_type.upper() == "DA":
            return self.tol
        else:
            return self.tol * mz / 1000000.0
        
    def open_raw(self, raw_path):
        if self.raw_reader: self.raw_reader.Close()
        self.raw_reader = RawFileReader(raw_path)
        
    def close_raw(self):
        if self.raw_reader: self.raw_reader.Close()
        self.raw_reader = None
        
    def __read_one_spec__(self, raw_scan):
        return self.raw_reader.GetCentroidMassListFromScanNum(raw_scan)
        
    def peak_hashing(self, peaks):
        max_mz, __ = max(peaks, key = lambda x: x[0])
        self.hash_table = [[] for i in range(int(max_mz / self.mz_bin_size)+10)]
        for i in range(len(peaks)):
            self.hash_table[int(peaks[i][0] / self.mz_bin_size)].append(i)
            
    def match_use_hashing(self, ions, charge, peaks):
        matched_peak_idx = [-1] * len(ions)
        matched_peak_inten = [0] * len(ions)
        for i in range(len(ions)):
            mass = ions[i] / charge + self.mass_proton
            if (mass / self.mz_bin_size) > len(self.hash_table) - 5: continue 
            min_tol = self.calc_tol(mass)
            hashed_peak_idlist = self.hash_table[int(mass / self.mz_bin_size)]
            for idx in hashed_peak_idlist:
                if abs(peaks[idx][0] - mass) <= min_tol:
                    matched_peak_idx[i] = idx
                    matched_peak_inten[i] = peaks[idx][1]
                    min_tol = abs(peaks[idx][0] - mass)
            
            # check the bound out of the bin
            hashed_peak_idlist = self.hash_table[int(mass / self.mz_bin_size)-1]
            for idx in hashed_peak_idlist:
                if abs(peaks[idx][0] - mass) <= min_tol:
                    matched_peak_idx[i] = idx
                    matched_peak_inten[i] = peaks[idx][1]
                    min_tol = abs(peaks[idx][0] - mass)
            
            hashed_peak_idlist = self.hash_table[int(mass / self.mz_bin_size)+1]
            for idx in hashed_peak_idlist:
                if abs(peaks[idx][0] - mass) <= min_tol:
                    matched_peak_idx[i] = idx
                    matched_peak_inten[i] = peaks[idx][1]
                    min_tol = abs(peaks[idx][0] - mass)
        
        return matched_peak_idx, matched_peak_inten
    
    def Plot_Predict(self, plotax, ions, predictions):
        valign = 'top'
        vmargin = -0.05
        
        predictions = -predictions*self.max_real_inten
        matched_inten = []
        for charged_ion_type, type_idx in zip(self.used_ion_types, self.ion_indices):
            ion_type = charged_ion_type.strip('+')
            charge = len(charged_ion_type) - len(ion_type)
            if charge > self.spec_charge: continue
            target_ions = ions[ion_type]
            for i in range(len(target_ions)):
                x = target_ions[i] / charge + self.mass_proton
                y = predictions[i, type_idx]
                if x > self.max_plot_mz or y > -self.min_plot_inten or x < 10: #x < 2Da for modloss outside the modsite
                    matched_inten.append(0)
                else:
                    if self.config.ion_terms[ion_type] == 'c': ion_idx = len(target_ions)-i
                    else: ion_idx = i+1
                    if self.show_plot:
                        plotax.plot( [x,x], [0, y], color=self.ion_color[ion_type], lw=2)
                        plotax.text( x, y + vmargin, charged_ion_type.format(str(ion_idx)), rotation = 90, color=self.ion_color[ion_type], horizontalalignment="center",verticalalignment=valign)
                    matched_inten.append(y)
        return matched_inten
    
    def Plot_Real(self, plotax, ions, scan):
        peaks = self.__read_one_spec__(scan)
        xmz = peaks[0,:]
        yint = peaks[1,:]
        yint = yint[xmz <= self.max_plot_mz]
        xmz = xmz[xmz <= self.max_plot_mz]
        max_inten = np.max(yint)
        yint /= max_inten
        peaks = list(zip(xmz, yint))
        self.peak_hashing(peaks)
        
        self.real_max_mz = min(np.max(xmz)+100, self.max_plot_mz)
        
        if self.show_plot: plotax.vlines(xmz, [0]*len(yint), yint, color = 'lightgray')
        
        valign = 'bottom'
        vmargin = 0.05
        
        for charge in range(1, self.spec_charge + 1):
            peak_idx, peak_inten = self.match_use_hashing([self.pepmass], charge, peaks)
            idx = peak_idx[0]
            if idx != -1:
                if self.show_plot: 
                    plotax.text( peaks[idx][0], peaks[idx][1] + vmargin, '{}({}+)'.format('M', charge),  color='gray', horizontalalignment="center",verticalalignment=valign)
        
        matched_inten = []
            # vmargin *= charge
        for charged_ion_type, type_idx in zip(self.used_ion_types, self.ion_indices):
            ion_type = charged_ion_type.strip('+')
            charge = len(charged_ion_type) - len(ion_type)
            if charge > self.spec_charge: continue
            target_ions = ions[ion_type]
            peak_idx, peak_inten = self.match_use_hashing(target_ions, charge, peaks)
            matched_inten.extend(peak_inten)
            for i in range(len(peak_idx)):
                idx = peak_idx[i]
                if idx != -1:
                    if peaks[idx][1] >= self.min_plot_inten:
                        if self.config.ion_terms[ion_type] == 'c': ion_idx = len(target_ions)-i
                        else: ion_idx = i+1
                        if self.show_plot: 
                            plotax.plot( [peaks[idx][0], peaks[idx][0]], [0, peaks[idx][1]], color=self.ion_color[ion_type], lw=2)
                            plotax.text( peaks[idx][0], peaks[idx][1] + vmargin, charged_ion_type.format(str(ion_idx)), rotation = 90, color=self.ion_color[ion_type], horizontalalignment="center",verticalalignment=valign)
        if self.show_plot: plotax.text(20,1.2,'x {:.2e}'.format(max_inten))
        return matched_inten
    
    def plot(self, raw_name, scan, peptide, modinfo, charge):
        self.raw_scan = "{}.{}".format(raw_name, scan)
        if self.raw_reader.LastSpectrumNumber < scan:
            print('no spec {} in spectrum file'.format(self.raw_scan)) 
            return None, None
        print('{}-{}-{}+ <-> {}'.format(peptide, modinfo, charge, self.raw_scan))
        
        modinfo = modinfo.strip(";")
        mod_cumsum, modloss_list, modname_list = ioncalc.calc_modification_mass(peptide, modinfo)
        
        ions = {}
        bions, pepmass = ioncalc.calc_b_ions_and_pepmass(peptide, mod_cumsum)
        print("parent m/z (%d+) = %.6f" %(charge, pepmass/charge + self.mass_proton))
        if 'b{}' in self.config.ion_types: ions['b{}'] = bions
        if 'y{}' in self.config.ion_types: ions['y{}'] = ioncalc.calc_y_from_b(bions, pepmass)
        if 'c{}' in self.config.ion_types: ions['c{}'] = ioncalc.calc_c_from_b(bions)
        if 'z{}' in self.config.ion_types: ions['z{}'] = ioncalc.calc_z_from_b(bions, pepmass)
        if 'b{}-ModLoss' in self.config.ion_types: ions['b{}-ModLoss'] = ioncalc.calc_Nterm_modloss(bions, modloss_list, modname_list)
        if 'y{}-ModLoss' in self.config.ion_types: ions['y{}-ModLoss'] = ioncalc.calc_Cterm_modloss(ions['y{}'], modloss_list, modname_list)
        
        if self.show_plot: 
            fig = Figure(figsize=(12,8), dpi=80)
            ax = get_ax(fig)
        else:
            fig = None
            ax = None
        self.spec_charge = charge
        self.pepmass = pepmass
        
        matched_inten1 = self.Plot_Real(ax, ions, scan)
        self.max_real_inten = max(matched_inten1)
        
        predictions = self.pdeep_prediction.GetIntensities(peptide, modinfo, charge)
        if predictions is None: return None, None
        matched_inten2 = self.Plot_Predict(ax, ions, predictions)
        
        if len(matched_inten1) < len(matched_inten2):
            matched_inten2 = matched_inten2[:len(matched_inten1)]
            print('[Warning] ion number is not equal')
        elif len(matched_inten1) > len(matched_inten2):
            matched_inten1 = matched_inten1[:len(matched_inten2)]
            print('[Warning] ion number is not equal')
        PCC = pearsonr(np.array(matched_inten1), -np.array(matched_inten2))[0]
        
        if self.show_plot:
            modinfo = modinfo.strip(";")
            plot_pep = peptide
            if modinfo:
                modlist = []
                for site_mod in modinfo.split(";"):
                    site, mod = site_mod.split(",")
                    modlist.append((int(site), mod))
                modlist.sort(reverse=True)
                
                for site, mod in modlist:
                    if site == 0: site = 1
                    elif site > len(peptide): site = len(peptide)
                    plot_pep = "%s[%s]%s"%(plot_pep[:site],mod[:3],plot_pep[site:])
                
            ax.text(200, 1.3, '{} ({}+), R = {:.2f}'.format(plot_pep,self.spec_charge,PCC), fontsize = 14)
            
            ax.set_xlim(left = 0, right = self.real_max_mz)
            ax.set_ylim(bottom = -1.2, top=1.4)
            ax.hlines([0], [0], [self.real_max_mz])
            ax.set_xlabel('m/z')
            ax.set_ylabel('Relative Abundance')
            ylabels = ax.get_yticks().tolist()
            # ylabels = ["{:.0f}%".format(abs(float(label)*100)) if float(label) >= 0 else '' for label in ylabels]
            ylabels = ["{:.0f}%".format(abs(float(label)*100)) for label in ylabels]
            #ylabels = ['' for label in ylabels]
            ax.set_yticklabels(ylabels)
            
            print('R = {:.2f}'.format(PCC))
        return (fig, PCC)
        
    def batch_save(self, spec_file, pep_list, save_dir):
        self.open_raw(spec_file)
        for pepinfo in pep_list:
            raw, scan, peptide, mod, charge = pepinfo[:5]
            fig, pcc = self.plot(raw, int(scan), peptide, mod, int(charge))
            if fig is None: plt.close()
            else:
                plt.tight_layout()
                # mng = plt.get_current_fig_manager()
                #mng.window.showMaximized()
                #mng.resize(*mng.window.maxsize())
                   
                
                plt.savefig(os.path.join(save_dir, "%s.%s-R=%.3f-%s.png"%(raw, scan, pcc, peptide)),format="png", dpi=120)
                plt.close()
        self.close_raw()
    
    def show(self, save_as = None):
        if self.show_plot and save_as is None: plt.show()
        elif save_as is not None: plt.savefig(save_as, dpi=120)
        
if __name__ == "__main__":
    from ..cmd.tune_and_predict import get_prediction
    matplotlib.use('TkAgg')
    raw_path = r"e:\DIAData\Specter\HEK_SpikeP100\DDA_data\CS20170922_SV_HEK_SpikeP100_108ng_DDA.raw"
    psm_path = r"e:\DIAData\Specter\HEK_SpikeP100\DDA_data\pFind3\result\pFind-Filtered.spectra.psm.txt"
    pdeep_prediction = get_prediction(psm_path)
    
    bbplot = b2bplot(pdeep_prediction, ['b{}', 'y{}', 'b{}-ModLoss', 'y{}-ModLoss'])
    bbplot.open_raw(raw_path)
    
    raw, scan, peptide, modinfo, charge = 'CS20170922_SV_HEK_SpikeP100_108ng_DDA	39575	GHVFEESQVAGTPMFVVK	14,Oxidation[M]	3'.split('\t')[:5]
    bbplot.plot(raw, int(scan), peptide, modinfo, int(charge))
    bbplot.show()#save_as = r'e:\DIAData\Specter\HEK_SpikeP100\DDA_data\bbplot1.png')
    
    raw, scan, peptide, modinfo, charge = 'CS20170922_SV_HEK_SpikeP100_108ng_DDA	39018	VLVEPDAGAGVAVMK		2	2302.39572'.split('\t')[:5]
    bbplot.plot(raw, int(scan), peptide, modinfo, int(charge))
    bbplot.show()#save_as = r'e:\DIAData\Specter\HEK_SpikeP100\DDA_data\bbplot2.png')
    
    bbplot.close_raw()
    
    
