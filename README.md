# pDeep3

## System requirements
python >= 3.5

tensorflow == 1.13.1 (tensorflow == 1.x, where x >= 13)

.NET Framework == 4.5.2 (or higher? to execute pLabel)

## Preparing fine-tuning data using pLabel
***pLabel now can directly access raw files by using Thermo RawFileReader.***

Run pLabel:
```
cd pLabel
pLabel.exe pLabel-sample.cfg
```
or in linux, run .NET applications using mono:
```
mono pLabel.exe pLabel-sample.cfg
```

Example of pLabel-sample.cfg:
```
psm_type = none
mode = pDeep
num_psm_file = 1
psm_file1 = D:\plabel\psm_sample.txt
ms2_type = raw
num_ms2_file = 2
ms2_file1 = D:\plabel\raw_sample1.raw
ms2_file2 = D:\plabel\raw_sample2.raw
output_folder = D:\plabel\output
NH3_loss = true
H2O_loss = true
Mod_loss = true
num_ion_type = 2
iontype1 = b|N_term|0
iontype2 = y|C_term|0
num_new_aa = 0
```
Example of input file psm_sample.txt:
```
raw_name	scan	peptide	modinfo	charge
raw_sample1	7932	TCEATHKTSTSPIVKSF	2,Carbamidomethyl[C];	2
raw_sample1	13419	KIDGMERQDGVLNSW		3
raw_sample2	10709	TCEATHKTSTSPIVKSF	16,Phospho[S];	2
```
***Check the names of modifications in pDeep/config/modification.py.***

Users now can prepare a engine-independent input psm file (psm_sample.txt) for pLabel.

After executing "pLabel.exe pLabel-sample.cfg", it will generate two result files at specified output_folder (D:\plabel\output in the example) named after the two raw files: raw_sample1.psmlabel and raw_sample2.psmlabel, these two files can be used to fine-tune/train and test the pDeep model.

## fine-tune and predict
Usage:
```
python -m pDeep.cmd.tune_and_predict tmp/predict/pDeep-tune.cfg
```
***There will be quite a few warnings when importing tensorflow, just ignore them.***

Example of tmp/predict/pDeep-tune.cfg:
```
model = tmp/model/pretrain-180921-modloss-mod8D.ckpt
threads = 4

###### predict ######
mod_no_check = Carbamidomethyl[C]
mod_check = Oxidation[M],Phospho[Y],Phospho[S],Phospho[T]
min_mod_check = 0
max_mod_check = 3
# format: peptide filename | instrument | NCE
predict_input = tmp/predict/peptide.txt | QE | 28


###### Data for fine-tuning, no tuning if it is empty. Files are seperated by '|' ######
tune_psmlabels = tmp\data\Olsen-Chymo-QE-28\raw0\20150708_QE3_UPLC8_DBJ_QC_HELA_39frac_Chymotrypsin_19.psmlabel
n_tune_per_psmlabel = 100

# Data for testing, no testing if it is empty. Files are seperated by '|'
test_psmlabels = tmp\data\Olsen-Chymo-QE-28\raw1\20150708_QE3_UPLC8_DBJ_QC_HELA_39frac_Chymotrypsin_23_2.psmlabel | tmp\data\Olsen-Chymo-QE-28\raw2\20150708_QE3_UPLC8_DBJ_QC_HELA_39frac_Chymotrypsin_24.psmlabel
```

See pDeep.cmd.tune_and_predict.py for details.

***Note that pDeep/cmd/tune_and_predict.py can also be imported and called by other python scripts.***
