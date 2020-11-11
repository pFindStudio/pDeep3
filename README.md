# pDeep3

## installation
```pip install .```
or
```pip install -e .```

## System requirements
python >= 3.5

tensorflow == 1.13.1 (tensorflow >= 1.13.0)

.NET Framework == 4.5.2 (or higher? to execute psmLabel)

## Tuning and predicting using get_prediction
***Now psmLabel, tuning and predicting are integrated in pDeep.cmd.tune_and_predict.py***
```
pDeep.cmd.tune_and_predict.get_prediction(input_peptides, tune_psm = None, raw = None, instrument = 'QE', ce = 27)

@param input_peptides, could be a peptide list [(sequence1, mod1, charge1), (seq2, mod2, charge2), ...] to be predicted, or a file containing tab seperated head "peptide, modinfo, charge, protein".
@param tune_psm evidence.txt (MaxQuant), .spectra (pFind) or *.psm.txt/*.txt (with tab seperated heads "raw_name, scan, peptide, modinfo, charge, RTInSeconds") file for tuning pDeep and pDeepRT. If it is None, the model will not be tuned (default None).
@param raw, raw file for tuning pDeep and pDeepRT (default None).
@param instrument, instrument type for prediction (default "QE").
@param ce, collision energy for prediction (default 27).
@return prediction, pDeep.prediction.pDeepPrediction object, 
# example: 
# ion_types = ['b','y','b-ModLoss','y-ModLoss'] or ['b','y']
# ion_indices, used_ion_types = prediction.GetIonTypeIndices(ion_types)
# print(used_ion_types)
# intensities = GetIntensitiesByIndices('ACDMNLK', '2,Carbamidomethyl[C];4,Oxidation[M]', 3, ion_indices)
```

Example of input file input_peptides file:
```
raw_name	scan	peptide	modinfo	charge	RTinSeconds
raw_sample1	7932	TCEATHKTSTSPIVKSF	2,Carbamidomethyl[C]	2	666.4
raw_sample1	13419	KIDGMERQDGVLNSW		3	1145.2
raw_sample1	13440	KIDGMERCDGVLNSW	5,Oxidation[M];8,Carbamidomethyl[C]	2	1147.0
raw_sample2	10709	TCEATHKTSTSPIVKSF	16,Phospho[S]	2	901.3
```

Run:
```
from pDeep.cmd.tune_and_predict import get_prediction

input_peptides = [('ACDMNLK', '2,Carbamidomethyl[C];4,Oxidation[M]', 3)]
ion_types = ['b','y']
prediction = get_prediction(input_peptides) # no fine-tuning
ion_indices, used_ion_types = prediction.GetIonTypeIndices(ion_types)
print(used_ion_types) 
# output: ['b+','b++','y+','y++']
print(prediction.GetIntensitiesByIndices(*input_peptides[0], ion_indices))
# output: 
[
[ 0.0000000e+00  0.0000000e+00  0.0000000e+00  0.0000000e+00]
[ 3.5514534e-01  3.9369406e-12  5.2654832e-03  0.0000000e+00]
[ 4.4627541e-01  1.0988828e-09  2.0606143e-02  0.0000000e+00]
[ 1.2166708e-02 -0.0000000e+00  1.4021127e-01  0.0000000e+00]
[ 2.1781624e-08 -0.0000000e+00  2.6448551e-01  8.2786764e-09]
[ 0.0000000e+00 -0.0000000e+00  9.9956471e-01  1.1667855e-06]
]
```

Fine-tuning using prepared psm file:
```
prediction = get_prediction(input_peptides, tune_psm=r"e:\DIAData\PECAN\tune.psm.txt", raw=r"e:\DIAData\PECAN\20141010_DIA_20x5mz_700to800.raw")
```

Or fine-tuning using pFind3's pFind.spectra or pFind-Filtered.spectra:
```
prediction = get_prediction(input_peptides, tune_psm=r"e:\DIAData\PECAN\pFind.spectra", raw=r"e:\DIAData\PECAN\20141010_DIA_20x5mz_700to800.raw")
```

Or fine-tuning using MaxQuant's evidence.txt:
```
prediction = get_prediction(input_peptides, tune_psm=r"e:\DDATools\MaxQuant_1.6.12.0\test_data\combined\txt\evidence.txt", raw=r"e:\DDATools\MaxQuant_1.6.12.0\test_data\20141010_DIA_20x5mz_700to800.raw")
```


## Run them seperately:

## Preparing fine-tuning data using psmLabel
***psmLabel now can directly access raw files by using Thermo RawFileReader.***

Run psmLabel:
```
cd psmLabel
psmLabel.exe psmLabel-sample.cfg
```
or in linux, run .NET applications using mono:
```
mono psmLabel.exe psmLabel-sample.cfg
```

Example of psmLabel-sample.cfg:
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
raw_sample1	7932	TCEATHKTSTSPIVKSF	2,Carbamidomethyl[C]	2
raw_sample1	13419	KIDGMERQDGVLNSW		3
raw_sample1	13440	KIDGMERCDGVLNSW	5,Oxidation[M];8,Carbamidomethyl[C]	2
raw_sample2	10709	TCEATHKTSTSPIVKSF	16,Phospho[S]	2
```
***Check the names of modifications in pDeep/config/modification.py.***

Users now can prepare a engine-independent input psm file (psm_sample.txt) for psmLabel.

After executing "psmLabel.exe psmLabel-sample.cfg", it will generate two result files at specified output_folder (D:\plabel\output in the example) named after the two raw files: raw_sample1.psmlabel and raw_sample2.psmlabel, these two files can be used to fine-tune/train and test the pDeep model.

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
```
import pDeep.cmd.tune_and_predict
pdeep_prediction = tune_and_predict.run("tmp/predict/pDeep-tune.cfg")
for peptide, intensities in pdeep_prediction.peptide_prediction_dict.items():
    print("b+ ions of %s ="%peptide, pdeep_prediction.GetIntensitiesByIonType(intensities, "b", 1)) # get b+ ions
for peptide, intensities in pdeep_prediction.peptide_prediction_dict.items():
    print("y+ ions of %s ="%peptide, pdeep_prediction.GetIntensitiesByIonType(intensities, "y", 1)) # get y+ ions
```

## Generating spectral libraries
Example:
```
python -m pDeep.cmd.generate_predicted_speclib --input xxx.fasta --target_proteins Q1234,Q6789 --output library.pqp --varmod Oxidation[M],Phospho[S] --instrument QE --ce 28 --min_intensity 0.1 --least_n_peaks 6
```
"--min_intensity 0.1 --least_n_peaks 6" means that if there are less than 6 peaks larger than 0.1, top-6 peaks will be kept, otherwise all peaks larger than 0.1 will be kept.
"--output library.pqp" for OpenSWATH, "--output library.dlib" for EncyclopeDIA, "--output library.csv" for Spectronaut.

Run
```
python -m pDeep.cmd.generate_predicted_speclib --help
```
to see detailed usage information.

## Using EThcD model
Example:
```
python -m pDeep.cmd.generate_predicted_speclib --input tmp/predict/peptide.txt --output xxxx.msp --model EThcD --ion_type b,y,c,z --instrument Lumos --ce 28 --min_intensity 0.0001
```
Here, for .msp file, all non-zero ions should be generated, hence "--min_intensity 0.0001".
***Note that EThcD model was only trained by ProteomeTools, so the instrument and ce parameters must be: "--instrument Lumos --ce 28" even you use Fusion or other instruments (28 is the NCE of HCD).***

## Supporting TensorFlow2
tf2 models have been trained for HCD and EThcD, and pDeep will automatically load the tf2 model by checking tf.__version__. However, there is a bug in fine-tuning when using tf2 (v2.1.0), it seems that it has been fixed in the nightly built version (see https://github.com/tensorflow/tensorflow/issues/34211).
