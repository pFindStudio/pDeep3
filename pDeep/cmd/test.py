import os
import sys
import time

import pDeep.config.pDeep_config as fconfig
import pDeep.evaluate as evaluate
import pDeep.load_data as load_data
import pDeep.similarity_calc as sim_calc
from pDeep.bucket import count_buckets

model_folder = 'tmp/model'
# model_name = 'pretrain-180921-modloss-mod8D.ckpt'
model_name = sys.argv[1]

argd = {}
if len(sys.argv) > 2:
    for i in range(2, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i + 1]
n = int(argd['-n'])

ion_types = ['b{}', 'y{}', 'b{}-ModLoss', 'y{}-ModLoss']

if model_name.endswith(".ckpt"):
    import pDeep.model_tf as model
else:
    print("Unknown pDeep model!!")
    sys.exit(-1)

mod_type = "CommonMod"

unmod_config = fconfig.HCD_Config()
unmod_config.time_step = 100
unmod_config.SetIonTypes(ion_types)
unmod_config.min_var_mod_num = 0
unmod_config.max_var_mod_num = 0

mod_config = fconfig.HCD_CommonMod_Config()
mod_config.time_step = 100
mod_config.SetIonTypes(ion_types)
mod_config.min_var_mod_num = 1
mod_config.max_var_mod_num = 2

pho_config = fconfig.HCD_pho_Config()
pho_config.time_step = 100
pho_config.SetIonTypes(ion_types)
pho_config.min_var_mod_num = 1
pho_config.max_var_mod_num = 3

QEHO = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Olsen-CellSys-2017/Hela-QE-28/plabel"
QE293O = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Olsen-CellSys-2017/293T-QE-28/plabel"
PT25_ts = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/zengwenfeng-ProteomeTools/plabel/HCD25/test"
PT30_ts = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/zengwenfeng-ProteomeTools/plabel/HCD30/test"
PT35_ts = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/zengwenfeng-ProteomeTools/plabel/HCD35/test"
VlsMilkK = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Kuster-Human-Nature-2014/milk-Velos-40/plabel"
EltStmK = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Kuster-Human-Nature-2014/stomach-Elite-30/plabel"
VlsFGP = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Pandey-Human-Nature-2014/Fetal_Gut_Gel_Velos_72_CE35/plabel"
VlsFBP = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Pandey-Human-Nature-2014/Fetal_Brain_Gel_Velos_16_CE39/plabel"
# VlsPhoSyn = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Mann-PhoSyn-NBT-2011-Velos-40/plabel"
QEHchyO = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Olsen-CellSys-2017/HelaChymo-QE-28/plabel"
QEHgluO = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Olsen-CellSys-2017/HelaGluC-QE-28/plabel"
QEHlysO = "/home/pfind/Documents/pDeep/pDeepMS2/datasets/Olsen-CellSys-2017/HelaLysC-QE-28/plabel"

strQE = "QE"
strVelos = "Velos"
strElite = "Elite"
strLumos = "Lumos"
strFusion = "Fusion"

# strUnknown = "XXXX"
# strQE = strUnknown
# strVelos = strUnknown
# strElite = strUnknown
# strLumos = strUnknown
# strFusion = strUnknown

pdeep = model.pDeepModel(unmod_config)
pdeep.batch_size = 1024 * 4

plot_folder = os.path.join(model_folder, 'log/plots/%s-%s' % (model_name, mod_type))
try:
    os.makedirs(plot_folder)
except:
    pass

pdeep.LoadModel(model_file=os.path.join(model_folder, model_name))

with open(os.path.join(model_folder, 'log/test_%s-%s.txt' % (model_name, mod_type)), 'w') as log_out:
    def test(folder, ce, ins, n, saveplot, phos=False):
        print('###################### Begin Unmod ######################', file=log_out)
        print("[D] " + folder, file=log_out)
        print("[T] Unmod PSMs:", file=log_out)
        buckets = load_data.load_folder_as_buckets(folder, unmod_config, nce=ce, instrument=ins, max_n_samples=n)
        print("[C] " + str(count_buckets(buckets)), file=log_out)
        output_buckets = pdeep.Predict(buckets)
        pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, buckets)
        sim_names = ['PCC', 'COS', 'SPC', 'KDT', 'SA']
        print("[A] " + str(evaluate.cum_plot([pcc, cos, spc, kdt, SA], sim_names, evaluate.thres_list,
                                             saveplot=os.path.join(plot_folder, saveplot + '.eps'),
                                             print_file=log_out)), file=log_out)
        print('####################### End Unmod #######################', file=log_out)
        print("", file=log_out)

        print('####################### Begin Mod #######################', file=log_out)
        print("[D] " + folder, file=log_out)
        if phos:
            print("[T] Phos PSMs:", file=log_out)
            config = pho_config
            mod = '-pho'
        else:
            print("[T] %s PSMs:" % mod_type, file=log_out)
            config = mod_config
            mod = '-mod'
        buckets = load_data.load_folder_as_buckets(folder, config, nce=ce, instrument=ins, max_n_samples=n)
        print("[C] " + str(count_buckets(buckets)), file=log_out)
        output_buckets = pdeep.Predict(buckets)
        pcc, cos, spc, kdt, SA = sim_calc.CompareRNNPredict_buckets(output_buckets, buckets)
        sim_names = ['PCC', 'COS', 'SPC', 'KDT', 'SA']
        print("[A] " + str(evaluate.cum_plot([pcc, cos, spc, kdt, SA], sim_names, evaluate.thres_list,
                                             saveplot=os.path.join(plot_folder, saveplot + mod + '.eps'),
                                             print_file=log_out)), file=log_out)
        print('######################## End Mod ########################', file=log_out)
        print("\n", file=log_out)


    start_time = time.perf_counter()

    ################# start one folder ##############################
    test_folder, ce, ins = QEHO, 28, strQE
    test(test_folder, ce, ins, n, "QE-H-O", phos=True)
    ################# end one folder ################################

    ################# start one folder ##############################
    test_folder, ce, ins = QE293O, 28, strQE
    test(test_folder, ce, ins, n, "QE-293-O", phos=True)
    ################# end one folder ################################

    ################# start one folder ##############################
    test_folder, ce, ins = VlsMilkK, 40, strVelos
    test(test_folder, ce, ins, n, "Velos-Milk-Ku", phos=False)
    ################# end one folder ################################

    ################# start one folder ##############################
    test_folder, ce, ins = EltStmK, 30, strElite
    test(test_folder, ce, ins, n, "Elite-Stm-Ku", phos=False)
    ################# end one folder ################################

    ################# start one folder ##############################
    test_folder, ce, ins = VlsFGP, 35, strVelos
    test(test_folder, ce, ins, n, "Velos-FtlGut-P", phos=False)
    ################# end one folder ################################

    ################# start one folder ##############################
    test_folder, ce, ins = VlsFBP, 39, strVelos
    test(test_folder, ce, ins, n, "Velos-FtlBrain-P", phos=False)
    ################# end one folder ################################

    ################# start one folder ##############################
    # test_folder,ce,ins = VlsPhoSyn,40,strVelos
    # test(test_folder, ce, ins, n, "Velos-PhosSyn", phos = True)
    ################# end one folder ################################

    ################# start one folder ##############################
    test_folder, ce, ins = QEHchyO, 28, strQE
    test(test_folder, ce, ins, n, "QE-Hchy-O", phos=False)
    ################# end one folder ################################

    ################# start one folder ##############################
    test_folder, ce, ins = QEHgluO, 28, strQE
    test(test_folder, ce, ins, n, "QE-Hglu-O", phos=False)
    ################# end one folder ################################

    ################# start one folder ##############################
    test_folder, ce, ins = QEHlysO, 28, strQE
    test(test_folder, ce, ins, n, "QE-Hlys-O", phos=False)
    ################# end one folder ################################

    ################### start one folder ##############################
    test_folder, ce, ins = PT25_ts, 25, strLumos
    test(test_folder, ce, ins, n, "PT25", phos=False)
    ################### end one folder ################################

    ################### start one folder ##############################
    test_folder, ce, ins = PT30_ts, 30, strLumos
    test(test_folder, ce, ins, n, "PT30", phos=False)
    ################### end one folder ################################

    ################### start one folder ##############################
    test_folder, ce, ins = PT35_ts, 35, strLumos
    test(test_folder, ce, ins, n, "PT35", phos=False)
    ################### end one folder ################################

    end_time = time.perf_counter()

    print("time = {:.3f}s".format(end_time - start_time))
