powershell run_osw.ps1 -raw_dir e:\DIAData\Specter\HEK_SpikeP100\108ng/ -output_dir e:\DIAData\Specter\HEK_SpikeP100\108ng/DDA -lib e:\DIAData\Specter\HEK_SpikeP100\osw_lib\HEK_P100_DDA_TDA.pqp -irt e:\DIATools\openswath\library\pDeep/cell_RT_proteins_fasta_QE27.tsv -use_window

powershell run_osw.ps1 -raw_dir e:\DIAData\Specter\HEK_SpikeP100\108ng/ -output_dir e:\DIAData\Specter\HEK_SpikeP100\108ng/phl -lib e:\DIATools\openswath\library\phl004_s32_TDA.pqp -irt e:\DIATools\openswath\library\phl004_s32_iRT_cells.tsv -use_window


powershell run_osw.ps1 -export -raw_dir e:\DIAData\Specter\HEK_SpikeP100\108ng/ -output_dir e:\DIAData\Specter\HEK_SpikeP100\108ng/DDA
powershell run_osw.ps1 -export -raw_dir e:\DIAData\Specter\HEK_SpikeP100\108ng/ -output_dir e:\DIAData\Specter\HEK_SpikeP100\108ng/phl