Set-Location -Path "e:/DIATools/openswath/OpenMS-2.4.0-nightly-2019-08-09/bin/"


# $output_dir="e:\DIAData\PECAN\5mz\test_spikein"
# $global:lib_TDA="e:/DIATools/openswath/library/phl004_s32_TDA.pqp"
# $global:irt="e:/DIATools/openswath/library/phl004_s32_iRT_cells.tsv"

$raw_dir="e:\DIAData\Specter\HEK_SpikeP100\108ng"
$output_dir="e:\DIAData\Specter\HEK_SpikeP100\108ng\osw"
$lib="e:\DIAData\Specter\HEK_SpikeP100\osw_lib\phl_spikein.tsv"
$global:lib_TDA="e:\DIAData\Specter\HEK_SpikeP100\osw_lib\phl_spikein.pqp"
$global:irt="e:/DIATools/openswath/library/pDeep/cell_RT_proteins_fasta_QE27.tsv"

# $output_dir="e:\DIAData\PECAN\plasma\peptide-list-test-reverse"
# $lib="e:/DIATools/openswath/library/pDeep/peptide_list_QE27.tsv"
# $global:lib_TDA="e:/DIATools/openswath/library/pDeep/peptide_list_QE27_reverse.pqp"
# $global:irt="e:/DIATools/openswath/library/pDeep/blood_RT_proteins_fasta_QE27.tsv"


$temp="E:/temp"
$ms1_tol_type='ppm' # or 'Th'
$ms1_tol=10
$ms2_tol_type='ppm' # or 'Th'
$ms2_tol=20
$batch=5000
$thread=6
$overlapDIA='false'
$use_window=1
$global:win="e:\DIAData\Specter\HEK_SpikeP100\DIAwindow-left.txt"
$global:win_left="e:\DIAData\Specter\HEK_SpikeP100\DIAwindow-left.txt"
$global:win_right="e:\DIAData\Specter\HEK_SpikeP100\DIAwindow-right.txt"

function Compare-File-Date($f1, $f2)
{
    return ([datetime](Get-ItemProperty -Path $f1 -Name LastWriteTime).lastwritetime).ToOADate() -lt ([datetime](Get-ItemProperty -Path $f2 -Name LastWriteTime).lastwritetime).ToOADate()
}

New-Item -Path $output_dir -ItemType Directory
if (!(Test-Path $global:lib_TDA -PathType leaf))
{
    Write-Host "========== Converting to PQP =========="
    $pqp=-join($lib.substring(0, $lib.Length-3), 'pqp')
    TargetedFileConverter -in $lib -out $pqp -threads 4
    Write-Host "========== Generating decoy =========="
    OpenSWATHDecoyGenerator -in $pqp -out $global:lib_TDA -method reverse -threads 4
}
function run_one($raw, $out)
{
    Write-Host $raw
    OpenSwathWorkflow -in $raw -tr $global:lib_TDA -sort_swath_maps -readOptions cache -tempDirectory $temp -batchSize $batch -out_osw $out -threads $thread -tr_irt $global:irt -rt_extraction_window 600 -mz_extraction_window_unit $ms2_tol_type -mz_extraction_window $ms2_tol -mz_extraction_window_ms1_unit $ms1_tol_type -mz_extraction_window_ms1 $ms1_tol -matching_window_only $overlapDIA -min_coverage 0.01 -min_rsq 0.95
    # OpenSwathWorkflow '-in' $raw -tr $global:lib_TDA -sort_swath_maps -readOptions cache -tempDirectory E:/Temp -batchSize 10000 -swath_windows_file $global:win -tr_irt $global:irt -out_osw $out -threads 4 -use_ms1_traces
    
    # no matter -use_ms1_traces or not, precursor ion will be matched in ms1? But if -use_ms1_traces, ms1 features will be used in scoring
}

function run_one_window($raw, $out)
{
    Write-Host $raw
    OpenSwathWorkflow -in $raw -tr $global:lib_TDA -sort_swath_maps -readOptions cache -tempDirectory $temp -batchSize $batch -out_osw $out -threads $thread -tr_irt $global:irt -rt_extraction_window 600 -mz_extraction_window_unit $ms2_tol_type -mz_extraction_window $ms2_tol -mz_extraction_window_ms1_unit $ms1_tol_type -mz_extraction_window_ms1 $ms1_tol -matching_window_only $overlapDIA -min_coverage 0.01 -min_rsq 0.95 -swath_windows_file $global:win -force
    # OpenSwathWorkflow '-in' $raw -tr $global:lib_TDA -sort_swath_maps -readOptions cache -tempDirectory E:/Temp -batchSize 10000 -swath_windows_file $global:win -tr_irt $global:irt -out_osw $out -threads 4 -use_ms1_traces
    
    # no matter -use_ms1_traces or not, precursor ion will be matched in ms1? But if -use_ms1_traces, ms1 features will be used in scoring
}

function run_pyprophet($dir, $subsample=0)
{   
    Set-Location -Path $dir
    $tmpl=-join("--template=", $global:lib_TDA)
    $osw_files = Get-ChildItem -Path . -Recurse -Include raw*.osw
    $sub_files = @()
    
    Foreach ($run in $osw_files)
    {
        $bak=-join($run,'.bak')
        Copy-Item -Path $run $bak
    }
    
    if (($subsample -eq 0) -or ($subsample -eq 1))
    {
        pyprophet merge $tmpl --out=score.model $osw_files
    }
    else
    {
        Write-Host "========== Sub-sampling =========="
        Write-Host "Sub-sample ratio =", $subsample
        Foreach ($run in $osw_files)
        {
            pyprophet subsample --in=$run --out=${run}.sub --subsample_ratio=$subsample
            $sub_files += -join($run, '.sub')
        }
        pyprophet merge $tmpl --out=score.model $sub_files
        Write-Host $sub_files
    }
    
    Write-Host "========== Scoring =========="
    pyprophet score --in=score.model --level=ms1ms2
    Foreach ($run in $osw_files)
    {
        pyprophet score --in=$run --apply_weights=score.model --level=ms1ms2
    }
    
    Write-Host "========== Reducing =========="
    $reduced_files = @()
    Foreach ($run in $osw_files)
    {
        $run_reduced = -join(${run},".reduce") # generates .oswr files
        pyprophet reduce --in=$run --out=$run_reduced
        $reduced_files += $run_reduced
    }
        
    Write-Host "========== Estimating FDR =========="
    pyprophet merge --template=score.model --out=global_fdr.model $reduced_files
    pyprophet peptide --context=global --in=global_fdr.model
    pyprophet protein --context=global --in=global_fdr.model
    
    Foreach ($run in $osw_files)
    {
        pyprophet backpropagate --in=$run --apply_scores=global_fdr.model
    }
}

$raw_files = Get-ChildItem -Path $raw_dir -Recurse -Include *.mzML
Write-Host $raw_files
$i = 1
Foreach ($raw in $raw_files)
{
    $out_path=Join-Path -Path $output_dir -ChildPath raw${i}.osw
    If ($use_window -eq 1)
    {
        $global:win=$global:win_left
        $out_path=Join-Path -Path $output_dir -ChildPath raw${i}.left.osw
        run_one_window $raw $out_path
        
        $global:win=$global:win_right
        $out_path=Join-Path -Path $output_dir -ChildPath raw${i}.right.osw
        run_one_window $raw $out_path
    }
    else
    {
        run_one $raw $out_path
    }
    $i++
}

$subsample=1/(@($raw_files).Length)
run_pyprophet $output_dir $subsample

