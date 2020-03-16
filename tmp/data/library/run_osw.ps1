Set-Location -Path "e:/DIATools/openswath/OpenMS-2.4.0-nightly-2019-08-09/bin/"


# $output_dir="e:\DIAData\PECAN\5mz\test-ms1-traces-phl"
# $global:lib_TDA="e:/DIATools/openswath/library/phl004_s32_TDA.pqp"
# $global:irt="e:/DIATools/openswath/library/phl004_s32_iRT_cells.tsv"

# $output_dir="e:\DIAData\PECAN\5mz\phl-pDeep-QE27-mixed-osw"
# $lib="e:/DIATools/openswath/library/pDeep/phl_pDeep_QE27.tsv"
# $global:lib_TDA="e:/DIATools/openswath/library/pDeep/phl_pDeep_QE27_mixed_TDA.pqp"
# $global:irt="e:/DIATools/openswath/library/pDeep/pDeep_iRT_cells.tsv"

$output_dir="e:\DIAData\PECAN\plasma\peptide-list-mixed-osw"
$lib="e:/DIATools/openswath/library/pDeep/peptide_list_QE27.tsv"
$global:lib_TDA="e:/DIATools/openswath/library/pDeep/peptide_list_QE27_TDA.pqp"
$global:irt="e:/DIATools/openswath/library/pDeep/blood_RT_proteins_fasta_QE27.tsv"

$raw_dir="e:\DIAData\PECAN\plasma"

$ms1_tol_type='ppm' # or 'Th'
$ms1_tol=10
$ms2_tol_type='ppm' # or 'Th'
$ms2_tol=20
$batch=5000
$thread=4
$overlapDIA='false'

# $global:win="e:/DIATools/openswath/sample_data/SWATHwindows_analysis.tsv"

function Compare-File-Date($f1, $f2)
{
    return ([datetime](Get-ItemProperty -Path $f1 -Name LastWriteTime).lastwritetime).ToOADate() -lt ([datetime](Get-ItemProperty -Path $f2 -Name LastWriteTime).lastwritetime).ToOADate()
}

New-Item -Path $output_dir -ItemType Directory
if (!(Test-Path $global:lib_TDA -PathType leaf))
{
    Write-Host "========== Generating decoy =========="
    $pqp=-join($lib.substring(0, $lib.Length-3), 'pqp')
    TargetedFileConverter.exe -in $lib -out $pqp
    Write-Host "========== Converting to PQP =========="
    OpenSWATHDecoyGenerator.exe -in $pqp -out $global:lib_TDA -method reverse
}
function run_one($raw, $out)
{
    OpenSwathWorkflow.exe -in $raw -tr $global:lib_TDA -sort_swath_maps -readOptions cache -tempDirectory E:/Temp -batchSize $batch -out_osw $out -threads $thread -tr_irt $global:irt -rt_extraction_window 600 -mz_extraction_window_unit $ms2_tol_type -mz_extraction_window $ms2_tol -mz_extraction_window_ms1_unit $ms1_tol_type -mz_extraction_window_ms1 $ms1_tol -matching_window_only -matching_window_only -$overlapDIA -min_coverage 0.01 -min_rsq 0.95
    # OpenSwathWorkflow.exe '-in' $raw -tr $global:lib_TDA -sort_swath_maps -readOptions cache -tempDirectory E:/Temp -batchSize 10000 -swath_windows_file $global:win -tr_irt $global:irt -out_osw $out -threads 4 -use_ms1_traces
    
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
    run_one $raw $out_path
    $i++
}

$subsample=1/(@($raw_files).Length)
run_pyprophet $output_dir $subsample

