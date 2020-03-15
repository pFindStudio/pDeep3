Set-Location -Path "e:/DIATools/openswath/OpenMS-2.4.0-nightly-2019-08-09/bin/"


# $output_dir="e:\DIAData\PECAN\20mz\pan-human-osw"
# $global:lib_TDA="e:/DIATools/openswath/library/phl004_s32_TDA.pqp"

$output_dir="e:\DIAData\PECAN\20mz\phl-pDeep-QE27-osw"
$lib="e:/DIATools/openswath/library/pDeep/phl_pDeep_QE27.tsv"
$global:lib_TDA="e:/DIATools/openswath/library/pDeep/phl_pDeep_QE27_TDA.pqp"

$raw_dir="e:\DIAData\PECAN\20mz"


$global:irt="e:/DIATools/openswath/library/phl004_s32_iRT_cells.tsv"
# $global:lib_TDA="e:/DIATools/openswath/sample_data/Mtb_TubercuList-R27_iRT_UPS_TDA.pqp"
# $global:irt="e:/DIATools/openswath/sample_data/iRTassays.TraML"
$global:win="e:/DIATools/openswath/sample_data/SWATHwindows_analysis.tsv"

New-Item -Path $output_dir -ItemType Directory
if (!(Test-Path $global:lib_TDA -PathType leaf))
{
    Write-Host "========== Generate decoy and convert to PQP =========="
    $pqp=-join($lib.substring(0, $lib.Length-3), 'pqp')
    TargetedFileConverter.exe -in $lib -out $pqp
    OpenSWATHDecoyGenerator.exe -in $pqp -out $global:lib_TDA
}
function run_one($raw, $out)
{
    OpenSwathWorkflow.exe '-in' $raw -tr $global:lib_TDA -sort_swath_maps -readOptions cache -tempDirectory E:/Temp -batchSize 10000 -out_osw $out -threads 4 -tr_irt $global:irt -min_coverage 0.01 -min_rsq 0.95
    # OpenSwathWorkflow.exe '-in' $raw -tr $global:lib_TDA -sort_swath_maps -readOptions cache -tempDirectory E:/Temp -batchSize 10000 -swath_windows_file $global:win -tr_irt $global:irt -out_osw $out -threads 4
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
        $run_reduced = -join(${run},"r") # generates .oswr files
        pyprophet reduce --in=$run --out=$run_reduced
        $reduced_files += $run_reduced
    }
        
    Write-Host "========== FDR =========="
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

