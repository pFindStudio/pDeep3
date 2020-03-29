param (
    [string]$mzml = $(throw "-mzml is required."),
    [string]$fasta = $(throw "-fasta is required."),
    [switch]$tune = $false,
    [string]$output = "",
    [string]$lib = "e:\DIATools\encyclopedia\library\phl_pdeep_QE27.dlib",
    [switch]$overlap = $false,
    [switch]$phosite = $false,
    [string]$enc_path = "e:\DIATools\encyclopedia\"
)

$thread = 6
$enc = 'encyclopedia-0.9.0-executable.jar'
$thr = 'thesaurus-0.9.4-executable.jar'
$mod = 'Phosphorylation'

$dlib = $lib
$elib = $output

Set-Location -Path $enc_path


# function run_enc_overlap()
# {
    # $res = -join($mzml,'.encyclopedia.txt')
    # Remove-Item $res
    # java -Xmx14g -jar $enc -i $mzml -l $dlib -f $fasta -numberOfThreadsUsed 6 -acquisition overlapping
    # if ($elib.Length -ne 0)
    # {
        # $saved_elib = -join($mzml,'.elib')
        # Rename-Item -Path $saved_elib -NewName $elib -Force
    # }
# }

function run_enc()
{
    # $res = -join($mzml,'.encyclopedia.txt')
    Write-Host delete $(-join($mzml,'.encyclopedia.txt'))
    Remove-Item -Path $(-join($mzml,'.encyclopedia.txt'))
    if ($overlap)
    {
        java -Xmx14g -jar $enc -i $mzml -l $dlib -f $fasta -numberOfThreadsUsed $thread -acquisition overlapping
    }
    else
    {
        java -Xmx14g -jar $enc -i $mzml -l $dlib -f $fasta -numberOfThreadsUsed $thread
    }
    if ($elib.Length -ne 0)
    {
        $saved_elib = -join($mzml,'.elib')
        Move-Item -Path $saved_elib -Destination $elib -Force
    }
    # $newelib = -join('raw=', $raw_id, '.', (Split-Path $dlib -leaf).split('.')[0], '.elib')
}

function run_thr()
{
    # $res_files = Get-ChildItem -Path . -Recurse -Include *thesaurus*
    # $res = -join($mzml,'.thesaurus.txt')
    # Write-Host delete $(-join($mzml,'.encyclopedia.txt'))
    # Remove-Item -Path $(-join($mzml,'.encyclopedia.txt'))
    Remove-Item $(-join($mzml, '*.thesaurus.txt.*.txt'))
    if ($overlap)
    {
        java -Xmx14g -jar $thr -i $mzml -l $dlib -f $fasta -localizationModification Phosphorylation -acquisition overlapping -numberOfThreadsUsed $thread
    }
    else
    {
        java -Xmx14g -jar $thr -i $mzml -l $dlib -f $fasta -localizationModification Phosphorylation -numberOfThreadsUsed $thread
    }
    if ($elib.Length -ne 0)
    {
        $saved_elib = -join($mzml,'.thesaurus.txt.Phosphorylation.thesaurus.elib')
        Move-Item -Path $saved_elib -Destination $elib -Force
    }
}

Write-Host "============== Start =============="
if ($phosite)
{
    run_thr
}
else
{
    run_enc
}
