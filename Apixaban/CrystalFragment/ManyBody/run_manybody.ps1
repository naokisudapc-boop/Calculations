# run_manybody.ps1 -- run all XTB2 single points under this folder (small jobs first).
# Jobs whose .out already ends with 'ORCA TERMINATED NORMALLY' are skipped. Run from this folder.
$env:XTBEXE = "C:\ORCA_6.1.1\xTB_bigstack\xtb.exe"
Get-ChildItem -Recurse -Filter '*_XTB2_SP.inp' | Sort-Object Length | ForEach-Object {
    $name = $_.Name
    $base = [IO.Path]::GetFileNameWithoutExtension($name)
    $out  = Join-Path $_.DirectoryName ($base + '.out')
    if ((Test-Path $out) -and (Select-String -Path $out -Pattern 'ORCA TERMINATED NORMALLY' -Quiet)) { "skip $name"; return }
    "run  $name"
    Push-Location $_.DirectoryName
    cmd /c "`"C:\ORCA_6.1.1\orca.exe`" `"$name`" > `"$base.out`""
    Pop-Location
}
