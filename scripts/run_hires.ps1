<#
Fine (0.025-degree) product and the fixed-radius / split-half fits of the resolution study.

Usage: pwsh scripts/run_hires.ps1 [-Work work] [-Python python] [-Workers 22]

Reads <Work>/upper_plate_declustered.npz (output of 03_decluster_nn.py) and writes the fit_*.npz
files read by 09_resolution_study.py, plus a log hires_runner.log, into <Work>.
About 75 minutes on 22 cores.
#>
param(
  [string]$Work = "work",
  [string]$Python = "python",
  [int]$Workers = 22
)
$ErrorActionPreference = "Continue"
$env:PYTHONUTF8 = "1"
$S = $PSScriptRoot
$W = (Resolve-Path $Work).Path
$log = "$W/hires_runner.log"
$src = "$W/upper_plate_declustered.npz"
$nw = "$Workers"
function Invoke-FitStep($name, $argv) {
  "=== $name  start $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Out-File -Append $log
  & $Python -I @argv 2>&1 | Out-File -Append $log
  "=== $name  end   $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')  exit $LASTEXITCODE" | Out-File -Append $log
}
Invoke-FitStep "fine 0.025" @("$S/04_fit_seismogenic_bottom.py", $src, "$W/fit_v3_0p025.npz", "--step","0.025","--workers",$nw,"--radii","3,5,7,10,15,20,30,40")
foreach ($r in @(2,3,5,10,20)) {
  Invoke-FitStep "full r=$r" @("$S/04_fit_seismogenic_bottom.py", $src, "$W/fit_full_r${r}_0p1.npz", "--step","0.1","--workers",$nw,"--radii","$r")
  foreach ($h in @(0,1)) {
    Invoke-FitStep "half$h r=$r" @("$S/04_fit_seismogenic_bottom.py", $src, "$W/fit_half${h}_r${r}_0p1.npz", "--step","0.1","--workers",$nw,"--radii","$r","--ess-min","15","--nmin","30","--half","$h")
  }
}
foreach ($box in @(@("kyushu","129.5,132.5,31,34"), @("tohoku_kanto","138,142,35,40"))) {
  foreach ($h in @(0,1)) {
    Invoke-FitStep "half$h fine $($box[0])" @("$S/04_fit_seismogenic_bottom.py", $src, "$W/fit_half${h}_0p025_$($box[0]).npz", "--step","0.025","--workers",$nw,"--radii","3,5,7,10,15,20,30,40","--ess-min","15","--nmin","30","--half","$h","--bbox",$box[1])
  }
}
"=== ALL DONE $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Out-File -Append $log
