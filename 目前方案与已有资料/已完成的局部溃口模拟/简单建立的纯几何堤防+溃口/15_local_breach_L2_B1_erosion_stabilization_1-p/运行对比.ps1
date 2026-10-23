$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = Split-Path -Parent $here
$old = Join-Path $root '15_local_breach_L2_B1_erosion_stabilization\05_run\local_breach_L2_B1_600s_res.slf'
$new = Join-Path $here '05_run\local_breach_L2_B1_1-p_600s_res.slf'
$out = Join-Path $here '06_comparison_1-p'
if (!(Test-Path $old)) { throw "原15结果文件不存在: $old" }
if (!(Test-Path $new)) { throw "新15结果文件不存在: $new" }
$python = Get-Command python -ErrorAction Stop
& $python.Source (Join-Path $here 'compare_15_vs_1p.py') --original $old --modified $new --output $out
if ($LASTEXITCODE -ne 0) { throw "绘图程序失败，退出代码: $LASTEXITCODE" }
Write-Host "对比完成。图片和表格位于: $out"
