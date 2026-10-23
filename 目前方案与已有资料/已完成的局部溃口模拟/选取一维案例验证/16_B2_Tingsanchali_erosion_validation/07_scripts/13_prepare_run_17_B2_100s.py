# -*- coding: utf-8 -*-
"""Create isolated stage 17, run five 100-s B-2 TELEMAC cases, compare A/B/C.
Run from configured openTELEMAC PowerShell:
 python .\07_scripts\13_prepare_run_17_B2_100s.py
"""
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime

BASE = Path(r'F:\Telemac_Projects\Tuanzho_New')
SRC = BASE / '16_B2_Tingsanchali_erosion_validation'
DST = BASE / '17_B2_Tingsanchali_100s_ABC_comparison'
CASES = [
 ('KD050', 'B2_Tingsanchali_erosion_kd_0p5e4.cas', 'B2_Tingsanchali_erosion_kd_0p5e4_res.slf'),
 ('KD100', 'B2_Tingsanchali_erosion_kd_1p0e4.cas', 'B2_Tingsanchali_erosion_kd_1p0e4_res.slf'),
 ('KD250', 'B2_Tingsanchali_erosion_kd_2p5e4.cas', 'B2_Tingsanchali_erosion_kd_2p5e4_res.slf'),
 ('KD500', 'B2_Tingsanchali_erosion_ref.cas', 'B2_Tingsanchali_erosion_ref_res.slf'),
 ('KD800', 'B2_Tingsanchali_erosion_kd_8p0e4.cas', 'B2_Tingsanchali_erosion_kd_8p0e4_res.slf'),
]

def fail(message):
    raise RuntimeError(message)

def replace_duration(content, name):
    # TELEMAC uses NUMBER OF TIME STEPS or DURATION, sometimes both.
    # Change every active occurrence and preserve all other steering parameters.
    lines = content.splitlines(keepends=True)
    nstep = 0
    nduration = 0
    dt_values = []
    out = []
    for line in lines:
        m = re.match(r'^(\s*)(NUMBER OF TIME STEPS|DURATION|TIME STEP)\s*=\s*([^/!\r\n]*)(.*)$', line, re.I)
        if not m or line.lstrip().startswith(('/', '!')):
            out.append(line)
            continue
        key = m.group(2).upper()
        raw = m.group(3).strip()
        if key == 'TIME STEP':
            dt_values.append(float(raw.replace('D','E').replace('d','e')))
            out.append(line)
        elif key == 'NUMBER OF TIME STEPS':
            nstep += 1
            out.append(f'{m.group(1)}{m.group(2)} = 10000{m.group(4)}')
        else:
            nduration += 1
            out.append(f'{m.group(1)}{m.group(2)} = 100.{m.group(4)}')
    if len(dt_values) != 1 or abs(dt_values[0] - 0.01) > 1e-10:
        fail(f'{name}: expected exactly one TIME STEP = 0.01, got {dt_values}')
    if nstep + nduration == 0:
        fail(f'{name}: no active NUMBER OF TIME STEPS or DURATION found')
    if nstep > 1 or nduration > 1:
        fail(f'{name}: ambiguous duration keywords')
    return ''.join(out), nstep, nduration

def prepare():
    if not SRC.is_dir():
        fail(f'Missing source folder: {SRC}')
    if DST.exists():
        print(f'Stage 17 already exists; will not overwrite: {DST}')
        return
    for _, cas, _ in CASES:
        if not (SRC / '05_run' / cas).is_file():
            fail(f'Missing source CAS: {cas}')
    if not (SRC / '07_scripts' / '12_compare_B2_ABC_measured.py').is_file():
        fail('Missing 12_compare_B2_ABC_measured.py in source 07_scripts')
    print('Creating stage 17 (no old SLF/results/logs copied) ...', flush=True)
    DST.mkdir()
    for sub in ['01_geometry', '02_mesh', '03_boundary', '04_input']:
        origin = SRC / sub
        if origin.is_dir():
            shutil.copytree(origin, DST / sub)
    for sub in ['05_run', '06_output', '07_scripts', '08_figures']:
        (DST / sub).mkdir(exist_ok=True)
    # Copy only input assets, not previous results, temporary files or old logs.
    for item in (SRC / '05_run').iterdir():
        if item.is_file() and item.suffix.lower() in {'.cas', '.f', '.f90', '.cli'}:
            shutil.copy2(item, DST / '05_run' / item.name)
        elif item.is_file() and item.suffix.lower() == '.slf' and 'res' not in item.stem.lower():
            shutil.copy2(item, DST / '05_run' / item.name)
    # Preserve existing scripts for traceability; do not copy generated plots.
    for item in (SRC / '07_scripts').glob('*.py'):
        shutil.copy2(item, DST / '07_scripts' / item.name)
    for label, cas, _ in CASES:
        p = DST / '05_run' / cas
        updated, ns, nd = replace_duration(p.read_text(encoding='utf-8-sig'), cas)
        p.write_text(updated, encoding='utf-8')
        print(f'  {label}: 100 s; steps={ns}, duration={nd}')
    # Make a stage-17-specific copy of the existing comparison script.
    p = DST / '07_scripts' / '12_compare_B2_ABC_measured.py'
    source = p.read_text(encoding='utf-8-sig')
    old = "16_B2_Tingsanchali_erosion_validation"
    if source.count(old) != 1:
        fail('Unexpected comparison-script ROOT; inspect manually')
    source = source.replace(old, '17_B2_Tingsanchali_100s_ABC_comparison')
    if source.count('xlim=(0,30)') != 2:
        fail('Unexpected comparison-script plot limits; inspect manually')
    source = source.replace('xlim=(0,30)', 'xlim=(0,100)')
    p.write_text(source, encoding='utf-8')
    print('Stage 17 prepared:', DST)

def run_cases():
    run = DST / '05_run'
    logs = DST / '06_output' / 'kd_run_logs'
    logs.mkdir(parents=True, exist_ok=True)
    check = subprocess.run(['powershell.exe','-NoProfile','-Command',
                            'Get-Command telemac2d.py | Select-Object -ExpandProperty Source'],
                           cwd=str(run), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, errors='replace')
    if check.returncode != 0 or not check.stdout.strip():
        fail('telemac2d.py unavailable in inherited PowerShell environment:\n' + check.stdout)
    for label, cas, result in CASES:
        path = run / result
        log = logs / f'{label}_100s_run.log'
        if path.exists():
            if log.is_file() and 'CORRECT END OF RUN' in log.read_text(encoding='utf-8',errors='replace'):
                print(f'{label}: previous completed result retained; skipping', flush=True)
                continue
            fail(f'{label}: result exists without verified log, refusing overwrite: {path}')
        print(f'[{datetime.now():%H:%M:%S}] Running {label}: {cas}', flush=True)
        ps = f"telemac2d.py '.\\{cas}'"
        proc = subprocess.run(['powershell.exe','-NoProfile','-Command',ps], cwd=str(run),
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, errors='replace')
        log.write_text(proc.stdout, encoding='utf-8')
        good = proc.returncode == 0 and 'CORRECT END OF RUN' in proc.stdout and path.is_file()
        print(f'  return={proc.returncode}, correct_end={"CORRECT END OF RUN" in proc.stdout}, slf={path.is_file()}', flush=True)
        if not good:
            print('\n'.join(proc.stdout.splitlines()[-35:]))
            fail(f'{label} FAILED. See {log}')
    print('Five 100-s runs completed. Running A/B/C comparison...', flush=True)
    script = DST / '07_scripts' / '12_compare_B2_ABC_measured.py'
    proc = subprocess.run([sys.executable, str(script)], cwd=str(DST))
    if proc.returncode:
        fail('Simulations completed but A/B/C comparison failed; rerun comparison separately')
    print('DONE: 08_figures (4 PNG); 06_output (2 CSV)')

if __name__ == '__main__':
    prepare()
    run_cases()
