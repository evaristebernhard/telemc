# -*- coding: utf-8 -*-
"""Add five tau_c=0.1 Pa / 100 s cases to stage 17, preserving old cases.
Run in configured openTELEMAC PowerShell.
"""
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime

ROOT = Path(r'F:\Telemac_Projects\Tuanzho_New\17_B2_Tingsanchali_100s_ABC_comparison')
RUN = ROOT / '05_run'
LOGS = ROOT / '06_output' / 'kd_run_logs'
CASES = [
    ('KD050', 'B2_Tingsanchali_erosion_kd_0p5e4.cas'),
    ('KD100', 'B2_Tingsanchali_erosion_kd_1p0e4.cas'),
    ('KD250', 'B2_Tingsanchali_erosion_kd_2p5e4.cas'),
    ('KD500', 'B2_Tingsanchali_erosion_ref.cas'),
    ('KD800', 'B2_Tingsanchali_erosion_kd_8p0e4.cas'),
]

def assignment(s, key):
    pat = re.compile(r"(?im)^\s*" + re.escape(key) + r"\s*=\s*['\"]?([^'\"\r\n/!]+)")
    matches = [m.group(1).strip() for m in pat.finditer(s)]
    if len(matches) != 1:
        raise RuntimeError(f'Expected one {key} entry, found {matches}')
    return matches[0]

def prepare():
    if not RUN.is_dir():
        raise RuntimeError(f'Missing stage 17: {RUN}')
    plans = []
    for label, source_cas in CASES:
        src = RUN / source_cas
        if not src.is_file():
            raise RuntimeError(f'Missing {src}')
        cas = src.read_text(encoding='utf-8-sig')
        if abs(float(assignment(cas, 'TIME STEP').replace('D','E')) - 0.01) > 1e-10:
            raise RuntimeError(f'{source_cas}: time step not 0.01')
        duration = assignment(cas, 'DURATION')
        if abs(float(duration.replace('D','E')) - 100) > 1e-8:
            raise RuntimeError(f'{source_cas}: duration not 100s ({duration})')
        fortran_name = Path(assignment(cas, 'FORTRAN FILE')).name
        result_name = Path(assignment(cas, 'RESULTS FILE')).name
        original_f = RUN / fortran_name
        if not original_f.is_file():
            raise RuntimeError(f'Missing Fortran: {original_f}')
        f = original_f.read_text(encoding='utf-8-sig')
        # Only change the active TAU_C assignment; leave HMIN_EROSION=0.01 untouched.
        pattern = re.compile(r'(?im)^(\s*TAU_C\s*=\s*)([0-9.]+(?:[dDeE][+-]?\d+)?)(\s*(?:!.*)?)$')
        matches = list(pattern.finditer(f))
        if len(matches) != 1:
            raise RuntimeError(f'{fortran_name}: expected exactly one active TAU_C assignment; found {len(matches)}')
        old_value = float(matches[0].group(2).replace('D','E').replace('d','e'))
        if abs(old_value - 0.01) > 1e-10:
            raise RuntimeError(f'{fortran_name}: TAU_C={old_value}, expected 0.01')
        new_f = pattern.sub(lambda m: m.group(1) + '0.1D0' + m.group(3), f, count=1)
        tag = f'TC0100_{label}'
        new_cas_name = f'B2_100s_{tag}.cas'
        new_f_name = f'user_condin_h_100s_{tag}.f'
        new_res_name = f'B2_100s_{tag}_res.slf'
        updated_cas = cas.replace(fortran_name, new_f_name).replace(result_name, new_res_name)
        if updated_cas == cas or assignment(updated_cas,'FORTRAN FILE') != new_f_name or assignment(updated_cas,'RESULTS FILE') != new_res_name:
            raise RuntimeError(f'Failed CAS rewriting for {label}')
        paths = [RUN / new_cas_name, RUN / new_f_name, RUN / new_res_name]
        plans.append((label, paths, updated_cas, new_f))
    # No overwrite: accept only an identical prepared input on a retry.
    for label, (cas_path, f_path, res_path), cas_text, f_text in plans:
        for p, expected in [(cas_path, cas_text), (f_path, f_text)]:
            if p.exists() and p.read_text(encoding='utf-8-sig') != expected:
                raise RuntimeError(f'Refusing to overwrite modified file: {p}')
    for label, (cas_path, f_path, res_path), cas_text, f_text in plans:
        if not cas_path.exists(): cas_path.write_text(cas_text, encoding='utf-8')
        if not f_path.exists(): f_path.write_text(f_text, encoding='utf-8')
        print(f'Prepared {label}: tau_c=0.1 Pa, 100 s; result={res_path.name}', flush=True)
    return plans

def run(plans):
    LOGS.mkdir(parents=True, exist_ok=True)
    check = subprocess.run(['powershell.exe','-NoProfile','-Command',
                            'Get-Command telemac2d.py | Select-Object -ExpandProperty Source'],
                           cwd=str(RUN), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, errors='replace')
    if check.returncode or not check.stdout.strip():
        raise RuntimeError('Cannot resolve telemac2d.py:\n' + check.stdout)
    for label, (cas_path, f_path, res_path), _, _ in plans:
        log_path = LOGS / f'TC0100_{label}_100s_run.log'
        if res_path.exists():
            if log_path.is_file() and 'CORRECT END OF RUN' in log_path.read_text(encoding='utf-8', errors='replace'):
                print(f'{label}: existing completed result preserved; skip', flush=True)
                continue
            raise RuntimeError(f'Existing unverified result, refusing overwrite: {res_path}')
        print(f'[{datetime.now():%H:%M:%S}] Running TC0100 {label}', flush=True)
        ps = f"telemac2d.py '.\\{cas_path.name}'"
        proc = subprocess.run(['powershell.exe','-NoProfile','-Command',ps], cwd=str(RUN),
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, errors='replace')
        log_path.write_text(proc.stdout, encoding='utf-8')
        success = proc.returncode == 0 and 'CORRECT END OF RUN' in proc.stdout and res_path.is_file()
        print(f'  return={proc.returncode}, correct_end={"CORRECT END OF RUN" in proc.stdout}, slf={res_path.is_file()}', flush=True)
        if not success:
            print('\n'.join(proc.stdout.splitlines()[-35:]))
            raise RuntimeError(f'Failed {label}; see {log_path}')
    print('SUCCESS: five tau_c=0.1 Pa / 100s cases finished; original tau_c=0.01 results untouched.')

if __name__ == '__main__':
    run(prepare())
