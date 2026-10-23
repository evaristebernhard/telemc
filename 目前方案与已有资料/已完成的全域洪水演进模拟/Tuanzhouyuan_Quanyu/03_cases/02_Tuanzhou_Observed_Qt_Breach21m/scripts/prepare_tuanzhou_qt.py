#!/usr/bin/env python3
"""Prepare a standalone TELEMAC-2D observed-Q(t) case, preserving source inputs."""
import csv, math, re, struct, shutil
from pathlib import Path

ROOT = Path(r'F:\Telemac_Projects\Tuanzhouyuan_Quanyu')
CASE = ROOT / '03_cases' / '05_Tuanzhou_Observed_Qt_Breach21m'
SRC_GEO = ROOT / '02_input' / 'tuanzhou_geo.slf'
SRC_CLI = ROOT / '02_input' / 'tuanzhou_test_226m.cli'
SRC_CAS = ROOT / '03_cases' / '01_Q2000_3h_baseline' / 'run' / 'tuanzhou_test_q2000_3h.cas'
BREACH_IDS = [10595, 14059, 4178, 10593, 10584, 10594, 7449]
DATA = '''0.29915 365.2174
0.47009 573.913
0.59829 747.8261
0.81197 973.913
0.98291 1130.435
1.11111 1304.348
1.32479 1513.043
1.49573 1669.565
1.62393 1826.087
1.83761 2034.783
1.96581 2139.13
2.13675 2313.043
2.39316 2521.739
2.5641 2660.87
2.64957 2765.217
2.77778 2869.565
2.90598 2973.913
3.16239 3165.217
3.46154 3391.304
3.67521 3547.826
3.97436 3739.13
4.18803 3878.261
4.44444 4017.391
4.74359 4173.913
4.95726 4313.043
5.25641 4434.783
5.47009 4521.739
5.76923 4660.87
6.28205 4834.783
6.49573 4886.957
6.88034 4973.913
7.30769 5026.087
7.73504 5078.261
8.16239 5113.043
8.58974 5130.435
9.01709 5095.652
9.44444 5043.478
9.87179 5008.696
10.29915 4939.13
10.64103 4921.739
11.02564 4834.783
11.41026 4782.609
11.83761 4730.435
12.26496 4695.652
12.69231 4660.87
13.07692 4660.87
13.2906 4556.522
13.58974 4452.174
13.97436 4365.217
14.23077 4278.261
14.52991 4191.304
14.82906 4104.348
15.08547 4000
15.29915 3895.652
15.47009 3739.13
15.68376 3565.217
15.89744 3408.696
16.06838 3269.565
16.28205 3113.043
16.45299 3008.696
16.5812 2869.565
16.75214 2713.043
16.88034 2556.522
17.05128 2382.609
17.22222 2208.696
17.35043 2052.174
17.39316 1878.261
17.64957 1773.913
17.86325 1634.783
18.07692 1495.652
18.16239 1339.13
18.37607 1234.783
18.63248 1060.87
18.88889 921.7391
19.10256 800
19.2735 695.6522
19.52991 608.6957
19.82906 486.9565
20.12821 382.6087'''

def records(f):
    start = f.tell(); head = f.read(4)
    if not head: return
    if len(head) != 4: raise ValueError('Truncated SLF record marker')
    n = struct.unpack('>i', head)[0]
    if n < 0 or n > 500000000: raise ValueError('Unexpected SLF record length or byte order')
    payload = f.tell(); f.seek(n, 1)
    tail = f.read(4)
    if len(tail) != 4 or struct.unpack('>i', tail)[0] != n:
        raise ValueError('SLF record marker mismatch')
    return payload, n

def patch_slf(source, target):
    shutil.copy2(source, target)
    with target.open('r+b') as f:
        def rd():
            rec = records(f)
            if rec is None: raise ValueError('Unexpected end of SLF')
            pos, n = rec; f.seek(pos); data = f.read(n); f.seek(pos+n+4)
            return pos,n,data
        _,_,title = rd()
        if b'SERAFIN' not in title.upper() and b'SELAFIN' not in title.upper():
            print('SLF title:', repr(title[:80]))
        _,n,b = rd(); assert n == 8, 'Unexpected variable count record'
        nbv1,nbv2 = struct.unpack('>2i',b); nvar=nbv1+nbv2
        names=[]
        for _ in range(nvar):
            _,n,b=rd(); names.append(b[:16].decode('latin1').strip().upper())
        _,n,b=rd(); iparam=struct.unpack('>'+str(n//4)+'i',b)
        if len(iparam)>9 and iparam[9] == 1: rd()
        _,n,b=rd(); nelem,npoin,ndp,*rest=struct.unpack('>'+str(n//4)+'i',b)
        rd(); rd(); rd(); rd() # IKLE, IPOBO, X, Y
        print('SLF:', 'NPOIN=',npoin,'NELEM=',nelem,'variables=',names)
        bottom_indices=[i for i,name in enumerate(names) if name in ('BOTTOM','FOND','BOTTOM ELEVATION','BOTTOM LEVEL') or name.startswith('BOTTOM')]
        if len(bottom_indices)!=1: raise ValueError('Expected exactly one BOTTOM variable; found '+str(bottom_indices))
        _,tn,tdata=rd() # first time
        if tn not in (4,8): raise ValueError('Unexpected time record')
        selected=bottom_indices[0]
        old=[]
        for vi in range(nvar):
            pos,n,b=rd()
            if vi != selected: continue
            if n == 4*npoin: code='f'; size=4
            elif n == 8*npoin: code='d'; size=8
            else: raise ValueError('Unexpected BOTTOM array length')
            for node_id in BREACH_IDS:
                if not 1<=node_id<=npoin: raise ValueError('Invalid node ID '+str(node_id))
                f.seek(pos+(node_id-1)*size)
                previous=struct.unpack('>'+code,f.read(size))[0]
                old.append((node_id,previous))
                f.seek(pos+(node_id-1)*size)
                f.write(struct.pack('>'+code,21.0))
            f.seek(pos+n+4)
        print('Modified BOTTOM (node ID, old value):',old)
        return old

def main():
    for p in (SRC_GEO,SRC_CLI,SRC_CAS):
        if not p.is_file(): raise FileNotFoundError(p)
    for name in ('input','run','results','scripts'):
        (CASE/name).mkdir(parents=True,exist_ok=True)
    geo=CASE/'input'/'tuanzhou_geo_breach21m.slf'
    if geo.exists(): raise FileExistsError('Refusing to overwrite existing modified geometry: '+str(geo))
    patch_slf(SRC_GEO,geo)
    shutil.copy2(SRC_CLI, CASE/'input'/SRC_CLI.name)
    points=[(0.0,0.0)]+[tuple(map(float,line.split())) for line in DATA.splitlines()]
    assert all(points[i][0]<points[i+1][0] for i in range(len(points)-1))
    points.append((20.5,points[-1][1]))
    liq=CASE/'input'/'breach_Qt.liq'
    with liq.open('w',encoding='ascii',newline='\n') as f:
        f.write('T Q(1)\ns m3/s\n')
        for hour,flow in points: f.write(f'{hour*3600:.6f} {flow:.6f}\n')
    with (CASE/'input'/'breach_Qt_source.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(['time_h','flow_m3_s']);w.writerows(points)
    cas=SRC_CAS.read_text(encoding='utf-8')
    replacements={
        'GEOMETRY FILE':"'../input/tuanzhou_geo_breach21m.slf'",
        'BOUNDARY CONDITIONS FILE':"'../input/tuanzhou_test_226m.cli'",
        'RESULTS FILE':"'../results/tuanzhou_observed_qt_res.slf'",
        'NUMBER OF TIME STEPS':'73800',
        'GRAPHIC PRINTOUT PERIOD':'1800',
        'LISTING PRINTOUT PERIOD':'1800',
    }
    for key,value in replacements.items():
        pattern=r'(?im)^\s*'+re.escape(key)+r'\s*=\s*[^\r\n]*'
        cas,count=re.subn(pattern,key+' = '+value,cas)
        if count!=1: raise ValueError(f'Expected exactly one {key}, found {count}')
    cas,count=re.subn(r'(?im)^\s*PRESCRIBED FLOWRATES\s*=\s*[^\r\n]*',
                       "PRESCRIBED FLOWRATES = 0.0\nLIQUID BOUNDARIES FILE = '../input/breach_Qt.liq'",cas)
    if count!=1: raise ValueError('PRESCRIBED FLOWRATES not found')
    cas=cas.replace('TUANZHOUYUAN TEST FIXED BREACH Q2000 3H','TUANZHOUYUAN OBSERVED QT BREACH 21M')
    (CASE/'run'/'tuanzhou_observed_qt.cas').write_text(cas,encoding='ascii')
    print('Prepared:',CASE)
    print('Boundary points:',len(points),'last time:',points[-1][0]*3600,'s')
    print('Run in openTELEMAC 9.1 PowerShell:')
    print('cd "'+str(CASE/'run')+'"')
    print('telemac2d.py tuanzhou_observed_qt.cas')
    print('IMPORTANT: confirm geometry vertical datum and Q(t) represent compatible elevations/flows.')

if __name__=='__main__': main()
