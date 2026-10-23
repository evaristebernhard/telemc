# -*- coding: utf-8 -*-
from pathlib import Path
import sys
import numpy as np

try:
    from data_manip.extraction.telemac_file import TelemacFile
except ImportError:
    print("ERROR: Run inside the configured openTELEMAC environment.")
    sys.exit(1)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
RESULT_FILE = PROJECT_DIR / "05_run" / "local_breach_L2_B1_test_res.slf"

XMIN, XMAX = -113.0, 113.0
YMIN, YMAX = -25.0, 35.0
ZMIN = 27.0
TOL = 1.0e-7
TARGETS = [0.0, 1.0, 10.0, 30.0, 60.0]

def get_bottom(res, rec):
    names=[str(v).strip() for v in res.varnames]
    for name in ("BOTTOM","FOND"):
        if name in names:
            return np.asarray(res.get_data_value(name,rec),dtype=float)
    raise RuntimeError("BOTTOM not found: "+str(names))

def nearest(times,target):
    return int(np.argmin(np.abs(times-target)))

def main():
    print("="*78)
    print("TUANZHOU L2-B1 60 s EROSION-STABILISATION CHECK")
    print("="*78)
    if not RESULT_FILE.exists():
        raise FileNotFoundError(RESULT_FILE)

    res=TelemacFile(str(RESULT_FILE))
    try:
        x=np.asarray(res.meshx,float); y=np.asarray(res.meshy,float)
        times=np.asarray(res.times,float)
        b0=get_bottom(res,0)
        zone=(x>=XMIN)&(x<=XMAX)&(y>=YMIN)&(y<=YMAX)
        initial20=np.abs(x)<=10.0
        hist=[]

        print(f"Nodes   : {len(x)}")
        print(f"Records : {len(times)}")
        print(f"Time    : {times[0]:.3f} -> {times[-1]:.3f} s\n")

        for target in TARGETS:
            rec=nearest(times,target); b=get_bottom(res,rec)
            d=b0-b; er=d>TOL; idx=int(np.argmax(d))
            hist.append((times[rec],float(np.max(d))))
            print(f"t={times[rec]:7.3f} s | record={rec}")
            print(f"  eroded nodes                  = {np.count_nonzero(er)}")
            print(f"  max cumulative lowering [m]   = {np.max(d):.9e}")
            print(f"  sum cumulative lowering [m]   = {np.sum(d[er]):.9e}")
            print(f"  max bed raising [m]            = {np.max(b-b0):.9e}")
            print(f"  strongest x,y [m]              = {x[idx]:.6f}, {y[idx]:.6f}")
            print(f"  eroded outside activation zone = {np.count_nonzero(er & ~zone)}")
            print(f"  eroded inside initial 20 m     = {np.count_nonzero(er & initial20)}")
            print(f"  eroded outside initial 20 m    = {np.count_nonzero(er & ~initial20)}")
            print(f"  nodes reaching Z_MIN           = {np.count_nonzero(er & (b <= ZMIN+1e-6))}\n")

        bf=get_bottom(res,len(times)-1); df=b0-bf; erf=df>TOL
        checks={
            "Reached 60 s": times[-1]>=60.0-1e-6,
            "Physical erosion occurred": np.max(df)>TOL,
            "No bed raising": np.max(bf-b0)<=TOL,
            "All erosion inside activation zone": np.count_nonzero(erf & ~zone)==0,
            "No BOTTOM below Z_MIN": np.count_nonzero(bf < ZMIN-1e-6)==0,
        }
        print("="*78)
        for k,v in checks.items(): print(f"{k:38s}: {v}")
        print(f"Final max cumulative lowering [m] : {np.max(df):.9e}")
        print(f"Final cumulative eroded nodes      : {np.count_nonzero(erf)}")
        print("="*78)
        print("L2-B1 60 s TECHNICAL CHECK PASSED" if all(checks.values())
              else "L2-B1 60 s TECHNICAL CHECK FAILED")
        print("="*78)
    finally:
        res.close()

if __name__=="__main__":
    main()
