# -*- coding: utf-8 -*-
"""Compare actual TELEMAC-2D SELAFIN outputs (original B1 vs B1 with 1-p).
Requires: numpy, matplotlib. Does not fabricate, interpolate in time, or modify SLF files.
"""
import argparse
import csv
import datetime
import json
import re
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import TwoSlopeNorm

class SelafinError(Exception):
    pass

class Selafin:
    def __init__(self, filename):
        self.path = Path(filename)
        self.f = self.path.open('rb')
        self.endian = self._detect_endian()
        title = self.record()
        self.title = title[:72].decode('latin1', errors='replace').strip()
        x = self.record()
        if len(x) != 8: raise SelafinError('Invalid variable count record')
        self.nvar1, self.nvar2 = struct.unpack(self.endian+'2i', x)
        self.nvar = self.nvar1 + self.nvar2
        if not 1 <= self.nvar <= 100: raise SelafinError('Unexpected variable count')
        self.names = [self.record()[:16].decode('latin1', errors='replace').strip() for _ in range(self.nvar)]
        self.units = []  # name records include units, not separately stored
        ip = self.record()
        self.iparam = struct.unpack(self.endian+'10i', ip)
        if self.iparam[9] == 1: self.record()  # date/time record
        mesh = self.record()
        self.nelem, self.npoin, self.ndp, self.nplan = struct.unpack(self.endian+'4i', mesh)
        if self.ndp != 3: raise SelafinError('Only triangular 2D SELAFIN mesh supported')
        conn = np.frombuffer(self.record(), dtype=self.endian+'i4')
        if conn.size != self.nelem * self.ndp: raise SelafinError('Connectivity size mismatch')
        self.triangles = conn.reshape(self.nelem, self.ndp).astype(np.int64)-1
        self.record()  # IPOBO
        xraw = self.record(); yraw = self.record()
        if len(xraw) == self.npoin*4: self.dtype = 'f4'
        elif len(xraw) == self.npoin*8: self.dtype = 'f8'
        else: raise SelafinError('Invalid coordinate size')
        self.x = np.frombuffer(xraw, dtype=self.endian+self.dtype).astype(float)
        self.y = np.frombuffer(yraw, dtype=self.endian+self.dtype).astype(float)
        self.frames = []
        while True:
            pos = self.f.tell()
            try: tdata = self.record(allow_eof=True)
            except SelafinError as e: raise SelafinError(f'Invalid time record at {pos}: {e}')
            if tdata is None: break
            if len(tdata) not in (4,8): raise SelafinError(f'Invalid time at {pos}')
            t = float(np.frombuffer(tdata, dtype=self.endian+('f4' if len(tdata)==4 else 'f8'))[0])
            var_positions=[]
            for _ in range(self.nvar):
                var_positions.append(self.f.tell())
                self.skip_record()
            self.frames.append((t,var_positions))
        if not self.frames: raise SelafinError('No time frames in file')
    def _detect_endian(self):
        b=self.f.read(4); self.f.seek(0)
        if len(b)!=4: raise SelafinError('Empty file')
        for e in ('>','<'):
            n=struct.unpack(e+'i',b)[0]
            if 40<=n<=1024: return e
        raise SelafinError('Not a recognized SELAFIN Fortran-unformatted file')
    def record(self, allow_eof=False):
        h=self.f.read(4)
        if not h and allow_eof:return None
        if len(h)!=4:raise SelafinError('Incomplete record header')
        n=struct.unpack(self.endian+'i',h)[0]
        if n<0 or n>2_000_000_000:raise SelafinError('Invalid record size')
        data=self.f.read(n); tail=self.f.read(4)
        if len(data)!=n or len(tail)!=4 or struct.unpack(self.endian+'i',tail)[0]!=n:
            raise SelafinError('Corrupt/incomplete Fortran record')
        return data
    def skip_record(self):
        h=self.f.read(4)
        if len(h)!=4:raise SelafinError('Incomplete variable record')
        n=struct.unpack(self.endian+'i',h)[0]
        if n<0 or n>self.npoin*16:raise SelafinError('Invalid variable record size')
        self.f.seek(n,1)
        t=self.f.read(4)
        if len(t)!=4 or struct.unpack(self.endian+'i',t)[0]!=n:raise SelafinError('Invalid variable record tail')
    def get(self, frame_idx, var_idx):
        self.f.seek(self.frames[frame_idx][1][var_idx])
        b=self.record()
        dtype = 'f4' if len(b)==self.npoin*4 else ('f8' if len(b)==self.npoin*8 else None)
        if dtype is None:raise SelafinError('Variable array size mismatch')
        return np.frombuffer(b,dtype=self.endian+dtype).astype(float)
    def close(self): self.f.close()

def get_bed_idx(names):
    for i,n in enumerate(names):
        v=re.sub(r'[^A-Z0-9]','',n.upper())
        if v in ('BOTTOM','FONDBED','BED','BEDLEVEL','BOTTOMELEVATION','BOTTOMLEVEL','ZF','COTEDUFOND'):
            return i
    return None

def get_depth_idx(names):
    for i,n in enumerate(names):
        v=re.sub(r'[^A-Z0-9]','',n.upper())
        if v in ('WATERDEPTH','HAUTEURDEAU','DEPTH','H'):return i
    return None

def areas(x,y,tri):
    a=tri[:,0]; b=tri[:,1]; c=tri[:,2]
    ar=0.5*np.abs((x[b]-x[a])*(y[c]-y[a])-(x[c]-x[a])*(y[b]-y[a]))
    nodal=np.zeros(len(x)); np.add.at(nodal,tri.ravel(),np.repeat(ar/3.,3))
    return nodal, ar

def mesh_equal(a,b):
    return a.npoin==b.npoin and a.nelem==b.nelem and np.array_equal(a.triangles,b.triangles) and np.allclose(a.x,b.x,atol=1e-5,rtol=0) and np.allclose(a.y,b.y,atol=1e-5,rtol=0)

def save_map(ax, tri, values, title, vmin=None, vmax=None, cmap='viridis', norm=None):
    pc=ax.tripcolor(tri,values,shading='gouraud',cmap=cmap,vmin=vmin if norm is None else None,vmax=vmax if norm is None else None,norm=norm,rasterized=True)
    ax.set_title(title); ax.set_xlabel('X (m)'); ax.set_ylabel('Y (m)');ax.set_aspect('equal',adjustable='box')
    plt.colorbar(pc,ax=ax,shrink=.85)

def save_fig(path):
    plt.savefig(path,dpi=220,bbox_inches='tight');plt.close()

def main():
    ap=argparse.ArgumentParser(description='Compare real TELEMAC original B1 and 1-p SELAFIN results')
    ap.add_argument('--original',required=True,type=Path)
    ap.add_argument('--modified',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    if not args.original.is_file() or not args.modified.is_file():
        raise FileNotFoundError(f'Input missing: original={args.original.exists()}, modified={args.modified.exists()}')
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    old=Selafin(args.original); new=Selafin(args.modified)
    try:
        print('Original variables:',old.names,'frames:',len(old.frames))
        print('Modified variables:',new.names,'frames:',len(new.frames))
        if not mesh_equal(old,new):raise SelafinError('Meshes are not identical. Refusing node-by-node difference/volume comparison.')
        bo=get_bed_idx(old.names); bn=get_bed_idx(new.names)
        if bo is None or bn is None:raise SelafinError(f'BOTTOM variable not found. Original: {old.names}; Modified: {new.names}')
        t_old=np.array([f[0] for f in old.frames]);t_new=np.array([f[0] for f in new.frames]);
        if abs(t_old[-1]-t_new[-1])>1e-3:raise SelafinError(f'Final times differ: {t_old[-1]} vs {t_new[-1]}')
        z0o=old.get(0,bo); z0n=new.get(0,bn)
        if not np.allclose(z0o,z0n,atol=1e-4,rtol=0):
            raise SelafinError('Initial bed fields differ. Cannot isolate 1-p effect; check inputs and first output frame.')
        zo=old.get(-1,bo);zn=new.get(-1,bn)
        if not np.isfinite(zo).all() or not np.isfinite(zn).all():raise SelafinError('Nonfinite bed elevation found')
        z0=z0o
        erosion_old=z0-zo; erosion_new=z0-zn
        delta=erosion_new-erosion_old
        node_area, _=areas(old.x,old.y,old.triangles)
        if np.any(node_area<=0):raise SelafinError('Zero area mesh node found')
        tri=mtri.Triangulation(old.x,old.y,old.triangles)
        tmax=float(t_old[-1]); label_time=f't = {tmax:g} s'
        plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.unicode_minus':False})
        bed_lo=float(min(z0.min(),zo.min(),zn.min()));bed_hi=float(max(z0.max(),zo.max(),zn.max()))
        er_hi=float(max(erosion_old.max(),erosion_new.max(),0.001))
        # 01 side-by-side final bed
        fig,axes=plt.subplots(1,2,figsize=(14,5),constrained_layout=True)
        save_map(axes[0],tri,zo,'Original B1: final bed ('+label_time+')',bed_lo,bed_hi,'terrain')
        save_map(axes[1],tri,zn,'B1 with 1-p: final bed ('+label_time+')',bed_lo,bed_hi,'terrain')
        save_fig(out/'01_final_bed_comparison.png')
        # 02 erosion depth, same scale
        fig,axes=plt.subplots(1,2,figsize=(14,5),constrained_layout=True)
        save_map(axes[0],tri,erosion_old,'Original B1: bed lowering (first to last)',0,er_hi,'inferno')
        save_map(axes[1],tri,erosion_new,'B1 with 1-p: bed lowering (first to last)',0,er_hi,'inferno')
        save_fig(out/'02_erosion_depth_comparison.png')
        # 03 difference
        lim=float(max(np.abs(delta).max(),1e-8))
        fig,ax=plt.subplots(figsize=(10,6),constrained_layout=True)
        save_map(ax,tri,delta,'Additional bed lowering: (1-p) minus original (m)',cmap='RdBu_r',norm=TwoSlopeNorm(vcenter=0,vmin=-lim,vmax=lim))
        save_fig(out/'03_erosion_difference_1p_minus_original.png')
        # 04 final profile across point with greatest cumulative erosion, X direction
        i=int(np.argmax(np.maximum(erosion_old,erosion_new)))
        ycut=float(old.y[i]); xx=np.linspace(old.x.min(),old.x.max(),800); yy=np.full_like(xx,ycut)
        fig,ax=plt.subplots(figsize=(12,5),constrained_layout=True)
        for z,label,style in [(z0,'Initial bed','k--'),(zo,'Original B1','b-'),(zn,'B1 with 1-p','r-')]:
            zline=mtri.LinearTriInterpolator(tri,z)(xx,yy)
            ax.plot(xx,zline,style,lw=1.5,label=label)
        ax.set_title(f'Bed profile at Y = {ycut:.2f} m (passes through strongest-scour node)')
        ax.set_xlabel('X (m)');ax.set_ylabel('Bed elevation (m)');ax.grid(alpha=.25);ax.legend()
        save_fig(out/'04_final_bed_profile.png')
        # 05 temporal comparisons: only common saved timestamps, never invent time points
        matches=[]
        for j,t in enumerate(t_new):
            k=int(np.argmin(np.abs(t_old-t)))
            if abs(t_old[k]-t)<1e-3:matches.append((k,j,float(t)))
        if not matches:raise SelafinError('No common saved time points')
        metrics=[]
        for k,j,t in matches:
            do=z0-old.get(k,bo);dn=z0-new.get(j,bn)
            # Positive scour volume = integral max(z0-z,0) dA, not exported sediment volume.
            vo=float(np.dot(np.maximum(do,0),node_area));vn=float(np.dot(np.maximum(dn,0),node_area))
            metrics.append({'time_s':t,'max_lowering_original_m':float(do.max()),'max_lowering_1p_m':float(dn.max()),
                'scour_volume_original_m3':vo,'scour_volume_1p_m3':vn,
                'net_bed_change_original_m3':float(np.dot(do,node_area)),
                'net_bed_change_1p_m3':float(np.dot(dn,node_area))})
        with (out/'05_time_series_metrics.csv').open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
        fig,axes=plt.subplots(2,1,figsize=(10,9),constrained_layout=True)
        ts=[r['time_s'] for r in metrics]
        axes[0].plot(ts,[r['max_lowering_original_m'] for r in metrics],'b-o',ms=3,label='Original B1')
        axes[0].plot(ts,[r['max_lowering_1p_m'] for r in metrics],'r-o',ms=3,label='B1 with 1-p')
        axes[0].set_ylabel('Maximum bed lowering (m)');axes[0].legend();axes[0].grid(alpha=.25)
        axes[1].plot(ts,[r['scour_volume_original_m3'] for r in metrics],'b-o',ms=3,label='Original B1')
        axes[1].plot(ts,[r['scour_volume_1p_m3'] for r in metrics],'r-o',ms=3,label='B1 with 1-p')
        axes[1].set_xlabel('Saved output time (s)');axes[1].set_ylabel('Positive scour volume (m³)');axes[1].legend();axes[1].grid(alpha=.25)
        save_fig(out/'05_erosion_time_series.png')
        # Summary and per-node csv
        with (out/'06_nodewise_erosion.csv').open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.writer(f);w.writerow(['node','x_m','y_m','control_area_m2','initial_bed_m','original_final_bed_m','1p_final_bed_m','original_lowering_m','1p_lowering_m','additional_lowering_m'])
            for i in range(old.npoin):w.writerow([i+1,old.x[i],old.y[i],node_area[i],z0[i],zo[i],zn[i],erosion_old[i],erosion_new[i],delta[i]])
        report={
          'original_file':str(args.original),'modified_file':str(args.modified),
          'original_variables':old.names,'modified_variables':new.names,
          'n_nodes':old.npoin,'n_triangles':old.nelem,'last_time_s':tmax,
          'common_saved_times':len(matches),
          'original_max_lowering_m':float(erosion_old.max()),
          'modified_max_lowering_m':float(erosion_new.max()),
          'additional_lowering_max_m':float(delta.max()),
          'additional_lowering_min_m':float(delta.min()),
          'original_positive_scour_volume_m3':float(np.dot(np.maximum(erosion_old,0),node_area)),
          'modified_positive_scour_volume_m3':float(np.dot(np.maximum(erosion_new,0),node_area)),
          'note':'All results calculated from supplied SELAFIN files. Positive scour volume is integrated positive bed lowering; it is not measured transported/exported sediment. First saved frame is the reference bed.'}
        (out/'07_summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
        (out/'README_results.txt').write_text('Actual SELAFIN comparison; no synthetic values.\n01: final bed maps\n02: bed lowering maps, same scale\n03: 1-p minus original lowering\n04: X-profile at Y of strongest-scour node\n05: time series using common saved frames only\n05_time_series_metrics.csv: time statistics\n06_nodewise_erosion.csv: node data\n07_summary.json: summary and provenance\nNote: Positive scour volume is not exported sediment mass.\n',encoding='utf-8')
        print('\nSUCCESS: comparison saved to',out)
        print(json.dumps({k:v for k,v in report.items() if k.endswith('_m') or k.endswith('_m3')},indent=2))
    finally:
        old.close();new.close()

if __name__=='__main__':
    try:main()
    except Exception as e:
        print('\nERROR:',type(e).__name__,str(e),file=sys.stderr)
        sys.exit(1)
