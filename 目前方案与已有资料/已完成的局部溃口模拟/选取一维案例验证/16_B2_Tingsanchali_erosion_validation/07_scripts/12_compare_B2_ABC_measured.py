# -*- coding: utf-8 -*-
"""Exploratory comparison only: Fig.5(b) calibration observations vs B-2 simulations.
Run with: python .\07_scripts\12_compare_B2_ABC_measured.py
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from data_manip.extraction.telemac_file import TelemacFile

ROOT = Path(r'F:\Telemac_Projects\Tuanzho_New\16_B2_Tingsanchali_erosion_validation')
RUN = ROOT / '05_run'
FIG = ROOT / '08_figures'
OUT = ROOT / '06_output'
FIG.mkdir(exist_ok=True, parents=True)
OUT.mkdir(exist_ok=True, parents=True)

# Raw digitization: seconds and centimetres, retained without smoothing.
RAW = {
 'A': (0.15, [(-0.03348,79.6714),(9.99841,79.6714),(19.90412,78.21722),(30.0622,77.20669),(40.22028,72.84416),(50.37837,70.97098),(60.41026,66.90421),(70.56834,62.83744),(80.24271,57.90803),(90.2746,52.95396),(100.1803,51.22867)]),
 'B': (0.65, [(-0.15967,72.99204),(9.99841,71.68574),(20.15649,69.9358),(30.18839,68.92527),(40.34647,66.90421),(50.37837,61.82691),(60.41026,58.32703)]),
 'C': (1.15, [(-0.03348,62.83744),(10.1246,62.54168),(20.15649,61.67903),(30.18839,59.7812),(40.22028,55.71444),(50.25218,52.95396),(60.28407,50.21814),(70.21081,52.68284),(80.24271,49.62661),(90.40079,48.46819),(100.3065,48.61608)])
}
CASES = [
 ('KD050', 0.5e-4, 'B2_Tingsanchali_erosion_kd_0p5e4_res.slf'),
 ('KD100', 1.0e-4, 'B2_Tingsanchali_erosion_kd_1p0e4_res.slf'),
 ('KD250', 2.5e-4, 'B2_Tingsanchali_erosion_kd_2p5e4_res.slf'),
 ('KD500', 5.0e-4, 'B2_Tingsanchali_erosion_ref_res.slf'),
 ('KD800', 8.0e-4, 'B2_Tingsanchali_erosion_kd_8p0e4_res.slf'),
]

def main():
    missing = [str(RUN / f) for _,_,f in CASES if not (RUN / f).is_file()]
    if missing:
        raise FileNotFoundError('Missing SLF files:\n' + '\n'.join(missing))
    # Preserve experimental data independently of the comparison.
    with (OUT / 'B2_Fig5b_ABC_measured_digitized.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f)
        writer.writerow(['station','x_m','time_s','bed_elevation_m','source'])
        for name,(x,pts) in RAW.items():
            for t,z_cm in pts:
                writer.writerow([name,x,t,z_cm/100,'Fig.5(b) Measured (calibration case, not B-2)'])
    sim={}
    for code,kd,filename in CASES:
        print('Reading',code,filename)
        res=TelemacFile(str(RUN/filename))
        try:
            # TELEMAC 9.1 returns (points, curvilinear_abscissa, values).
            # 24 segments -> exactly 0.05 m increments along the centreline.
            points,s,values=res.get_timeseries_on_polyline(
                'BOTTOM', [[0.0,0.5],[1.2,0.5]], [24])
            s=np.asarray(s,dtype=float)
            values=np.asarray(values,dtype=float)
            times=np.asarray(res.times,dtype=float)
            if values.shape != (len(s),len(times)):
                raise RuntimeError(f'Unexpected polyline array shape {values.shape}')
            sim[code]={}
            for name,(x,_) in RAW.items():
                idx=int(np.argmin(abs(s-x)))
                if abs(s[idx]-x)>1e-5:
                    raise RuntimeError(f'Station {name} mismatch: {s[idx]} vs {x}')
                sim[code][name]=(times,values[idx,:])
        finally:
            res.close()
    # Compare at experimental timestamps in the common simulation interval.
    records=[]
    for name,(x,raw) in RAW.items():
        exp=np.array(raw,dtype=float)
        te=exp[:,0]
        ze=exp[:,1]/100.0
        for code,kd,_ in CASES:
            ts,zs=sim[code][name]
            valid=(te>=ts.min()) & (te<=ts.max())
            # Fig.5(b) first point is slightly negative due to digitization.
            # Treat it as nominal t=0 only for plotting; exclude it from metrics.
            n=int(valid.sum())
            if n:
                pred=np.interp(te[valid],ts,zs)
                diff=pred-ze[valid]
                rmse=float(np.sqrt(np.mean(diff**2)))
                mae=float(np.mean(abs(diff)))
                bias=float(np.mean(diff))
            else:
                rmse=mae=bias=float('nan')
            records.append([name,code,kd,n,rmse,mae,bias])
    with (OUT/'B2_Fig5b_ABC_exploratory_metrics.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f)
        w.writerow(['station','case','kd_m_per_Pa_s','n_points','RMSE_m','MAE_m','bias_sim_minus_obs_m'])
        w.writerows(records)
    for name,(x,raw) in RAW.items():
        exp=np.array(raw,dtype=float)
        fig,ax=plt.subplots(figsize=(8.2,5.4))
        # Clip slightly negative digitized time to zero for display only.
        ax.scatter(np.maximum(exp[:,0],0),exp[:,1]/100,color='black',marker='o',s=32,
                   label='Fig. 5(b) Measured (different experiment)',zorder=6)
        for code,kd,_ in CASES:
            ts,zs=sim[code][name]
            ax.plot(ts,zs,lw=1.8,label=f'{code}: kd={kd:.1e}')
        ax.set(xlim=(0,30),xlabel='Time (s)',ylabel='Bed elevation (m)',
               title=f'Station {name} (x={x:.2f} m): exploratory comparison')
        ax.grid(alpha=.25);ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG/f'B2_Fig5b_station_{name}_exploratory.png',dpi=250)
        plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),sharex=True)
    for ax,(name,(x,raw)) in zip(axes,RAW.items()):
        exp=np.array(raw,dtype=float)
        ax.scatter(np.maximum(exp[:,0],0),exp[:,1]/100,color='black',s=27,zorder=6,label='Measured')
        for code,kd,_ in CASES:
            ts,zs=sim[code][name]
            ax.plot(ts,zs,lw=1.5,label=code)
        ax.set(title=f'Station {name}: x={x:.2f} m',xlabel='Time (s)',xlim=(0,30))
        ax.grid(alpha=.25)
    axes[0].set_ylabel('Bed elevation (m)')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=6,bbox_to_anchor=(.5,-.03))
    fig.suptitle('Exploratory only: Fig. 5(b) calibration data vs B-2 simulation')
    fig.tight_layout(rect=[0,.06,1,.94])
    fig.savefig(FIG/'B2_Fig5b_ABC_exploratory_combined.png',dpi=250,bbox_inches='tight')
    plt.close(fig)
    print('\nStation  Case      N      RMSE(m)     MAE(m)      Bias(m)')
    for name,code,kd,n,rmse,mae,bias in records:
        print(f'{name:7s}  {code:6s} {n:3d}  {rmse:11.6f} {mae:11.6f} {bias:11.6f}')
    print('\nSaved: 4 PNG in 08_figures; 2 CSV in 06_output.')
    print('CAUTION: cross-experiment exploratory errors are NOT calibration/validation scores.')

if __name__=='__main__':
    main()
