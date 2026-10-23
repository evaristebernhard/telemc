# -*- coding: utf-8 -*-
"""Stage 17: compare 10 completed 100-s B-2 runs with Fig.5(b) stations.
Exploratory cross-experiment comparison only, NOT calibration/validation.
Run from configured openTELEMAC PowerShell:
    python .\07_scripts\15_compare_17_tauc_kd_ABC.py
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from data_manip.extraction.telemac_file import TelemacFile

ROOT = Path(__file__).resolve().parent.parent
RUN, OUT, FIG = ROOT/'05_run', ROOT/'06_output', ROOT/'08_figures'
OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
RAW = {
 'A': (0.15, [(-0.03348,79.6714),(9.99841,79.6714),(19.90412,78.21722),(30.0622,77.20669),(40.22028,72.84416),(50.37837,70.97098),(60.41026,66.90421),(70.56834,62.83744),(80.24271,57.90803),(90.2746,52.95396),(100.1803,51.22867)]),
 'B': (0.65, [(-0.15967,72.99204),(9.99841,71.68574),(20.15649,69.9358),(30.18839,68.92527),(40.34647,66.90421),(50.37837,61.82691),(60.41026,58.32703)]),
 'C': (1.15, [(-0.03348,62.83744),(10.1246,62.54168),(20.15649,61.67903),(30.18839,59.7812),(40.22028,55.71444),(50.25218,52.95396),(60.28407,50.21814),(70.21081,52.68284),(80.24271,49.62661),(90.40079,48.46819),(100.3065,48.61608)])
}
BASE = {'KD050':'B2_Tingsanchali_erosion_kd_0p5e4_res.slf',
        'KD100':'B2_Tingsanchali_erosion_kd_1p0e4_res.slf',
        'KD250':'B2_Tingsanchali_erosion_kd_2p5e4_res.slf',
        'KD500':'B2_Tingsanchali_erosion_ref_res.slf',
        'KD800':'B2_Tingsanchali_erosion_kd_8p0e4_res.slf'}
KD = {'KD050':0.5e-4,'KD100':1e-4,'KD250':2.5e-4,'KD500':5e-4,'KD800':8e-4}
TAUS = [(0.01,BASE), (0.1,{k:f'B2_100s_TC0100_{k}_res.slf' for k in KD})]
COLORS = dict(zip(KD,plt.get_cmap('tab10').colors[:5]))

def main():
    expected=[RUN/f for _,files in TAUS for f in files.values()]
    missing=[str(f) for f in expected if not f.is_file()]
    if missing: raise FileNotFoundError('Missing SLF:\n'+'\n'.join(missing))
    sim={}; diagnostics=[]
    for tau,files in TAUS:
        for case,filename in files.items():
            print(f'Reading tau_c={tau:g} {case}: {filename}',flush=True)
            res=TelemacFile(str(RUN/filename))
            try:
                points,s,values=res.get_timeseries_on_polyline('BOTTOM',[[0.,0.5],[1.2,0.5]],[24])
                s=np.asarray(s,float); z=np.asarray(values,float); t=np.asarray(res.times,float)
                if z.shape!=(len(s),len(t)): raise ValueError(f'BOTTOM shape {z.shape} != {(len(s),len(t))}')
                if len(t)<2 or np.any(np.diff(t)<=0): raise ValueError('Invalid output times')
                if t[-1]<99.99: raise ValueError(f'Not a 100-s result: {filename}, end={t[-1]}')
                if not np.isfinite(z).all(): raise ValueError(f'Nonfinite BOTTOM: {filename}')
                data={}
                for station,(x,_) in RAW.items():
                    idx=int(np.argmin(np.abs(s-x)))
                    if abs(s[idx]-x)>1e-5: raise ValueError(f'Station {station} not on polyline')
                    data[station]=z[idx].copy()
                sim[(tau,case)]=(t,data)
                # These checks concern the three stations, not the whole 2D mesh.
                for station,series in data.items():
                    dz=np.diff(series)
                    diagnostics.append([tau,case,station,float(t[-1]),float(series[0]),float(series[-1]),
                                        float(series.min()),float(series.max()),float(np.min(dz)),
                                        int(np.sum(dz>1e-6)),int(np.sum(series<-1e-6))])
            finally: res.close()
    with (OUT/'B2_100s_tauc_kd_ABC_station_diagnostics.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['tau_c_Pa','case','station','end_time_s','initial_bottom_m','final_bottom_m','min_bottom_m','max_bottom_m','min_step_change_m','number_of_rising_output_steps','number_of_negative_bottom_outputs']);w.writerows(diagnostics)
    metrics=[]
    for station,(x,raw) in RAW.items():
        obs=np.asarray(raw,float); ot=obs[:,0]; oz=obs[:,1]/100
        for tau,_ in TAUS:
            for case in KD:
                t,dat=sim[(tau,case)]; z=dat[station]
                valid=(ot>=t[0])&(ot<=t[-1])
                residual=np.interp(ot[valid],t,z)-oz[valid]
                metrics.append([station,tau,case,KD[case],len(residual),
                    float(np.sqrt(np.mean(residual**2))),float(np.mean(np.abs(residual))),float(np.mean(residual))])
    with (OUT/'B2_100s_tauc_kd_ABC_exploratory_metrics.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['station','tau_c_Pa','case','kd_m_per_Pa_s','n','RMSE_m','MAE_m','bias_sim_minus_measured_m']);w.writerows(metrics)
    # Three rows (stations) x two columns (tau_c); keep original digitized times.
    fig,axs=plt.subplots(3,2,figsize=(14,13),sharex=True,sharey=False,constrained_layout=True)
    for i,(station,(x,raw)) in enumerate(RAW.items()):
        obs=np.asarray(raw,float)
        for j,(tau,_) in enumerate(TAUS):
            ax=axs[i,j]
            ax.scatter(np.maximum(obs[:,0],0),obs[:,1]/100,color='black',s=26,zorder=10,label='Measured (Fig.5b)')
            for case in KD:
                t,data=sim[(tau,case)]
                ax.plot(t,data[station],label=case,color=COLORS[case],lw=1.55)
            ax.set_xlim(0,100);ax.grid(alpha=.22)
            ax.set_title(f'Station {station}, x={x:.2f} m | tau_c={tau:g} Pa')
            ax.set_ylabel('Bed elevation (m)')
            if i==2:ax.set_xlabel('Time (s)')
            if i==0:ax.legend(fontsize=8,ncol=2)
    fig.suptitle('100-s B-2 simulations vs Fig.5(b) measured (different experimental forcing)')
    fig.savefig(FIG/'B2_100s_tauc_kd_ABC_6panel.png',dpi=190);plt.close(fig)
    # Five figures isolate the effect of tau_c at each k_d.
    for case in KD:
        fig,axs=plt.subplots(1,3,figsize=(16,4.8),constrained_layout=True)
        for ax,(station,(x,raw)) in zip(axs,RAW.items()):
            obs=np.asarray(raw,float)
            ax.scatter(np.maximum(obs[:,0],0),obs[:,1]/100,color='black',s=24,label='Measured (different experiment)',zorder=10)
            for tau,style in [(0.01,'-'),(0.1,'--')]:
                t,data=sim[(tau,case)]
                ax.plot(t,data[station],style,lw=1.8,label=f'tau_c={tau:g} Pa')
            ax.set(title=f'Station {station} (x={x:.2f} m)',xlabel='Time (s)',ylabel='Bed elevation (m)',xlim=(0,100))
            ax.grid(alpha=.22);ax.legend(fontsize=8)
        fig.suptitle(f'{case}: effect of tau_c (100 s, cross-experiment exploratory comparison)')
        fig.savefig(FIG/f'B2_100s_tauc_effect_{case}_ABC.png',dpi=180);plt.close(fig)
    print('\nStation  tau_c    Case      N     RMSE(m)    MAE(m)     Bias(m)')
    for station,tau,case,kd,n,rmse,mae,bias in metrics:
        print(f'{station:7s} {tau:6.2f}  {case:6s} {n:3d}  {rmse:10.6f} {mae:10.6f} {bias:10.6f}')
    suspicious=[d for d in diagnostics if d[-2] or d[-1]]
    print(f'\nStation diagnostics: {len(suspicious)} case/station entries have rising steps or negative bottom values; inspect CSV.')
    print('Saved 6 PNG in 08_figures and 2 new CSV in 06_output.')
    print('NOTE: station checks do NOT establish whole-domain physical validity.')
    print('CAUTION: different experimental forcing; errors are NOT calibration/validation scores.')

if __name__=='__main__':main()
