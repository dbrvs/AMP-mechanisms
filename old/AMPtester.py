#!/usr/bin/env python
#code the simulates AMP trials

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr,spearmanr
import AMPmodule

import resource
resource.setrlimit(resource.RLIMIT_NOFILE, (1000,-1)) #allow many plots

rez=600

#simple test for control ppt
p_ctl=AMPmodule.participant(name=1, et=50, tF=560, stop2=False, A0=1, vdp=VD[1,:], pkp=[1e-5,0,0,0], pdp=[1e5,1], v_flg=True)
t,sol=AMPmodule.p_ctl.vd_sim()
AMPmodule.plot_ppt('figures/ex_ctl',t,sol,p_ctl)

#simple tests for VRC01 ppt
VD=np.array(pd.DataFrame.from_csv('data/viral_dynamics.csv'))
PK=np.array(pd.DataFrame.from_csv('data/PK10.csv'))
PD=np.array(pd.DataFrame.from_csv('data/PDall.csv'))
p_VRC=AMPmodule.participant(name=1, et=50, tF=560, stop2=False, A0=1, vdp=VD[1,:], pkp=PK[1,:], pdp=PD[1,:], v_flg=True)
t,sol=AMPmodule.p_VRC.vd_sim()
AMPmodule.plot_ppt('figures/ex_VRC01',t,sol,p_VRC)

#test all viral dynamics parameter sets in control setting in a single dose interval
VD=np.array(pd.DataFrame.from_csv('data/viral_dynamics.csv'))
etf=0 #how to space out for plotting
for i in range(len(VD)):
    p=AMPmodule.participant(name=1, et=i*etf, tF=8*7+1, stop2=False, A0=1, 
                  vdp=VD[i,:], pkp=[1e-5,0,0,0], pdp=[1e5,1], v_flg=False)
    t,sol=AMPmodule.p.vd_sim()
    plt.semilogy(t/7,sol[:,6]/p.vol*1e3)
    plt.scatter(p.obs_times/7,p.obs+0.7**i)
plt.xlabel('time (weeks)')
plt.ylabel('viral load (copies/mL)')
plt.ylim([1,1e8])
plt.yticks(np.logspace(0,8,9))
plt.axhline(30,color='k')
plt.tight_layout()
plt.xlim([-0.5,8.5])
plt.tight_layout()
plt.savefig('figures/ex_allsets.pdf')

#simple test of trial including VRC01
tr=trial(name='test', nps=10, tF=560, dose=10, clade='B', A0=1, iv_phi=1, maxIC50=50, v_flg=True, keep_flg=True)
plt.figure(figsize=(5,3),dpi=rez)
for ip in range(tr.nps):
    plt.scatter(obs_times,np.log10(tr.obz[ip]+2e-4*1.1**ip))#,c=plt.cm.jet(ip))#1.4**ip),alpha=0.7)
plt.axhline(np.log10(30),color='k',ls='--')
plt.ylim([-4,8])
plt.ylabel('viral load log10(copies/mL)')
plt.xlabel('time (weeks)')
plt.title('breakthroughs: '+str(len(tr.brks))+'/'+str(tr.nps))
plt.tight_layout()
plt.savefig('figures/ex_trial.pdf')

#boxplot to check on parameters
plt.figure(figsize=(7,3),dpi=rez)
plt.subplot(131)
plt.boxplot(tr.VD)
plt.semilogy()
plt.xticks(np.arange(7)+1,
           [r'$\mathcal{R}_0$',r'$\delta_S$',r'$\tau$',r'$\kappa$',r'$\delta_E$',r'$E_{50}$',r'$\omega$'])
plt.subplot(132)
plt.boxplot(tr.PK)
plt.semilogy()
plt.xticks(np.arange(4)+1,[r'$\mathcal{Y}_1$',r'$k_1$',r'$\mathcal{Y}_2$',r'$k_2$'])
plt.title('dose: ' + str(tr.dose))
plt.subplot(133)
plt.boxplot(tr.PD)
plt.semilogy()
plt.xticks(np.arange(2)+1,[r'$IC_{50}$',r'$h$'])
plt.title('clade: ' + str(tr.clade))
plt.tight_layout()
plt.savefig('figures/param_boxes.pdf')

