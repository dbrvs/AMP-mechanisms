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

### sensitivity analyses, run many trials ###
N=1000 #number of participants in each trial

#need the control trial
tr_ctl=AMPmodule.trial(name='ctl', nps=N, tF=560, dose=10, clade='all', 
                        A0=1, iv_phi=100000, maxIC50=50, v_flg=False, keep_flg=True)

#run local sensitivity
sens_iv=AMPmodule.sensitivity(N=N,dose=10,clade='all',A0=np.array([1]),iv=np.logspace(0,4,20),mx=np.array([50]))
sens_A0=AMPmodule.sensitivity(N=N,dose=10,clade='all',A0=np.round(np.logspace(0,2,10)),iv=1,mx=50)
sens_mx=AMPmodule.sensitivity(N=N,dose=10,clade='all',A0=1,iv=1,mx=np.logspace(1,5,20))


nts=len(sens)
trials=np.arange(nts)

med_ttd=np.zeros([nts]); std_ttd=np.zeros([nts])
med_fv=np.zeros([nts]); std_fv=np.zeros([nts])

brks4=np.zeros([nts]); brks8=np.zeros([nts]); clrs=np.zeros([nts]); occls=np.zeros([nts])

med_vlz4=np.zeros([nts]); med_vlz8=np.zeros([nts]); 
std_vlz4=np.zeros([nts]); std_vlz8=np.zeros([nts]); 

med_ic50=np.zeros([nts]); std_ic50=np.zeros([nts]); #breakthru ic50s
med_ic50f=np.zeros([nts]); std_ic50f=np.zeros([nts]); #breakthru ic50 factors


for it in trials:
    tr=sens[it]
    med_ttd[it]=np.median(tr.ttds)
    std_ttd[it]=np.std(tr.ttds)
    
    brks4[it]=np.sum(np.array(tr.fp_ints)==1)
    brks8[it]=np.sum(np.array(tr.fp_ints)==2)
    occls[it]=len(tr.occls)
    
    vlz=np.array(tr.fp_vls)
    
    med_fv[it]=np.median(np.log10(vlz[vlz>0]))
    std_fv[it]=np.std(np.log10(vlz[vlz>0]))

    med_vlz4[it]=np.median(np.log10(vlz[np.array(tr.fp_ints)==1]))
    med_vlz8[it]=np.median(np.log10(vlz[np.array(tr.fp_ints)==2]))

    std_vlz4[it]=np.std(np.log10(vlz[np.array(tr.fp_ints)==1]))
    std_vlz8[it]=np.std(np.log10(vlz[np.array(tr.fp_ints)==2]))
    
    med_ic50[it]=np.median(tr.bt_ic50s)
    std_ic50[it]=np.median(tr.bt_ic50s)

    med_ic50f[it]=np.median(tr.bt_ic50fs)
    std_ic50f[it]=np.median(tr.bt_ic50fs)


# In[171]:

# INPUTS
# In vivo : in vitro IC50 ratio (1-100)
# Proportion of isolates with IC50 >500 ug/mL (0-1)
# OUTPUTS
# Median time to detection
# Variance of time to detection
# # of episodes detected at the 4 week sample (# detected at week 4 / # detected at week 4 + # detected at week 8 + # that would only be detected after week 8)
# # of episodes detected at the 8 week sample (# detected at week 8 / # detected at week 4 + # detected at week 8 + # that would only be detected after week 8)
# # or episodes that may only be detected with serology (# that would only be detected after week 8 / # detected at week 4 + # detected at week 8 + # that would only be detected after week 8)
# Median log10 viral load at detection (this is as far as we?ll need to go with viral dynamics)
# Variance in log10 viral load at detection
# Median log10 viral load at detection among week 4 isolates (this is as far as we?ll need to go with viral dynamics)
# Median log10 viral load at detection among week 8 isolates (this is as far as we?ll need to go with viral dynamics)
# Variance in log10 viral load at detection among week 4 isolates (this is as far as we?ll need to go with viral dynamics)
# Variance in log10 viral load at detection among week 8 isolates (this is as far as we?ll need to go with viral dynamics)
# Median in vitro IC50 of breakthrough viruses
# Variance of in vitro IC50 of breakthrough viruses
# Median ratio of VRC01 at time of breakthrough : breakthrough IC50
# Variance in ratio of VRC01 at time of breakthrough : breakthrough IC50

#tr.exp_ints,  #exposure intervals
#tr.fp_ints,   #fp interval
#tr.fp_vls,    #first pos viral load
#tr.occls,     #occl infection times
#tr.brks,      #number of breakthrus
#tr.ttds,      #time to detection
#tr.bt_ic50s,  #ic50s of breakthru
#tr.bt_ic50fs  #ic50 ratio of breakthru


# In[183]:

#prevention efficacy plots

brks=np.zeros([nts])
phi=np.zeros([nts])
for it in trials:    
    phi[it]=sens[it].iv_phi
    brks[it]=len(sens[it].brks)
    
plt.plot(phi,1-brks/len(tr_ctl.brks),color='k',marker='D',ls='')
plt.xlim([0,50])
plt.ylim([0,1])
plt.ylabel('prevention efficacy, PE')
plt.xlabel('in vivo IC50 multiplier $\phi$')

plt.tight_layout()


# In[184]:

#time to detection plots
for it in trials:    
    ttds=np.array(sens[it].ttds)
    
    ttds=ttds[ttds>0]/7
    
    plt.subplot(121)
    s1=plt.scatter(np.zeros(len(ttds))+sens[it].iv_phi,ttds,alpha=0.01)
    s2=plt.plot(sens[it].iv_phi,np.median(ttds),color='k',marker='D',ls='')
    plt.xlim([0,50])
    plt.ylim([0,10])
    plt.ylabel('time to detection (weeks)')
    plt.xlabel('in vivo IC50 multiplier, $\phi$')
    
    plt.subplot(122)
    plt.plot(sens[it].iv_phi,np.std(ttds),color='k',marker='D',ls='')
    plt.xlim([0,100])
    plt.ylabel('std dev time to detection (weeks)')
    plt.xlabel('in vivo IC50 multiplier, $\phi$')

#plt.subplot(121)
#plt.legend([s1,s2],['individuals','median'])
plt.tight_layout()


# In[188]:

#viral load at detection plots
for it in trials:    
    vlz=np.array(sens[it].fp_vls)
    
    vlzA=np.log10(vlz[vlz>0])
    vlz4=np.log10(vlz[np.array(tr.fp_ints)==1])
    vlz8=np.log10(vlz[np.array(tr.fp_ints)==2])
    
    plt.subplot(131)
    plt.scatter(np.zeros(len(vlzA))+sens[it].iv_phi,vlzA,alpha=0.01)
    plt.plot(sens[it].iv_phi,np.mean(vlzA),color='k',marker='D',ls='')
    plt.xlim([0,50])
    #plt.ylim([0,10])
    plt.ylabel('viral load at first positive log10(copies/mL)')
    plt.xlabel('in vivo IC50 multiplier, $\phi$')
    plt.title('any time')

    plt.subplot(132)
    plt.scatter(np.zeros(len(vlz4))+sens[it].iv_phi,vlz4,alpha=0.01)
    plt.plot(sens[it].iv_phi,np.mean(vlz4[vlz4>0]),color='k',marker='D',ls='')
    plt.xlim([0,50])
    #plt.ylim([0,10])
    #plt.ylabel('viral load at first positive log10(copies/mL)')
    plt.xlabel('in vivo IC50 multiplier, $\phi$')
    plt.title('week 0-4')
    
    plt.subplot(133)
    plt.scatter(np.zeros(len(vlz8))+sens[it].iv_phi,vlz8,alpha=0.01)
    plt.plot(sens[it].iv_phi,np.mean(vlz8[vlz8>0]),color='k',marker='D',ls='')
    plt.xlim([0,50])
    #plt.ylim([0,10])
    #plt.ylabel('viral load at first positive log10(copies/mL)')
    plt.xlabel('in vivo IC50 multiplier, $\phi$')
    plt.title('week 4-8')

plt.tight_layout()


# In[189]:

plt.plot(trials,brks4)
plt.plot(trials,brks8)
#get clears!
#plt.plot(trials,tr.nps-brks4-brks8)
plt.ylabel('counts')
plt.xlabel('trial')
plt.legend(['week 0-4','week 4-8','occult'])
plt.tight_layout()


# In[192]:

plt.subplot(121)
plt.errorbar(trials,med_ic50,std_ic50)
plt.ylabel('IC50 of breakthrough strain ($\mu$g/mL)')
plt.xlabel('trial')

plt.subplot(122)
plt.errorbar(trials,med_ic50f,std_ic50f)
plt.ylabel('ratio of VRC01 concentration \n to IC50 of breakthrough strain')
plt.xlabel('trial')

plt.tight_layout()


# In[ ]:



