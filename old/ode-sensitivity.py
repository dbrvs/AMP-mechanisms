
# coding: utf-8

# In[1]:

#!/usr/bin/env python
#get_ipython().magic('matplotlib inline')

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import pandas as pd
import pyDOE
from scipy.stats import pearsonr,spearmanr
from scipy.interpolate import interp1d

rez=600


# In[2]:

#whole AB model
def Ab_model(t,A0,phi,r,s,hill,IC50):
    Ab_t = A0*(np.exp(-r*t)+phi*np.exp(-s*t))
    if A0>0:
        E_t = 1/(1+(Ab_t/IC50)**-hill)
    else:
        E_t=0
    return Ab_t,E_t

def model(X,t,aS,dS,Bt0,tau,lam,dI,thL,xi,k,aE,dE,E50,w,p,g,A0,phi,r,s,hill,IC50):
    
    dY = np.zeros(7);
    Ab_t,E_t = Ab_model(t,A0,phi,r,s,hill,IC50)    

    S=X[0]; AU=X[1]; AP=X[2]; LU=X[3]; LP=X[4]; E=X[5]; V=X[6];    
    
    Bt=Bt0*(1-E_t)
    
    dY[0] = aS - dS*S - Bt*S*V                         #susceptible cells
    dY[1] = (1-tau)*(1-lam)*Bt*S*V - dI*AU - k*E*AU + xi*LU    #active unproductively infected
    dY[2] = tau*(1-lam)*Bt*S*V - dI*AP - k*E*AP + xi*LP        #active productively infected
    dY[3] = (1-tau)*lam*Bt*S*V + thL*LU                #latent unproductively infected
    dY[4] = tau*lam*Bt*S*V + thL*LP                    #latent productively infected
    dY[5] = w*E*(AP+AU)/(E+E50) + aE - dE*E;           #adaptive immune system
    dY[6] = p*AP - g*V - Bt*S*V                        #virus
    return dY

def run_model(tt,X0,aS,dS,Bt,tau,lam,dI,thL,xi,k,aE,dE,E50,w,p,g,A0,phi,r,s,hill,IC50):
    
    infpoz,maxpoz,t_fp,t_max=0,0,0,0
    
    sol=spi.odeint(model, X0, tt, 
                   (aS,dS,Bt,tau,lam,dI,thL,xi,k,aE,dE,E50,w,p,g,A0,phi,r,s,hill,IC50),
                   mxstep=10000)    

    vll=np.log10(sol[:,6]/vol*1e3) #viral load as usual
    
    
    if (vll>-10).all() & (vll>2).any():        
        infpoz=np.where(vll>2)[0][0] #index of first positive
        maxpoz=np.where(vll==max(vll))[0][0]
        t_fp=tt[infpoz]
        t_max=tt[maxpoz]
    
    return sol,vll,infpoz,t_fp,t_max


# In[3]:

#fit Robb VL
robb=pd.read_csv('data/robb_ranges.csv',names=['logVL','days']); robbHI=robb[:15]; robbLO=robb[15:]

interpLO=interp1d(robbLO.days,robbLO.logVL); interpHI=interp1d(robbHI.days,robbHI.logVL)


# In[4]:

# parameters for viral dynamics model
vol = 1      # volume of blood [uL]
aS  = 70*vol;   #constant growth rate of susceptibles [cells/uL/day]
dS  = 0.3;   #susceptible death rate [1/day] 
Bt0 = 1e-4/vol  # infection rate of T-cells [uL/cells-day]/[uL]
dA  = 1.0       # active death rate [1/day]
p   = 5e4       # burst rate of virus from cells [virions/cell]
g   = 23        # virus clearance rate [1/day]
tau = 0.05      # productive infection probability []
lam = 1e-4      # latency probability []
thL = 5.2e-4    # latent clearance rate [1/day]
aL  = 0.015;    # latent proliferation rate [1/day] (Tcm)
xi  = 1e-5;     # latent activation rate [1/day]
dL  = aL-thL-xi # latent death rate
k   = 0.3/vol;  #immune cell killing rate [uL/cell-day]/[uL]
w   = 1.6;     #immune cell multiplier [1/day]
aE  = 1e-4*vol;   #initial E cell concentration [cells/uL]*[uL]
dE  = 0.002;  #immune death rate [1/day]
E50 = 250*vol;   #50 pct max E cell concentration [cells/uL]*[uL]

R0=aS*Bt0*tau*(1-lam)*p/g/dS/dA; print(R0) # basic reproductive number

#PK parameters
A0,phi,r,s = 1, 0, 0.348, 0.0538 #low 30 or high 10 dose range

#PD parameters
IC50=20 # IC50 [ug/uL]
hill=1  # hill coefficient []

V0=1e-3 #1 virion in body!

X0=np.array([aS/dS,0,0,0,0,aE/dE,V0])

tY=np.linspace(0,7*8,1e3) #time for antibody dose
Ab_t,E_t  = Ab_model(tY,A0,phi,r,s,hill,IC50) #Ab model

tc=7 #first contact at 1 week

in_tc=np.where(tY>tc)[0][0] #index of first contact

tV=tY[in_tc:] #time for virus

sol,vll,infpoz,t_fp,t_max=run_model(
    tV,X0,aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g,A0,phi,r,s,hill,IC50) #first positive at 1 week

#calculate R0(t)
#dIdt=np.diff(sol[:,2])
#I=sol[1:,2]
#dt=tt[1]
#Rt=1+dIdt/(I*dt*dA)
Rt=Bt0*tau*(1-lam)*(1-E_t)/dA*p/g*aS/dS

fig,axarr=plt.subplots(3,2,figsize=(6,5),dpi=rez,sharex=True)
titz=['S','$A_U$','$A_P$','$L_U$','$L_P$','$E$','$V$']

axarr[0][0].plot(tV/7,sol[:,0],color='gray')
axarr[0][0].set_ylabel('cells/$\mu$L')
axarr[0][0].legend(titz[0],fontsize=8,loc=2)

axarr[1][0].plot(tV/7,sol[:,2],color='red',ls='--')
axarr[1][0].plot(tV/7,sol[:,1],color='red')
axarr[1][0].set_ylabel('cells/$\mu$L')
axarr[1][0].legend(titz[1:3],fontsize=8,loc=2)

axarr[2][0].plot(tV/7,sol[:,3],color='blue',ls='--')
axarr[2][0].plot(tV/7,sol[:,4],color='blue')
axarr[2][0].legend(titz[3:],fontsize=8,loc=2)
axarr[2][0].set_ylabel('cells/$\mu$L')

axarr[0][0].set_xticks(range(9))
axarr[2][0].set_xlabel('time (weeks)')

#plot viral stuff
axarr[0][1].plot(tY/7,Ab_t,color='gold')
axarr[0][1].set_ylabel('bnAbs ($\mu$g/$\mu$L)')
    
axarr[1][1].plot(tY,Rt,color='black') #allow burn in for some R0
axarr[1][1].axhline(1,color='k',ls='--')
axarr[1][1].set_xticks(range(9))
axarr[1][1].set_ylabel(r'$\mathcal{R}(t)$')
#axarr[1].set_yticks(np.linspace(0,9,4))

axarr[2][1].plot(tV/7,vll,color='purple') #virus/mL
#axarr[2][1].plot((t_fp_shift-tc/7,vll_fp_shift,color='purple') #virus/mL
#axarr[2][1].fill_between((tY+t_fp)/7, interpLO(tY), interpHI(tY),color='gray',alpha=0.3)
#axarr[2][1].fill_between(tY/7, interpLO(tY), interpHI(tY),color='gray',alpha=0.3)
#axarr[2][1].set_ylim([np.log10(30),9])
axarr[2][1].set_xlim([0,8])
#axarr[2][1].set_yticks(np.linspace(2,8,7))
axarr[2][1].set_ylabel('viral load \n $\log_{10}$(copies/mL)')
axarr[2][1].set_xlabel('time (weeks)')

plt.tight_layout()
plt.savefig('figures/det_model.pdf')


# In[5]:

#now see which parameter ranges fit in here

sampz=10**4

tt=np.linspace(0,100,1e4)

guessP=np.array([aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g])

goodP_list=[]

bdz=1.#extra bounds for fitting shaded area

fposT=7

for i in range(sampz):
    aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr=np.random.lognormal(0,2,[15])*guessP
    
    if aSr*Bt0r*taur*(1-lamr)*pr/gr/dSr/dAr<10: # basic reproductive number

        sol,vll,infpoz,t_fp,t_max=run_model(
            tt,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)
        
        if (vll>-10).all() & (vll>2).any():        

            t_shifted = tt[infpoz:]-tt[infpoz]
            
            loV,hiV = interpLO(t_shifted),interpHI(t_shifted) #the bounds

            #first positive must be before 3 weeks and t_max must be after 3 days since first positive
            if (t_fp<21.) & (t_max-t_fp>3.):

                if ~ ((vll[infpoz:]>hiV*bdz).any() | (vll[infpoz:]<loV/bdz).any()):
                    goodP_list.append([aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr])#,x0r])
                    
                    #plt.plot(tt[infpoz:],vll[infpoz:])


# In[6]:

plt.figure(figsize=(4,3),dpi=rez)

loV,hiV = interpLO(tt),interpHI(tt) #the bounds
for i in range(len(goodP_list)):
    aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr=goodP_list[i]

    sol,vll,infpoz,t_fp,t_max=run_model(
        tt,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)

    plt.plot(tt[infpoz:]-tt[infpoz],vll[infpoz:])#,color='k',alpha=0.5)

plt.fill_between(tt, loV, hiV, color='gray',alpha=0.3)
plt.xlabel('time (days after first positive test)')
plt.ylabel('viral load $\log_{10}$(copies/mL)')
plt.tight_layout()
plt.savefig('figures/robb_examples.pdf')


# In[7]:

#boxplots from robb

#days to nadir
[17.878483268976264, 27.917385628280154, 30.934789888203735, 40.682980986119745, 41.7276975706514]

#days to EIA
[7.950909781228649,12.88320520610375,14.101772311072896, 18.04749340369393, 32.84449492559828]

#days to peak
[6.147116993275322, 6.320967513712819, 9.918987565458853, 13.052445835379505, 14.096931925353058, 17.868744873897068]

#peak log10 VL
[4.565059970202178, 6.35257087867727, 6.849094572557447, 7.663367820014483, 8.59685797501567]

#nadir log10 VL
[1.780762540024421, 3.6079697335034737, 4.422294201973418, 4.9188435063600515, 6.368654276730486]



# In[8]:

pnames=[r'$\alpha_S$',r'$\delta_S$',r'$\beta_0$',r'$\tau$',r'$\lambda$',
        r'$\delta_I$',r'$\theta_L$',r'$\xi$',r'$\kappa$',r'$\alpha_E$',
        r'$\delta_E$',r'$E_{50}$',r'$\omega$',r'$\pi$',r'$\gamma$']
        
fig,axarr=plt.subplots(3,5,sharey=False,figsize=(7,6),dpi=rez)
for j in range(15):
    ax=axarr[j%3][int(j/3.)]
    
    logpz=np.log10(np.array(goodP_list))
    
    #ax.violinplot(np.array(goodP_list)[:,j])#/guessP[j])
    ax.scatter(np.random.normal(1,0.05,len(goodP_list)),logpz[:,j],color='gray',alpha=0.5)#/guessP[j])
    ax.boxplot(logpz[:,j])#/guessP[j])
    ax.set_xlabel(pnames[j])
    ax.set_xticks([])
    ax.set_xlim([0.8,1.2])
    #ax.set_yscale('log')
#axarr[1][0].set_ylabel('deviation from \n median value')
plt.tight_layout()
plt.savefig('figures/robb_fitparams.pdf')


# In[9]:

pnames_unix=['aS','dS','Bt0','tau','lam','dA','thL','xi','k','aE','dE','E50','w','p','g']

gf=pd.DataFrame(goodP_list,columns=pnames_unix)

plt.figure(figsize=(7,3),dpi=rez)
plt.subplot(141)
plt.boxplot(gf.aS/gf.dS)
plt.yscale('log')
plt.xticks([1],[r'$\alpha_S/\delta_S$'])

plt.subplot(142)
plt.boxplot(gf.aE/gf.dE)
plt.yscale('log')
plt.xticks([1],[r'$\alpha_E/\delta_E$'])

plt.subplot(143)
plt.boxplot(gf.p/gf.g)
plt.yscale('log')
plt.xticks([1],[r'$\pi/\gamma$'])

plt.subplot(144)
plt.boxplot(gf.tau*(1-gf.lam)*gf.aS/gf.dS*gf.p/gf.g*gf.Bt0/gf.dA)
#ax.scatter(np.random.normal(1,0.1,len(gf.tau)),gf.tau*(1-gf.lam)*gf.aS/gf.dS*gf.p/gf.g*gf.Bt0/gf.dA)#/guessP[j])
#plt.yscale('log')
plt.xticks([1],[r'$\mathcal{R}_0$'])

plt.tight_layout()
plt.savefig('figures/robb_fitparams_identifiability.pdf')


# In[10]:

fig,axarr=plt.subplots(15,15,sharey=True,sharex=True,figsize=(10,10),dpi=rez)

cz_range=np.linspace(0,1,11)

for i in range(15):
    for j in range(15):
        ax=axarr[i][j]

        #only plot bottom diagonal
        if j>i-1:
            ax.axis('off')
        else:
            pcorr=spearmanr(np.log10(gf.iloc[:,i]/guessP[i]),np.log10(gf.iloc[:,j]/guessP[j]))
            #if pcorr[1]<0.05:
            #    cz='crimson'
            #else:
            #    cz='green'
            ax.scatter(np.log10(gf.iloc[:,i]/guessP[i]),np.log10(gf.iloc[:,j]/guessP[j]),s=10, color=[0,0,pcorr[1]],alpha=0.3)

    axarr[i][0].set_ylabel(pnames[i], fontsize=12)
    axarr[14][i].set_xlabel(pnames[i], fontsize=12)

ax.set_ylim([-2,2])
ax.set_xlim([-2,2])
ax.set_xticks([-1,0,1])
ax.set_yticks([-1,0,1])
ax.set_xticklabels([-1,0,1],fontsize=8)
ax.set_yticklabels([-1,0,1],fontsize=8)
fig.subplots_adjust(wspace=0.01, hspace=0.01)
fig.savefig('figures/robb_fitparams_varcorrs.pdf')


# In[11]:

################################################################################
#local sensitivity analysis function

#variables to check sensitivity
variables = ['first contact time ' + r'$t_c$',
             'bnAb conc phase 1 ' + r'$\mathcal{Y}_1$', 
             'bnAb conc phase 2 ' + r'$\mathcal{Y}_2$',
             'decay rate 1 ' + r'$k_1$',
             'decay rate 2 ' + r'$k_2$',
             'Hill coefficient ' + r'$h$',
             '50% effectiveness ' + r'$IC_{50}$']

#from data_analysis fitting:
#[86  2.5 0.15 0.007] +/- [ 32.7 32.6 0.12 0.29 ]
#[240 2.9 0.086 -0.011] +/- [ 112 116 0.081 0.71]

#ranges for best (b) and worst (w) cases
b_case=np.array([0 ,600,300,0  ,0  ,5  ,0.2])
w_case=np.array([30,0  ,0  ,0.3,0.7,0.8,20 ])
                            
tY=np.linspace(0,7*8,1e4) #time for antibody dose 8 week time series in days

fig1,axarr1 = plt.subplots(2,4,sharex=True,sharey=True,figsize=(9,4),dpi=rez)
fig2,axarr2 = plt.subplots(2,4,sharex=True,sharey=True,figsize=(9,4),dpi=rez)
fig3,axarr3 = plt.subplots(2,4,sharex=True,sharey=True,figsize=(9,4),dpi=rez)
#looping over the number of variables to check sensitivity

num_ppts=25
for k in range(num_ppts): 
    aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr=goodP_list[k] #choose a random model
    for i in range(len(variables)):
        ax=axarr1[i%2][int(i/2.)]
        axL=axarr2[i%2][int(i/2.)]
        axE=axarr3[i%2][int(i/2.)]

        lo=np.array([7,200,50,0.15,0.007,1,1]); hi=np.array([7,200,50,0.15,0.007,1,1]) #set both at typical values

        #change the i-th parameter
        lo[i]=w_case[i]; hi[i]=b_case[i]

        #evaluate model twice
        tc,A0,phi,r,s,hill,IC50=lo #choose low range of params
        in_tc=np.where(tY>tc)[0][0] #index of first contact
        tV=tY[in_tc:] #time for virus
        sol,vll,infpoz,t_fp,t_max=run_model(
            tV,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)

        loV=vll #for logplot later
        ax.plot(tV/7,loV,alpha=0.7,color='darkred')
        
        axL.plot(tV/7,sol[:,3],color='darkred',ls='--')
        axL.plot(tV/7,sol[:,4],color='darkred')
        
        axE.plot(tV/7,sol[:,5],color='darkred')

        tc,A0,phi,r,s,hill,IC50=hi #choose high range of params
        in_tc=np.where(tY>tc)[0][0] #index of first contact
        tV=tY[in_tc:] #time for virus
        sol,vll,infpoz,t_fp,t_max=run_model(
            tV,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)

        hiV=vll
        ax.plot(tV/7,hiV,alpha=0.7,color='deepskyblue')
        axL.plot(tV/7,sol[:,3]*1e6/vol,color='deepskyblue',ls='--')
        axL.plot(tV/7,sol[:,4]*1e6/vol,color='deepskyblue')
        
        axE.plot(tV/7,sol[:,5],color='deepskyblue')

        #ax.fill_between(lot/7,lo_sol[:,6]/vol*1e3,hi_sol[:,6]/vol*1e3,alpha=0.6)

        #ax.set_ylim([0.1,100])
        ax.set_title(variables[i],fontsize=10)# + str(rangez[i][0]) + r'$\to$' +str(rangez[i][1]),fontsize=10)
        axL.set_title(variables[i],fontsize=10)# + str(rangez[i][0]) + r'$\to$' +str(rangez[i][1]),fontsize=10)
        axE.set_title(variables[i],fontsize=10)# + str(rangez[i][0]) + r'$\to$' +str(rangez[i][1]),fontsize=10)
        ax.axhline(2,color='k',ls='--')

ax.set_xticks(np.linspace(0,8,9))
#ax.set_yticks(np.linspace(-4,8,1))
ax.set_ylim([-4,8])
axL.set_ylim([0,20])
axE.set_ylim([0,1e3])

axarr1[0][0].set_ylabel('VL $\log_{10}$(copies/mL)',fontsize=10)
axarr1[1][0].set_ylabel('VL $\log_{10}$(copies/mL)',fontsize=10)
axarr1[1][3].axis('off')

axarr2[0][0].set_ylabel('latent (cells)',fontsize=10)
axarr2[1][0].set_ylabel('latent (cells)',fontsize=10)
axarr2[1][3].axis('off')

axarr3[0][0].set_ylabel('adaptive (cells/$\mu$L)',fontsize=10)
axarr3[1][0].set_ylabel('adaptive (cells/$\mu$L)',fontsize=10)
axarr3[1][3].axis('off')

for i in range(4):
    axarr1[1][i].set_xlabel('time (weeks)')
    axarr2[1][i].set_xlabel('time (weeks)')
    axarr3[1][i].set_xlabel('time (weeks)')

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig1.savefig('figures/sensitivity_local_V.pdf')
fig2.savefig('figures/sensitivity_local_L.pdf')
fig3.savefig('figures/sensitivity_local_E.pdf')




# In[ ]:

################################################################################
#global sensitivity analysis function

sampz=10**3

lhs_samples = pyDOE.lhs( int(len(variables)), samples=sampz, criterion='center')

tt=np.linspace(0,56,1e4) #8 week time series in days

var_vals=np.zeros(len(variables))
var_vals_list=np.zeros([len(variables),sampz])

outcomes=np.zeros([sampz,5])

outcome_names=['first pos time','peak time','peak VL','max $E$','max $L$']

if sampz<30:
    fig, ax = plt.subplots(1,2,figsize=(4,3),dpi=rez,sharex=True)

rangez=np.array([[0,30],[0,600],[0,300],[0,0.3],[0,0.7],[0.8,5],[0.2,20]])

ppt_list=[]
#looping over the number of LHS samples
for i in range(sampz):

    #I3_0    = 10**(np.log10(rangez[0][0]) + lhs_samples[i][0]*(np.log10(rangez[0][1])-np.log10(rangez[0][0])))

    #make the variable values using the LHS sampling
    var_vals = rangez[:,0]+lhs_samples[i,:]*(rangez[:,1]-rangez[:,0])
    
    var_vals_list[:,i]=var_vals
    
    rand_ppt=np.random.randint(0,len(goodP_list))
    
    ppt_list.append(rand_ppt)
    
    aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr=goodP_list[rand_ppt] #choose a random model
    tc,A0,phi,r,s,hill,IC50=var_vals #choose LHS range of params

    in_tc=np.where(tY>tc)[0][0] #index of first contact
    tV=tY[in_tc:] #time for virus
    sol,vll,infpoz,t_fp,t_max=run_model(
        tV,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)

    outcomes[i,:]=t_fp, t_max, np.max(vll), np.max(sol[:,5]), np.max(sol[:,3]+sol[:,4])
        
    if sampz<30:
        #make plots showing variation              
        ax[0].plot(tV/7,vll,alpha=0.3)
        ax[0].set_ylim([0,8])
        ax[0].set_ylabel('VL $\log_{10}$(copies/mL)')
        ax[0].set_xlabel('time (weeks)')

        Ab_t,E_t = Ab_model(tt,A0,phi,r,s,hill,IC50)    #refreshes every 8 weeks
        ax[1].plot(tt/7, R0*(1-E_t),alpha=0.3)
        ax[1].set_ylabel('$\mathcal{R}(Ab,t)$')
        ax[0].set_xlabel('time (weeks)')
        ax[0].set_xticks(np.linspace(0,8,9))

if sampz<30:
    plt.tight_layout()
    plt.savefig('figures/global_sensistivity_examples.pdf')


# In[ ]:


corrz=np.zeros([len(variables),len(outcomes[0,:])])
pvalz=np.zeros([len(variables),len(outcomes[0,:])])

fig,axarr = plt.subplots(len(outcomes[0,:]),len(variables),figsize=(15,10),dpi=rez)

cz=['teal','orchid','coral','olive','slategray']

for i in range(len(variables)):

    for j in range(len(outcomes[0,:])):
        #compute correlations (ranked via spearman)
        corrz[i,j]=spearmanr(var_vals_list[i,:],outcomes[:,j])[0]
        pvalz[i,j]=spearmanr(var_vals_list[i,:],outcomes[:,j])[1]

        ax=axarr[j][i]
        #ax.scatter(var_vals_list[i,:]/np.mean(rangez,1)[i],outcomes[:,j],s=10,alpha=0.5,color=cz[j])
        ax.scatter(var_vals_list[i,:],outcomes[:,j],s=10,alpha=0.5,color=cz[j])
        ax.set_xticks(np.linspace(rangez[i,0],rangez[i,1],3))
        axarr[j][0].set_ylabel(outcome_names[j],fontsize=10)
    
    axarr[len(outcome_names)-1][i].set_xlabel(variables[i],fontsize=10)

plt.tight_layout()
plt.savefig('figures/sensitivity_global_scatters.pdf')


# In[ ]:


#correlation coefficient bar plots
fig,axarr=plt.subplots(1,len(outcome_names),figsize=(8,3),sharey=True,sharex=True,dpi=rez)
bspace=np.arange(0,len(variables)*2,2)

#sort the variable names for y axis labels
#indz=cum_corrs[:,0].argsort()
#sorted_vars=[]
#for kk in range(len(variables)):
#    sorted_vars.append(variables[indz[kk]])

for j in range(len(outcome_names)):
    ax=axarr[j]
    ax.barh(bspace,corrz[:,j],1,color=cz[j])
    
    #plot significance
    for i in range(len(variables)):
        if pvalz[i,j]<0.05:
            ax.scatter(corrz[i,j]+np.sign(corrz[i,j])*0.1,bspace[i]+0.5,color=cz[j],marker='*',s=20)
        if pvalz[i,j]<0.01:
            ax.scatter(corrz[i,j]+np.sign(corrz[i,j])*0.1,bspace[i],color=cz[j],marker='*',s=20)
    
    ax.set_title(outcome_names[j],fontsize=10)
    ax.set_xlabel('spearman '+ r'$\rho$',fontsize=10)
    #ax.axvline(-1.1,color='k',lw=0.8)
        
ax.set_xlim([-0.75,0.75])
ax.set_xticks([-0.5,0,0.5])

plt.yticks(bspace+0.4,variables,fontsize=8)
fig.tight_layout()


fig.savefig('figures/sensitivity_correlation_bars.pdf')
    


# In[ ]:


#correlation coefficient bar plots

plt.figure(figsize=(4,6),dpi=rez)
bspace=np.arange(0,len(variables)*2,2)

bw=0.2 #barwidth

for j in range(len(outcome_names)):

    plt.barh(bspace+j*bw,corrz[:,j],bw,color=cz[j])
    
    #plot significance
    for i in range(len(variables)):
        if pvalz[i,j]<0.05:
            plt.scatter(corrz[i,j]+np.sign(corrz[i,j])*0.1,bspace[i]+j*bw+0.1,color=cz[j],marker='*',s=10)
        if pvalz[i,j]<0.01:
            plt.scatter(corrz[i,j]+np.sign(corrz[i,j])*0.1,bspace[i]+j*bw,color=cz[j],marker='*',s=10)
    
plt.xlabel('spearman '+ r'$\rho$',fontsize=10)
plt.xlim([-0.65,0.65])
plt.xticks([-0.5,0,0.5])

plt.yticks(bspace+0.4,variables,fontsize=10)
plt.tight_layout()


plt.savefig('figures/sensitivity_correlation_bars_ALL.pdf')
    


# In[ ]:

#plot correlation between Ab level at first contact and other outcomes

plt.figure(figsize=(10,3),dpi=rez)
for i in range(sampz):
    tc,A0,phi,r,s,hill,IC50=var_vals_list[:,i]
    
    Ab_t,E_t = Ab_model(tt,A0,phi,r,s,hill,IC50)    #refreshes every 8 weeks

    aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g=goodP_list[ppt_list[i]]
    Rt=Bt0*tau*(1-lam)*(1-E_t)/dA*p/g*aS/dS
    for j in range(len(outcome_names)):
        plt.subplot(1,len(outcome_names),j+1)
        plt.scatter(np.log10(Ab_t[0]),outcomes[i,j],s=10,alpha=0.5,color=cz[j])
        #plt.xscale('log')
        plt.xticks(range(1,6))
        plt.ylabel(outcome_names[j],fontsize=10)
        plt.xlabel('$\log_{10}\mathcal{Y}(0)$')
plt.tight_layout()


# In[ ]:

#plot correlation between R0 first contact and other outcomes

plt.figure(figsize=(10,3),dpi=rez)
for i in range(sampz):
    tc,A0,phi,r,s,hill,IC50=var_vals_list[:,i]
    
    Ab_t,E_t = Ab_model(tt,A0,phi,r,s,hill,IC50)    #refreshes every 8 weeks

    aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g=goodP_list[ppt_list[i]]
    R0=Bt0*tau*(1-lam)*(1-E_t[0])/dA*p/g*aS/dS
    for j in range(len(outcome_names)):
        plt.subplot(1,len(outcome_names),j+1)
        plt.scatter(R0,outcomes[i,j],s=10,alpha=0.5,color=cz[j])
        #plt.xscale('log')
        #plt.xticks(range(1,6))
        plt.ylabel(outcome_names[j],fontsize=10)
        plt.xlabel('$\mathcal{R}_0$')
plt.tight_layout()


# In[ ]:

#plot correlation between AUC bnAb and outcomes

plt.figure(figsize=(10,3),dpi=rez)
for i in range(sampz):
    tc,A0,phi,r,s,hill,IC50=var_vals_list[:,i]
    
    Ab_t,E_t = Ab_model(tt,A0,phi,r,s,hill,IC50)    #refreshes every 8 weeks

    for j in range(len(outcome_names)):
        plt.subplot(1,len(outcome_names),j+1)
        plt.scatter(np.log10(np.sum(Ab_t)),outcomes[i,j],s=10,alpha=0.5,color=cz[j])
        #plt.xscale('log')
        #plt.xticks(range(1,6))
        plt.ylabel(outcome_names[j],fontsize=10)
        plt.xlabel('$\log_{10}\sum_t \mathcal{Y}(t)$')
plt.tight_layout()


# coding: utf-8

# In[1]:

#!/usr/bin/env python
#get_ipython().magic('matplotlib inline')

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import pandas as pd
import pyDOE
from scipy.stats import pearsonr,spearmanr
from scipy.interpolate import interp1d

rez=600


# In[2]:

#whole AB model
def Ab_model(t,A0,phi,r,s,hill,IC50):
    Ab_t = A0*(np.exp(-r*t)+phi*np.exp(-s*t))
    if A0>0:
        E_t = 1/(1+(Ab_t/IC50)**-hill)
    else:
        E_t=0
    return Ab_t,E_t

def model(X,t,aS,dS,Bt0,tau,lam,dI,thL,xi,k,aE,dE,E50,w,p,g,A0,phi,r,s,hill,IC50):
    
    dY = np.zeros(7);
    Ab_t,E_t = Ab_model(t,A0,phi,r,s,hill,IC50)    

    S=X[0]; AU=X[1]; AP=X[2]; LU=X[3]; LP=X[4]; E=X[5]; V=X[6];    
    
    Bt=Bt0*(1-E_t)
    
    dY[0] = aS - dS*S - Bt*S*V                         #susceptible cells
    dY[1] = (1-tau)*(1-lam)*Bt*S*V - dI*AU - k*E*AU + xi*LU    #active unproductively infected
    dY[2] = tau*(1-lam)*Bt*S*V - dI*AP - k*E*AP + xi*LP        #active productively infected
    dY[3] = (1-tau)*lam*Bt*S*V + thL*LU                #latent unproductively infected
    dY[4] = tau*lam*Bt*S*V + thL*LP                    #latent productively infected
    dY[5] = w*E*(AP+AU)/(E+E50) + aE - dE*E;           #adaptive immune system
    dY[6] = p*AP - g*V - Bt*S*V                        #virus
    return dY

def run_model(tt,X0,aS,dS,Bt,tau,lam,dI,thL,xi,k,aE,dE,E50,w,p,g,A0,phi,r,s,hill,IC50):
    
    infpoz,maxpoz,t_fp,t_max=0,0,0,0
    
    sol=spi.odeint(model, X0, tt, 
                   (aS,dS,Bt,tau,lam,dI,thL,xi,k,aE,dE,E50,w,p,g,A0,phi,r,s,hill,IC50),
                   mxstep=10000)    

    vll=np.log10(sol[:,6]/vol*1e3) #viral load as usual
    
    
    if (vll>-10).all() & (vll>2).any():        
        infpoz=np.where(vll>2)[0][0] #index of first positive
        maxpoz=np.where(vll==max(vll))[0][0]
        t_fp=tt[infpoz]
        t_max=tt[maxpoz]
    
    return sol,vll,infpoz,t_fp,t_max


# In[3]:

#fit Robb VL
robb=pd.read_csv('data/robb_ranges.csv',names=['logVL','days']); robbHI=robb[:15]; robbLO=robb[15:]

interpLO=interp1d(robbLO.days,robbLO.logVL); interpHI=interp1d(robbHI.days,robbHI.logVL)


# In[4]:

# parameters for viral dynamics model
vol = 1      # volume of blood [uL]
aS  = 70*vol;   #constant growth rate of susceptibles [cells/uL/day]
dS  = 0.3;   #susceptible death rate [1/day] 
Bt0 = 1e-4/vol  # infection rate of T-cells [uL/cells-day]/[uL]
dA  = 1.0       # active death rate [1/day]
p   = 5e4       # burst rate of virus from cells [virions/cell]
g   = 23        # virus clearance rate [1/day]
tau = 0.05      # productive infection probability []
lam = 1e-4      # latency probability []
thL = 5.2e-4    # latent clearance rate [1/day]
aL  = 0.015;    # latent proliferation rate [1/day] (Tcm)
xi  = 1e-5;     # latent activation rate [1/day]
dL  = aL-thL-xi # latent death rate
k   = 0.3/vol;  #immune cell killing rate [uL/cell-day]/[uL]
w   = 1.6;     #immune cell multiplier [1/day]
aE  = 1e-4*vol;   #initial E cell concentration [cells/uL]*[uL]
dE  = 0.002;  #immune death rate [1/day]
E50 = 250*vol;   #50 pct max E cell concentration [cells/uL]*[uL]

R0=aS*Bt0*tau*(1-lam)*p/g/dS/dA; print(R0) # basic reproductive number

#PK parameters
A0,phi,r,s = 1, 0, 0.348, 0.0538 #low 30 or high 10 dose range

#PD parameters
IC50=20 # IC50 [ug/uL]
hill=1  # hill coefficient []

V0=1e-3 #1 virion in body!

X0=np.array([aS/dS,0,0,0,0,aE/dE,V0])

tY=np.linspace(0,7*8,1e3) #time for antibody dose
Ab_t,E_t  = Ab_model(tY,A0,phi,r,s,hill,IC50) #Ab model

tc=7 #first contact at 1 week

in_tc=np.where(tY>tc)[0][0] #index of first contact

tV=tY[in_tc:] #time for virus

sol,vll,infpoz,t_fp,t_max=run_model(
    tV,X0,aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g,A0,phi,r,s,hill,IC50) #first positive at 1 week

#calculate R0(t)
#dIdt=np.diff(sol[:,2])
#I=sol[1:,2]
#dt=tt[1]
#Rt=1+dIdt/(I*dt*dA)
Rt=Bt0*tau*(1-lam)*(1-E_t)/dA*p/g*aS/dS

fig,axarr=plt.subplots(3,2,figsize=(6,5),dpi=rez,sharex=True)
titz=['S','$A_U$','$A_P$','$L_U$','$L_P$','$E$','$V$']

axarr[0][0].plot(tV/7,sol[:,0],color='gray')
axarr[0][0].set_ylabel('cells/$\mu$L')
axarr[0][0].legend(titz[0],fontsize=8,loc=2)

axarr[1][0].plot(tV/7,sol[:,2],color='red',ls='--')
axarr[1][0].plot(tV/7,sol[:,1],color='red')
axarr[1][0].set_ylabel('cells/$\mu$L')
axarr[1][0].legend(titz[1:3],fontsize=8,loc=2)

axarr[2][0].plot(tV/7,sol[:,3],color='blue',ls='--')
axarr[2][0].plot(tV/7,sol[:,4],color='blue')
axarr[2][0].legend(titz[3:],fontsize=8,loc=2)
axarr[2][0].set_ylabel('cells/$\mu$L')

axarr[0][0].set_xticks(range(9))
axarr[2][0].set_xlabel('time (weeks)')

#plot viral stuff
axarr[0][1].plot(tY/7,Ab_t,color='gold')
axarr[0][1].set_ylabel('bnAbs ($\mu$g/$\mu$L)')
    
axarr[1][1].plot(tY,Rt,color='black') #allow burn in for some R0
axarr[1][1].axhline(1,color='k',ls='--')
axarr[1][1].set_xticks(range(9))
axarr[1][1].set_ylabel(r'$\mathcal{R}(t)$')
#axarr[1].set_yticks(np.linspace(0,9,4))

axarr[2][1].plot(tV/7,vll,color='purple') #virus/mL
#axarr[2][1].plot((t_fp_shift-tc/7,vll_fp_shift,color='purple') #virus/mL
#axarr[2][1].fill_between((tY+t_fp)/7, interpLO(tY), interpHI(tY),color='gray',alpha=0.3)
#axarr[2][1].fill_between(tY/7, interpLO(tY), interpHI(tY),color='gray',alpha=0.3)
#axarr[2][1].set_ylim([np.log10(30),9])
axarr[2][1].set_xlim([0,8])
#axarr[2][1].set_yticks(np.linspace(2,8,7))
axarr[2][1].set_ylabel('viral load \n $\log_{10}$(copies/mL)')
axarr[2][1].set_xlabel('time (weeks)')

plt.tight_layout()
plt.savefig('figures/det_model.pdf')


# In[5]:

#now see which parameter ranges fit in here

sampz=10**4

tt=np.linspace(0,100,1e4)

guessP=np.array([aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g])

goodP_list=[]

bdz=1.#extra bounds for fitting shaded area

fposT=7

for i in range(sampz):
    aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr=np.random.lognormal(0,2,[15])*guessP
    
    if aSr*Bt0r*taur*(1-lamr)*pr/gr/dSr/dAr<10: # basic reproductive number

        sol,vll,infpoz,t_fp,t_max=run_model(
            tt,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)
        
        if (vll>-10).all() & (vll>2).any():        

            t_shifted = tt[infpoz:]-tt[infpoz]
            
            loV,hiV = interpLO(t_shifted),interpHI(t_shifted) #the bounds

            #first positive must be before 3 weeks and t_max must be after 3 days since first positive
            if (t_fp<21.) & (t_max-t_fp>3.):

                if ~ ((vll[infpoz:]>hiV*bdz).any() | (vll[infpoz:]<loV/bdz).any()):
                    goodP_list.append([aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr])#,x0r])
                    
                    #plt.plot(tt[infpoz:],vll[infpoz:])


# In[6]:

plt.figure(figsize=(4,3),dpi=rez)

loV,hiV = interpLO(tt),interpHI(tt) #the bounds
for i in range(len(goodP_list)):
    aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr=goodP_list[i]

    sol,vll,infpoz,t_fp,t_max=run_model(
        tt,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)

    plt.plot(tt[infpoz:]-tt[infpoz],vll[infpoz:])#,color='k',alpha=0.5)

plt.fill_between(tt, loV, hiV, color='gray',alpha=0.3)
plt.xlabel('time (days after first positive test)')
plt.ylabel('viral load $\log_{10}$(copies/mL)')
plt.tight_layout()
plt.savefig('figures/robb_examples.pdf')


# In[7]:

#boxplots from robb

#days to nadir
[17.878483268976264, 27.917385628280154, 30.934789888203735, 40.682980986119745, 41.7276975706514]

#days to EIA
[7.950909781228649,12.88320520610375,14.101772311072896, 18.04749340369393, 32.84449492559828]

#days to peak
[6.147116993275322, 6.320967513712819, 9.918987565458853, 13.052445835379505, 14.096931925353058, 17.868744873897068]

#peak log10 VL
[4.565059970202178, 6.35257087867727, 6.849094572557447, 7.663367820014483, 8.59685797501567]

#nadir log10 VL
[1.780762540024421, 3.6079697335034737, 4.422294201973418, 4.9188435063600515, 6.368654276730486]



# In[8]:

pnames=[r'$\alpha_S$',r'$\delta_S$',r'$\beta_0$',r'$\tau$',r'$\lambda$',
        r'$\delta_I$',r'$\theta_L$',r'$\xi$',r'$\kappa$',r'$\alpha_E$',
        r'$\delta_E$',r'$E_{50}$',r'$\omega$',r'$\pi$',r'$\gamma$']
        
fig,axarr=plt.subplots(3,5,sharey=False,figsize=(7,6),dpi=rez)
for j in range(15):
    ax=axarr[j%3][int(j/3.)]
    
    logpz=np.log10(np.array(goodP_list))
    
    #ax.violinplot(np.array(goodP_list)[:,j])#/guessP[j])
    ax.scatter(np.random.normal(1,0.05,len(goodP_list)),logpz[:,j],color='gray',alpha=0.5)#/guessP[j])
    ax.boxplot(logpz[:,j])#/guessP[j])
    ax.set_xlabel(pnames[j])
    ax.set_xticks([])
    ax.set_xlim([0.8,1.2])
    #ax.set_yscale('log')
#axarr[1][0].set_ylabel('deviation from \n median value')
plt.tight_layout()
plt.savefig('figures/robb_fitparams.pdf')


# In[9]:

pnames_unix=['aS','dS','Bt0','tau','lam','dA','thL','xi','k','aE','dE','E50','w','p','g']

gf=pd.DataFrame(goodP_list,columns=pnames_unix)

plt.figure(figsize=(7,3),dpi=rez)
plt.subplot(141)
plt.boxplot(gf.aS/gf.dS)
plt.yscale('log')
plt.xticks([1],[r'$\alpha_S/\delta_S$'])

plt.subplot(142)
plt.boxplot(gf.aE/gf.dE)
plt.yscale('log')
plt.xticks([1],[r'$\alpha_E/\delta_E$'])

plt.subplot(143)
plt.boxplot(gf.p/gf.g)
plt.yscale('log')
plt.xticks([1],[r'$\pi/\gamma$'])

plt.subplot(144)
plt.boxplot(gf.tau*(1-gf.lam)*gf.aS/gf.dS*gf.p/gf.g*gf.Bt0/gf.dA)
#ax.scatter(np.random.normal(1,0.1,len(gf.tau)),gf.tau*(1-gf.lam)*gf.aS/gf.dS*gf.p/gf.g*gf.Bt0/gf.dA)#/guessP[j])
#plt.yscale('log')
plt.xticks([1],[r'$\mathcal{R}_0$'])

plt.tight_layout()
plt.savefig('figures/robb_fitparams_identifiability.pdf')


# In[10]:

fig,axarr=plt.subplots(15,15,sharey=True,sharex=True,figsize=(10,10),dpi=rez)

cz_range=np.linspace(0,1,11)

for i in range(15):
    for j in range(15):
        ax=axarr[i][j]

        #only plot bottom diagonal
        if j>i-1:
            ax.axis('off')
        else:
            pcorr=spearmanr(np.log10(gf.iloc[:,i]/guessP[i]),np.log10(gf.iloc[:,j]/guessP[j]))
            #if pcorr[1]<0.05:
            #    cz='crimson'
            #else:
            #    cz='green'
            ax.scatter(np.log10(gf.iloc[:,i]/guessP[i]),np.log10(gf.iloc[:,j]/guessP[j]),s=10, color=[0,0,pcorr[1]],alpha=0.3)

    axarr[i][0].set_ylabel(pnames[i], fontsize=12)
    axarr[14][i].set_xlabel(pnames[i], fontsize=12)

ax.set_ylim([-2,2])
ax.set_xlim([-2,2])
ax.set_xticks([-1,0,1])
ax.set_yticks([-1,0,1])
ax.set_xticklabels([-1,0,1],fontsize=8)
ax.set_yticklabels([-1,0,1],fontsize=8)
fig.subplots_adjust(wspace=0.01, hspace=0.01)
fig.savefig('figures/robb_fitparams_varcorrs.pdf')


# In[11]:

################################################################################
#local sensitivity analysis function

#variables to check sensitivity
variables = ['first contact time ' + r'$t_c$',
             'bnAb conc phase 1 ' + r'$\mathcal{Y}_1$', 
             'bnAb conc phase 2 ' + r'$\mathcal{Y}_2$',
             'decay rate 1 ' + r'$k_1$',
             'decay rate 2 ' + r'$k_2$',
             'Hill coefficient ' + r'$h$',
             '50% effectiveness ' + r'$IC_{50}$']

#from data_analysis fitting:
#[86  2.5 0.15 0.007] +/- [ 32.7 32.6 0.12 0.29 ]
#[240 2.9 0.086 -0.011] +/- [ 112 116 0.081 0.71]

#ranges for best (b) and worst (w) cases
b_case=np.array([0 ,600,300,0  ,0  ,5  ,0.2])
w_case=np.array([30,0  ,0  ,0.3,0.7,0.8,20 ])
                            
tY=np.linspace(0,7*8,1e4) #time for antibody dose 8 week time series in days

fig1,axarr1 = plt.subplots(2,4,sharex=True,sharey=True,figsize=(9,4),dpi=rez)
fig2,axarr2 = plt.subplots(2,4,sharex=True,sharey=True,figsize=(9,4),dpi=rez)
fig3,axarr3 = plt.subplots(2,4,sharex=True,sharey=True,figsize=(9,4),dpi=rez)
#looping over the number of variables to check sensitivity

num_ppts=25
for k in range(num_ppts): 
    aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr=goodP_list[k] #choose a random model
    for i in range(len(variables)):
        ax=axarr1[i%2][int(i/2.)]
        axL=axarr2[i%2][int(i/2.)]
        axE=axarr3[i%2][int(i/2.)]

        lo=np.array([7,200,50,0.15,0.007,1,1]); hi=np.array([7,200,50,0.15,0.007,1,1]) #set both at typical values

        #change the i-th parameter
        lo[i]=w_case[i]; hi[i]=b_case[i]

        #evaluate model twice
        tc,A0,phi,r,s,hill,IC50=lo #choose low range of params
        in_tc=np.where(tY>tc)[0][0] #index of first contact
        tV=tY[in_tc:] #time for virus
        sol,vll,infpoz,t_fp,t_max=run_model(
            tV,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)

        loV=vll #for logplot later
        ax.plot(tV/7,loV,alpha=0.7,color='darkred')
        
        axL.plot(tV/7,sol[:,3],color='darkred',ls='--')
        axL.plot(tV/7,sol[:,4],color='darkred')
        
        axE.plot(tV/7,sol[:,5],color='darkred')

        tc,A0,phi,r,s,hill,IC50=hi #choose high range of params
        in_tc=np.where(tY>tc)[0][0] #index of first contact
        tV=tY[in_tc:] #time for virus
        sol,vll,infpoz,t_fp,t_max=run_model(
            tV,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)

        hiV=vll
        ax.plot(tV/7,hiV,alpha=0.7,color='deepskyblue')
        axL.plot(tV/7,sol[:,3]*1e6/vol,color='deepskyblue',ls='--')
        axL.plot(tV/7,sol[:,4]*1e6/vol,color='deepskyblue')
        
        axE.plot(tV/7,sol[:,5],color='deepskyblue')

        #ax.fill_between(lot/7,lo_sol[:,6]/vol*1e3,hi_sol[:,6]/vol*1e3,alpha=0.6)

        #ax.set_ylim([0.1,100])
        ax.set_title(variables[i],fontsize=10)# + str(rangez[i][0]) + r'$\to$' +str(rangez[i][1]),fontsize=10)
        axL.set_title(variables[i],fontsize=10)# + str(rangez[i][0]) + r'$\to$' +str(rangez[i][1]),fontsize=10)
        axE.set_title(variables[i],fontsize=10)# + str(rangez[i][0]) + r'$\to$' +str(rangez[i][1]),fontsize=10)
        ax.axhline(2,color='k',ls='--')

ax.set_xticks(np.linspace(0,8,9))
#ax.set_yticks(np.linspace(-4,8,1))
ax.set_ylim([-4,8])
axL.set_ylim([0,20])
axE.set_ylim([0,1e3])

axarr1[0][0].set_ylabel('VL $\log_{10}$(copies/mL)',fontsize=10)
axarr1[1][0].set_ylabel('VL $\log_{10}$(copies/mL)',fontsize=10)
axarr1[1][3].axis('off')

axarr2[0][0].set_ylabel('latent (cells)',fontsize=10)
axarr2[1][0].set_ylabel('latent (cells)',fontsize=10)
axarr2[1][3].axis('off')

axarr3[0][0].set_ylabel('adaptive (cells/$\mu$L)',fontsize=10)
axarr3[1][0].set_ylabel('adaptive (cells/$\mu$L)',fontsize=10)
axarr3[1][3].axis('off')

for i in range(4):
    axarr1[1][i].set_xlabel('time (weeks)')
    axarr2[1][i].set_xlabel('time (weeks)')
    axarr3[1][i].set_xlabel('time (weeks)')

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig1.savefig('figures/sensitivity_local_V.pdf')
fig2.savefig('figures/sensitivity_local_L.pdf')
fig3.savefig('figures/sensitivity_local_E.pdf')




# In[ ]:

################################################################################
#global sensitivity analysis function

sampz=10**3

lhs_samples = pyDOE.lhs( int(len(variables)), samples=sampz, criterion='center')

tt=np.linspace(0,56,1e4) #8 week time series in days

var_vals=np.zeros(len(variables))
var_vals_list=np.zeros([len(variables),sampz])

outcomes=np.zeros([sampz,5])

outcome_names=['first pos time','peak time','peak VL','max $E$','max $L$']

if sampz<30:
    fig, ax = plt.subplots(1,2,figsize=(4,3),dpi=rez,sharex=True)

rangez=np.array([[0,30],[0,600],[0,300],[0,0.3],[0,0.7],[0.5,1.5],[0.2,20]])

ppt_list=[]
#looping over the number of LHS samples
for i in range(sampz):

    #I3_0    = 10**(np.log10(rangez[0][0]) + lhs_samples[i][0]*(np.log10(rangez[0][1])-np.log10(rangez[0][0])))

    #make the variable values using the LHS sampling
    var_vals = rangez[:,0]+lhs_samples[i,:]*(rangez[:,1]-rangez[:,0])
    
    var_vals_list[:,i]=var_vals
    
    rand_ppt=np.random.randint(0,len(goodP_list))
    
    ppt_list.append(rand_ppt)
    
    aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr=goodP_list[rand_ppt] #choose a random model
    tc,A0,phi,r,s,hill,IC50=var_vals #choose LHS range of params

    in_tc=np.where(tY>tc)[0][0] #index of first contact
    tV=tY[in_tc:] #time for virus
    sol,vll,infpoz,t_fp,t_max=run_model(
        tV,X0,aSr,dSr,Bt0r,taur,lamr,dAr,thLr,xir,kr,aEr,dEr,E50r,wr,pr,gr,A0,phi,r,s,hill,IC50)

    outcomes[i,:]=t_fp, t_max, np.max(vll), np.max(sol[:,5]), np.max(sol[:,3]+sol[:,4])
        
    if sampz<30:
        #make plots showing variation              
        ax[0].plot(tV/7,vll,alpha=0.3)
        ax[0].set_ylim([0,8])
        ax[0].set_ylabel('VL $\log_{10}$(copies/mL)')
        ax[0].set_xlabel('time (weeks)')

        Ab_t,E_t = Ab_model(tt,A0,phi,r,s,hill,IC50)    #refreshes every 8 weeks
        ax[1].plot(tt/7, R0*(1-E_t),alpha=0.3)
        ax[1].set_ylabel('$\mathcal{R}(Ab,t)$')
        ax[0].set_xlabel('time (weeks)')
        ax[0].set_xticks(np.linspace(0,8,9))

if sampz<30:
    plt.tight_layout()
    plt.savefig('figures/global_sensistivity_examples.pdf')


# In[ ]:


corrz=np.zeros([len(variables),len(outcomes[0,:])])
pvalz=np.zeros([len(variables),len(outcomes[0,:])])

fig,axarr = plt.subplots(len(outcomes[0,:]),len(variables),figsize=(15,10),dpi=rez)

cz=['teal','orchid','coral','olive','slategray']

for i in range(len(variables)):

    for j in range(len(outcomes[0,:])):
        #compute correlations (ranked via spearman)
        corrz[i,j]=spearmanr(var_vals_list[i,:],outcomes[:,j])[0]
        pvalz[i,j]=spearmanr(var_vals_list[i,:],outcomes[:,j])[1]

        ax=axarr[j][i]
        #ax.scatter(var_vals_list[i,:]/np.mean(rangez,1)[i],outcomes[:,j],s=10,alpha=0.5,color=cz[j])
        ax.scatter(var_vals_list[i,:],outcomes[:,j],s=10,alpha=0.5,color=cz[j])
        ax.set_xticks(np.linspace(rangez[i,0],rangez[i,1],3))
        axarr[j][0].set_ylabel(outcome_names[j],fontsize=10)
    
    axarr[len(outcome_names)-1][i].set_xlabel(variables[i],fontsize=10)

plt.tight_layout()
plt.savefig('figures/sensitivity_global_scatters.pdf')


# In[ ]:


#correlation coefficient bar plots
fig,axarr=plt.subplots(1,len(outcome_names),figsize=(8,3),sharey=True,sharex=True,dpi=rez)
bspace=np.arange(0,len(variables)*2,2)

#sort the variable names for y axis labels
#indz=cum_corrs[:,0].argsort()
#sorted_vars=[]
#for kk in range(len(variables)):
#    sorted_vars.append(variables[indz[kk]])

for j in range(len(outcome_names)):
    ax=axarr[j]
    ax.barh(bspace,corrz[:,j],1,color=cz[j])
    
    #plot significance
    for i in range(len(variables)):
        if pvalz[i,j]<0.05:
            ax.scatter(corrz[i,j]+np.sign(corrz[i,j])*0.1,bspace[i]+0.5,color=cz[j],marker='*',s=20)
        if pvalz[i,j]<0.01:
            ax.scatter(corrz[i,j]+np.sign(corrz[i,j])*0.1,bspace[i],color=cz[j],marker='*',s=20)
    
    ax.set_title(outcome_names[j],fontsize=10)
    ax.set_xlabel('spearman '+ r'$\rho$',fontsize=10)
    #ax.axvline(-1.1,color='k',lw=0.8)
        
ax.set_xlim([-0.75,0.75])
ax.set_xticks([-0.5,0,0.5])

plt.yticks(bspace+0.4,variables,fontsize=8)
fig.tight_layout()


fig.savefig('figures/sensitivity_correlation_bars.pdf')
    


# In[ ]:


#correlation coefficient bar plots

plt.figure(figsize=(4,6),dpi=rez)
bspace=np.arange(0,len(variables)*2,2)

bw=0.2 #barwidth

for j in range(len(outcome_names)):

    plt.barh(bspace+j*bw,corrz[:,j],bw,color=cz[j])
    
    #plot significance
    for i in range(len(variables)):
        if pvalz[i,j]<0.05:
            plt.scatter(corrz[i,j]+np.sign(corrz[i,j])*0.1,bspace[i]+j*bw+0.1,color=cz[j],marker='*',s=10)
        if pvalz[i,j]<0.01:
            plt.scatter(corrz[i,j]+np.sign(corrz[i,j])*0.1,bspace[i]+j*bw,color=cz[j],marker='*',s=10)
    
plt.xlabel('spearman '+ r'$\rho$',fontsize=10)
plt.xlim([-0.65,0.65])
plt.xticks([-0.5,0,0.5])

plt.yticks(bspace+0.4,variables,fontsize=10)
plt.tight_layout()

plt.savefig('figures/sensitivity_correlation_bars_ALL.pdf')
    


# In[ ]:

#plot correlation between Ab level at first contact and other outcomes

plt.figure(figsize=(10,3),dpi=rez)
for i in range(sampz):
    tc,A0,phi,r,s,hill,IC50=var_vals_list[:,i]
    
    Ab_t,E_t = Ab_model(tt,A0,phi,r,s,hill,IC50)    #refreshes every 8 weeks

    aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g=goodP_list[ppt_list[i]]
    Rt=Bt0*tau*(1-lam)*(1-E_t)/dA*p/g*aS/dS
    for j in range(len(outcome_names)):
        plt.subplot(1,len(outcome_names),j+1)
        plt.scatter(np.log10(Ab_t[0]),outcomes[i,j],s=10,alpha=0.5,color=cz[j])
        #plt.xscale('log')
        plt.xticks(range(1,6))
        plt.ylabel(outcome_names[j],fontsize=10)
        plt.xlabel('$\log_{10}\mathcal{Y}(0)$')
plt.tight_layout()
plt.savefig('figures/corr_AbTC.pdf')


# In[ ]:

#plot correlation between R0 first contact and other outcomes

plt.figure(figsize=(10,3),dpi=rez)
for i in range(sampz):
    tc,A0,phi,r,s,hill,IC50=var_vals_list[:,i]
    
    Ab_t,E_t = Ab_model(tt,A0,phi,r,s,hill,IC50)    #refreshes every 8 weeks

    aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g=goodP_list[ppt_list[i]]
    R0=Bt0*tau*(1-lam)*(1-E_t[0])/dA*p/g*aS/dS
    for j in range(len(outcome_names)):
        plt.subplot(1,len(outcome_names),j+1)
        plt.scatter(R0,outcomes[i,j],s=10,alpha=0.5,color=cz[j])
        #plt.xscale('log')
        #plt.xticks(range(1,6))
        plt.ylabel(outcome_names[j],fontsize=10)
        plt.xlabel('$\mathcal{R}_0$')
plt.tight_layout()
plt.savefig('figures/corr_R0.pdf')


# In[ ]:

#plot correlation between AUC bnAb and outcomes

plt.figure(figsize=(10,3),dpi=rez)
for i in range(sampz):
    tc,A0,phi,r,s,hill,IC50=var_vals_list[:,i]
    
    Ab_t,E_t = Ab_model(tt,A0,phi,r,s,hill,IC50)    #refreshes every 8 weeks

    for j in range(len(outcome_names)):
        plt.subplot(1,len(outcome_names),j+1)
        plt.scatter(np.log10(np.sum(Ab_t)),outcomes[i,j],s=10,alpha=0.5,color=cz[j])
        #plt.xscale('log')
        #plt.xticks(range(1,6))
        plt.ylabel(outcome_names[j],fontsize=10)
        plt.xlabel('$\log_{10}\sum_t \mathcal{Y}(t)$')
plt.tight_layout()
plt.savefig('figures/corr_AbAUC.pdf')

