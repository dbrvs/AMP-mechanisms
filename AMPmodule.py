#!/usr/bin/env python
#holds the classes for AMP sims

### TO DO
#add in 2,4,6,8,12,24 for timing
#diversity model!
#bias and precision of timing
#incomplete neutralization instead of max IC50?
#remove aL and ksi

import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv

#class to keep information about a single participant
class participant:

    #initialize participant class
    def __init__(self, name, et, tF, stop2, A0, vdp, pkp, pdp, v_flg):
        
        self.v_flg = v_flg #Boolean flag, if True prints stuff out during sim

        self.dose_interval = 56 #how often bnabs are dosed (days)
        self.obs_interval = 28 #how often measurements are observed (days)
        
        self.tF=tF #maximum time of sim, gets changed later...
        self.obs=np.zeros([int(np.ceil(self.tF/self.obs_interval)+2)]) #list for observed viral loads 21 visits including day 0, extra for positive on 21st visit
        self.obs_times = np.arange(0,len(self.obs))*self.obs_interval

        self.name = name #ppt name, usually a number when looping in a trial
        self.ex_t = et #exposure time            
                
        self.stop2 = stop2 #make the sim stop after 2 positives if True
        self.A0 = A0 #initial number of infected cells

        self.ct=0; self.bt=0 #initialize clearance and breakthru times
        self.fp_t=0 #first positive time (days)
        self.fp_v=0 #first positive viral load
        
        self.rl=np.zeros(20); self.T=np.zeros([7,20]);  #initialize state/rate matrices
        self.dt=0.01 #timestep for tauleap
        
        self.set_parameters(vdp,pkp,pdp) #set to standard values first
            
        self.pk0,self.In0=self.ab_sim(self.ex_t); #calculate concentration at transmission

        self.rs0=self.pk0/self.IC50
        
    #enables user supplied parameter set
    def set_parameters(self, vdp, pkp, pdp):
    
        #start at control case, no bnabs at all
        #self.Y1,self.Y2,self.k1,self.k2,self.h,self.IC50 = 1e-10,0,0,0,1,1e10 #note h=1 otherwise weird 0^0 problem
        
        #set all values at a good model parameter set
        vol=1e8
        self.vol = vol        #volume of blood [uL]
        self.aS  = 70*vol;    #constant growth rate of susceptibles [cells/uL/day]
        self.Bt0 = 1e-4/vol   #infection rate of T-cells [uL/cells-day]/[uL]
        self.dA  = 0.8        #active death rate [1/day]
        self.p   = 5e4        #burst rate of virus from cells [virions/cell]
        self.g   = 23         #virus clearance rate [1/day]
        self.lam = 1e-4       #latency probability []
        self.thL = 5.2e-4     #latent clearance rate [1/day]
        self.aE  = 1e-5*vol;  #initial E cell concentration [cells/uL]*[uL]
        self.aL  = 0.015;     #latent proliferation rate [1/day] (Tcm)
        self.ksi = 1e-5;      #latent activation rate [1/day]
        self.dL  = self.aL-self.thL-self.ksi #latent death rate

        #initial values, will be reset by supplied list
        self.dS  = 0.3;       #susceptible death rate [1/day] 
        self.tau = 0.06       #productive infection probability []
        self.kap = 1/vol;     #immune cell killing rate [uL/cell-day]/[uL]
        self.dE  = 0.002;     #immune death rate [1/day]
        self.E50 = 250*vol;   #50 pct max E cell concentration [cells/uL]*[uL]
        self.w   = 1.6;       #immune cell multiplier [1/day]

        R0,self.dS,self.tau,kap,self.dE,E50,self.w=vdp 
        self.kap=kap/self.vol
        self.E50=E50*self.vol
        
        self.Y1,self.k1,self.Y2,self.k2=pkp
        self.IC50,self.h=pdp #note swap

    #calculate Ab model at any time (note dosing interval mod)
    def ab_sim(self,t):
        Y_t = self.Y1*np.exp(-self.k1*(t%self.dose_interval))+self.Y2*np.exp(-self.k2*(t%self.dose_interval)) #VRC01 concentration
        nu_t = 1/(1+(self.IC50/Y_t)**self.h) #neutralization
        return Y_t,nu_t
    
    #function that updates stochastic rates/transitions
    def update_rates(self, ti, xi):

        S,AU,AP,LU,LP,E,V=xi; #i-th state, easier to read if variables actually spelled out
    
        [vdp, pkp, pdp]=self.return_param_list()
        aS,dS,Bt0,tau,lam,dA,thL,aL,dL,ksi,kap,aE,dE,E50,w,p,g=vdp
        Y1,k1,Y2,k2=pkp
        IC50,h=pdp 

        Y_t,nu_t=self.ab_sim(ti); Bt=Bt0*(1-nu_t) #recalculate infectivity based on bnabs
        
        pr=np.random.poisson(p) #burst size randomized a bit

        self.rl[0] = aS;                     self.T[:,0] =[1,0,0,0,0,0,0];  #constant production 
        self.rl[1] = dS*S;                   self.T[:,1] =[-1,0,0,0,0,0,0]  #density dependent susceptible death
        self.rl[2] = (1-tau)*(1-lam)*Bt*S*V; self.T[:,2] =[-1,1,0,0,0,0,-1] #unproductive active infection
        self.rl[3] = tau*(1-lam)*Bt*S*V;     self.T[:,3] =[-1,0,1,0,0,0,-1] #productive active infection
        self.rl[4] = (1-tau)*lam*Bt*S*V;     self.T[:,4] =[-1,0,0,1,0,0,-1] #unproductive latent infection
        self.rl[5] = tau*lam*Bt*S*V;         self.T[:,5] =[-1,0,0,0,1,0,-1] #productive latent infection
        self.rl[6] = dA*AU;                  self.T[:,6] =[0,-1,0,0,0,0,0]  #unproductive active death
        self.rl[7] = dA*AP;                  self.T[:,7] =[0,0,-1,0,0,0,pr]  #productive active death and virus burst
        self.rl[8] = dL*LU;                  self.T[:,8] =[0,0,0,-1,0,0,0]  #unproductive latent death
        self.rl[9] = dL*LP;                  self.T[:,9] =[0,0,0,0,-1,0,0]  #productive latent death
        self.rl[10] = aL*LU;                 self.T[:,10]=[0,0,0,1,0,0,0]  #unproductive latent proliferation
        self.rl[11] = aL*LP;                 self.T[:,11]=[0,0,0,0,1,0,0]  #productive latent proliferation
        self.rl[12] = ksi*LU;                self.T[:,12]=[0,1,0,-1,0,0,0]  #unproductive latent activation
        self.rl[13] = ksi*LP;                self.T[:,13]=[0,0,1,0,-1,0,0]  #productive latent activation
        self.rl[14] = kap*E*AU;              self.T[:,14]=[0,-1,0,0,0,0,0]  #unproductive active immune removal
        self.rl[15] = kap*E*AP;              self.T[:,15]=[0,0,-1,0,0,0,0]  #productive active immune removal
        self.rl[16] = w*E*(AP+AU)/(E+E50);   self.T[:,16]=[0,0,0,0,0,1,0]  #immune cell recruitment
        self.rl[17] = aE;                    self.T[:,17]=[0,0,0,0,0,1,0]  #immune cell birth
        self.rl[18] = dE*E;                  self.T[:,18]=[0,0,0,0,0,-1,0]  #immune cell clearance
        self.rl[19] = g*V;                   self.T[:,19]=[0,0,0,0,0,0,-1]  #innate viral clearance

    #function that solves stochastically using tau-leap method
    def vd_sim(self):
                
        if self.v_flg:
            print('ppt:',self.name,', exposed:',self.ex_t,', transmission ratio (conc/IC50):',self.pk0/self.IC50) #print info!
        
        tl=[0, self.ex_t]; #initial time list
        
        #initial state list
        S0=self.aS/self.dS
        E0=self.aE/self.dE
        xl=[np.array([S0,0,0,0,0,E0,0]),np.array([S0,0,self.A0,0,0,E0,0])]; #initial state matrix

        obs_ind=int(np.floor(self.ex_t/self.obs_interval)) #set index at most recent observation time pre-exposure
        fp_ind=0 #index for first observed viral load (if any)
        
        t_ind=1; ti=tl[t_ind] #counter for dt, start on second when et happens
        
        #loop over time 
        while ti<self.tF:
            xi=xl[t_ind]; ti=tl[t_ind]
            
            #breaking criteria for clearance, all sell states and virus = 0
            if (sum(xi[1:5])+xi[6])==0:
                self.ct=ti
                if self.v_flg:
                    print('clearance:',self.ct)
                break
                        
            #breaking criteria for breakthrough infection
            vl=xi[6]/self.vol*1e3; #viral load in copies/mL of blood (implies "instant" lymph trafficking to blood)
            if vl>30 and self.bt==0: 
                self.bt=ti; #breakthru time
                if self.v_flg:
                    print('breakthrough:',self.bt)
                if self.stop2:
                    self.tF=4*7*(np.ceil(self.bt/(4*7))+1)+self.dt; #change the end point (days) to be after 2 visits (as required, first pos and confirmatory)

            #do the observations!
            if ti>(obs_ind*self.obs_interval):
                self.obs[obs_ind]=vl
                obs_ind+=1
                if vl>30 and fp_ind==0:
                    self.fp_t=ti #first positive week, rounded up
                    self.fp_v=vl #first positive viral load
                    fp_ind+=1
            
            self.update_rates(ti,xi) #make new rate/transition matrices vector
            events=np.random.poisson(self.rl*self.dt) #calculate events
            xi=xi+np.sum(self.T*events,1); xi[xi<0]=0 #update state variable and make sure no negative numbers
            ti=ti+self.dt
            xl.append(xi); tl.append(ti) #the list of states
            t_ind+=1

        return np.array(tl),np.array(xl)
            
    #helper because want lists sometimes
    def return_param_list(self):
        vd_params=(self.aS,self.dS,self.Bt0,self.tau,self.lam,self.dA,self.thL,self.aL,self.dL,
                        self.ksi,self.kap,self.aE,self.dE,self.E50,self.w,self.p,self.g)
        pk_params=self.Y1,self.k1,self.Y2,self.k2
        pd_params=self.IC50,self.h
        return vd_params, pk_params, pd_params


#class for simulating trials
class trial:
    def __init__(self, name, nps, tF, dose, clade, rfrac, A0, iv_phi, maxIC50, v_flg, keep_flg):
        
        self.name=name; self.nps=nps; self.tF=tF; self.dose=dose; self.clade=clade; #trial information

        self.v_flg=v_flg #verbose flag, usually False
        self.keep_flg = keep_flg #keep all observations if 1, keep all ppts if 2
        
        self.A0=A0; self.iv_phi=iv_phi; self.maxIC50=maxIC50 #variables for sensitivity analysis

        #import parameters and update for sensitivity
<<<<<<< HEAD
        self.VD=np.array(pd.DataFrame.read_csv('data/viral_dynamics.csv'))
        if dose!=0:
            self.PK=np.array(pd.DataFrame.read_csv('data/PK'+str(self.dose)+'.csv'))
=======
        self.VD=np.array(read_csv('data/viral_dynamics.csv',usecols=range(1,8)))

        if dose!=0:
            self.PK=np.array(read_csv('data/PK'+str(self.dose)+'.csv',usecols=range(1,5)))

>>>>>>> monolix_modeling

            if clade=='bimodal':
                #pick 1000 strains, rfrac of which are resistant, i.e. ic50 is high, >50
                nn=1000
                nh=np.int(rfrac*nn) #number resistant
                
                ICH=10**np.random.normal(3,1,[nh]) #high IC50
                ICL=10**np.random.normal(-2,1,[nn-nh]) #low IC50
                
                IC2=np.append(ICL,ICH)
                PD=np.random.normal(1,0.1,[nn,2]) #make the hill slopes
                PD[PD<0]=1 #just to double check
                PD[:,0]=IC2
            else:
<<<<<<< HEAD
                PD=np.array(pd.DataFrame.read_csv('data/PD'+self.clade+'.csv'))
=======
                PD=np.array(read_csv('data/PD'+self.clade+'.csv',usecols=range(1,3)))
>>>>>>> monolix_modeling
                                
            #update PD for sensitivity
            if self.maxIC50 > 50:
                PD[PD[:,0]>49,0]=50+np.random.random([sum(PD[:,0]>49)])*self.maxIC50 #reset to max IC50
            if self.iv_phi > 1:
                PD[:,0]=PD[:,0]*self.iv_phi            
            self.PD=PD #set vals        
        
        else:
            self.PK='no VRC01!'
            self.PD='no VRC01!'
            
        #containers for output variables
        self.ex_ts=np.zeros([nps]); #time of exposure (days)
        self.ex_cs=np.zeros([nps]); #exposure concentration (\mug/mL)
        self.ex_rs=np.zeros([nps]); #exposure concentration ratio (c/IC50)
        self.fp_ts=np.zeros([nps]); #first pos detected index (0,1,2 for never interval1 and interval 2 respectively)
        self.fp_vs=np.zeros([nps]); #first pos detected viral load
        self.brks=np.zeros([nps]);  self.clrs=np.zeros([nps]); #time of brks and clrs(days)
        self.ic50s=np.zeros([nps]);  #IC50s 
        
        #self.ex_rs=np.zeros([nps]);  #ratios of true IC50 to concentration at exposure 
        
        self.obz=[[],[]] #for t and obs lists, only filled if keep_flg>0
        self.sim=[[],[]] #for t and sol lists, only filled if keep_flg>1
        self.run_trial() #run!
        
    #functino to run the trial
    def run_trial(self):
        #loop ppts with random samples from PK, PD, and viral dynamics data, random t0
        for n in range(self.nps):
            t0=np.random.rand()*self.tF #random time of crossing mucosal barrierr

            vdp=self.VD[np.random.randint(len(self.VD)),:]

            #allow for control trial
            if self.dose==0:
                pkp=[1e-5,0,0,0];
                pdp=[1e5,1];        
            else:
                pkp=self.PK[np.random.randint(len(self.PK)),:]
                pdp=self.PD[np.random.randint(len(self.PD)),:]
            
            ppt=participant(name=n, A0=self.A0, et=t0, tF=self.tF, stop2=True, vdp=vdp, pkp=pkp, pdp=pdp, v_flg=False)
            t,sol=ppt.vd_sim() #simulate ppt transmission event

            #verbose print outs
            if self.v_flg:
                print('name:',ppt.name,
                      'conc ratio:',ppt.pk0/ppt.IC50,
                      'first pos time:',ppt.fp_t,
                      'first pos VL:',ppt.fp_v)
            
            self.ex_ts[n]  = ppt.ex_t #exposure time
            self.ex_cs[n]  = ppt.pk0 #exposure time
            self.ex_rs[n]  = ppt.rs0 #exposure time
            self.fp_ts[n]  = ppt.fp_t #first positive time
            self.fp_vs[n]  = ppt.fp_v #first positive viral loads
            self.brks[n]   = ppt.bt #true breakthrough time
            self.clrs[n]   = ppt.ct #true clearance time
            self.ic50s[n]  = ppt.IC50 #all IC50s
                
            if self.keep_flg>0:
                self.obz[0].append(ppt.obs_times)
                self.obz[1].append(ppt.obs)
                if self.keep_flg>1:
                    self.sim[0].append(t)
                    self.sim[1].append(np.log10(sol[:,6]/ppt.vol*1e3+1e-3))

#function for plotting stochastic stuff
def simplot(t,sol,tAb,ppt,fn):
    
    plt.figure(figsize=(9,5))
    
    plt.subplot(231)
    plt.step(t/7,(sol[:,0]-sol[0,0])/sol[0,0]*100,lw=2,color='gray')
    plt.xlabel('time (weeks)')
    plt.ylabel('susceptible cells \n (% deviation from $S_0$)')
    plt.xticks(ppt.obs_times[::2]/7)
    plt.xlim([0,80])
    
    plt.subplot(232)
    plt.step(t/7,sol[:,1]+sol[:,2],lw=2,where='post',color='tomato') #all infected
    #plt.ylim([0,100])
    plt.xlabel('time (weeks)')
    plt.ylabel('total infected cells')
    plt.xticks(ppt.obs_times[::2]/7)
    plt.xlim([0,80])
    
    plt.subplot(233)
    plt.step(t/7,(sol[:,3]+sol[:,4])/sol[0,0]*1e6,lw=2,color='navy') #divide by S0 (CD4+ T cells) and then multiply by 1 million
    plt.step(t/7,sol[:,4]/sol[0,0]*1e6,lw=2,color='teal') #divide by S0 (CD4+ T cells) and then multiply by 1 million
    plt.legend(['total','infectious'],fontsize=10)
    plt.ylim([-1,150])
    plt.xlabel('time (weeks)')
    plt.ylabel('latent cells (per mil CD4+)')
    plt.xticks(ppt.obs_times[::2]/7)
    plt.xlim([0,80])
    
    plt.subplot(234)
    plt.step(t/7,sol[:,5]*ppt.kap/ppt.dA,color='lawngreen',lw=2)
    plt.xlabel('time (weeks)')
    plt.ylabel('adaptive killing ratio ($\kappa E/\delta_I$)')
    plt.xticks(ppt.obs_times[::2]/7)
    plt.xlim([0,80])
    
    plt.subplot(235)
    vll=np.log10(sol[:,6]/ppt.vol*1e3+1e-3)
    plt.step(t/7,vll,lw=2,color='violet',alpha=0.7)
    plt.scatter(ppt.obs_times/7,np.log10(ppt.obs+30),facecolors='none',edgecolors='k',s=40)
    plt.axhline(np.log10(30),color='k',ls='--')
    plt.ylim([0,9])
    plt.xlabel('time (weeks)')
    plt.ylabel('viral load log10(copies/mL)')
    plt.xticks(ppt.obs_times[::2]/7)
    plt.xlim([0,80])
    
    Y_t,I_t=ppt.ab_sim(tAb)
    plt.subplot(236)
    plt.plot(tAb/7,I_t*100,lw=2,color='k')
    plt.ylim([0,100])
    plt.xlabel('time (weeks)')
    plt.ylabel('VRC01 neutralization (%)')
    plt.xticks(ppt.obs_times[::2]/7)
    plt.xlim([0,80])
    
    plt.tight_layout()
    
    plt.savefig('figures/'+fn+'.pdf',dpi=600)
