#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

#the functions to model the same process stochastically with the hybrid stochastic model

#whole Ab model
def Ab_model(t,ab_ps):
    Y1,Y2,k1,k2,h,IC50=ab_ps
    Y_t = Y1*np.exp(-k1*t)+Y2*np.exp(-k2*t)
    I_t = 1/(1+(Y_t/IC50)**-h)
    return Y_t,I_t

#updating the rates of events
num_rates=21; num_states=7; rl=np.zeros(num_rates); T=np.zeros([num_states,num_rates])
def update_rates(X,ti,v_ps,ab_ps):
        
    S,AU,AP,LU,LP,E,V=X; #easier to read if variables actually spelled out

    aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g=v_ps #viral dynamics parameters
    Y1,Y2,k1,k2,h,IC50=ab_ps #bnab parameters

    Bt=Bt0*(1-1/(1+((Y1*np.exp(-k1*ti)+Y2*np.exp(-k2*ti))/IC50)**-h)) #compute Ab effect on infectivity
                
    pr=np.random.poisson(p)    #kp=0.1; pr=np.random.gamma(kp,p)

    rl[0] = aS;                     T[:,0] =[1,0,0,0,0,0,0];  #constant production 
    rl[1] = dS*S;                   T[:,1] =[-1,0,0,0,0,0,0]  #density dependent susceptible death
    rl[2] = (1-tau)*(1-lam)*Bt*S*V; T[:,2] =[-1,1,0,0,0,0,-1] #unproductive active infection
    rl[3] = tau*(1-lam)*Bt*S*V;     T[:,3] =[-1,0,1,0,0,0,-1] #productive active infection
    rl[4] = (1-tau)*lam*Bt*S*V;     T[:,4] =[-1,0,0,1,0,0,-1] #unproductive latent infection
    rl[5] = tau*lam*Bt*S*V;         T[:,5] =[-1,0,0,0,1,0,-1] #productive latent infection
    rl[6] = dA*AU;                  T[:,6] =[0,-1,0,0,0,0,0]  #unproductive active death
    rl[7] = dA*AP;                  T[:,7] =[0,0,-1,0,0,0,pr]  #productive active death and virus burst
    rl[8] = dL*LU;                  T[:,8] =[0,0,0,-1,0,0,0]  #unproductive latent death
    rl[9] = dL*LP;                  T[:,9] =[0,0,0,0,-1,0,0]  #productive latent death
    rl[10] = aL*LU;                 T[:,10]=[0,0,0,1,0,0,0]  #unproductive latent proliferation
    rl[11] = aL*LP;                 T[:,11]=[0,0,0,0,1,0,0]  #productive latent proliferation
    rl[12] = xi*LU;                 T[:,12]=[0,1,0,-1,0,0,0]  #unproductive latent activation
    rl[13] = xi*LP;                 T[:,13]=[0,0,1,0,-1,0,0]  #productive latent activation
    rl[14] = k*E*AU;                T[:,14]=[0,-1,0,0,0,0,0]  #unproductive active immune removal
    rl[15] = k*E*AP;                T[:,15]=[0,0,-1,0,0,0,0]  #productive active immune removal
    rl[16] = w*E*(AP+AU)/(E+E50);   T[:,16]=[0,0,0,0,0,1,0]  #immune cell recruitment
    rl[17] = aE;                    T[:,17]=[0,0,0,0,0,1,0]  #immune cell birth
    rl[18] = dE*E;                  T[:,18]=[0,0,0,0,0,-1,0]  #immune cell clearance
    rl[19] = g*V;                   T[:,19]=[0,0,0,0,0,0,-1]  #innate viral clearance
    
    return rl,T

#function that solves stochastically using tau-leap method
def tauleap_simulate(X0,t,inc_time,v_ps,ab_ps):

    dt=t[1]; x=X0; y=[] #initialize
    inf_flag=0 #flag to add the infected cell at inc_time
    
    for ti in t:
        
        y.append(x) #the list of states
        
        rl,T = update_rates(x,ti,v_ps,ab_ps) #make new rate vector
        
        events = np.random.poisson(rl*dt) #calculate events
        
        dx = np.sum(T*events,1) #calculate change in state
        
        x=x+dx #update state variable
        
        x[x<0]=0 #make sure no negative numbers
        
        #add a single infected cell
        if ti>=inc_time and inf_flag==0:
            x[2] = x[2] + 1 #add 1 active productively infected cell at the contact time
            inf_flag=1
               
    return np.array(y)
    
# parameters for viral dynamics model
vol = 1e6      # volume of blood [uL]
aS  = 70*vol;   #constant growth rate of susceptibles [cells/uL/day]
dS  = 0.3;   #susceptible death rate [1/day] 
Bt0 = 1e-4/vol  # infection rate of T-cells [uL/cells-day]/[uL]
dA  = 0.8       # active death rate [1/day]
p   = 5e4       # burst rate of virus from cells [virions/cell]
g   = 23        # virus clearance rate [1/day]
tau = 0.07      # productive infection probability []
lam = 1e-4      # latency probability []
thL = 5.2e-4    # latent clearance rate [1/day]
aL  = 0.015;    # latent proliferation rate [1/day] (Tcm)
xi  = 1e-5;     # latent activation rate [1/day]
dL  = aL-thL-xi # latent death rate
k   = 1/vol;  #immune cell killing rate [uL/cell-day]/[uL]
w   = 2.9;     #immune cell multiplier [1/day]
aE  = 1e-4*vol;   #initial E cell concentration [cells/uL]*[uL]
dE  = 0.002;  #immune death rate [1/day]
E50 = 250*vol;   #50 pct max E cell concentration [cells/uL]*[uL]

R0=aS*Bt0*tau*(1-lam)*p/g/dS/dA; print(R0) # basic reproductive number

t0=1 #time of crossing mucosal barrier
Y1=100
Y2=1
k1=1
k2=0.1
h=1
IC50=10

t=np.linspace(0,500,1e4) #time for antibody dose

X0=np.round([aS/dS,0,0,0,0,aE/dE,0]) #always true because virus gets added at inc_time

v_ps = aS,dS,Bt0,tau,lam,dA,thL,xi,k,aE,dE,E50,w,p,g #viral dynamics parameters
ab_ps = Y1,Y2,k1,k2,h,IC50 #bnab parameters

num_sims=2

plt.figure(figsize=(7,4),dpi=600)
for i in range(num_sims):
    sol=tauleap_simulate(X0,t,t0,v_ps,ab_ps)
    vll=np.log10(sol[:,6]/vol*1e3+1e-3) #viral load as usual units copies/mL

    plt.subplot(231)
    plt.plot(t/7,sol[:,0],'gray',lw=2)
    plt.ylabel('susceptible \n cells per $\mu$L')

    plt.subplot(232)
    plt.plot(t/7,sol[:,1],'red',lw=2)
    plt.plot(t/7,sol[:,2],'orange',lw=2)
    plt.ylabel('actively infected \n cells per $\mu$L')
    plt.legend(['$U$','$P$'],fontsize=8)

    plt.subplot(233)
    plt.plot(t/7,sol[:,3]*1000,'blue',lw=2) #divide by 1000 CD4+ T cells and then multiply by 1 million
    plt.plot(t/7,sol[:,4]*1000,'teal',lw=2)
    plt.ylabel('latently infected \n cells (IUPM CD4+)')
    plt.legend(['$U$','$P$'],fontsize=8)

    plt.subplot(234)
    te=0
    It=sol[te:,1]+sol[te:,2]
    Et=sol[te:,5]
    plt.plot(t[te:]/7,Et*k/dA,'green',lw=2)
    plt.ylabel('adaptive immunity \n (killing ratio)')

    plt.subplot(235)
    plt.plot(t/7,vll,'purple',lw=2)
    plt.axhline(np.log10(30),color='k')
    plt.annotate('detect lim',xy=[10,0],fontsize=8)
    plt.ylabel('viral load \n log10(copies/mL)')
    plt.yticks(range(1,8,1))

    plt.subplot(236)
    plt.plot(t/7,np.cumsum(sol[:,2])/1e6,'tan',lw=2)
    plt.plot(t/7,t/7*300/52,'k--',lw=1)
    plt.ylabel('diversity (arbitrary units)')
    #plt.yticks(range(1,8,1))

    for i in range(6):
        plt.subplot(231+i)
        if i>2:
            plt.xlabel('time (weeks)')
            plt.xticks(range(1,10,1))
        else:
            plt.xticks(range(1,10,1),[])
        
plt.tight_layout()
plt.savefig('sto_model.pdf')
