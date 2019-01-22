
#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import pylab as pl
from timeit import default_timer as timer

# parameters for viral dynamics model
aS  = 350    # constant growth rate of suseceptible cells
dS  = 0.3    # death rate of suseceptible cells
Bt0 = 1e-4   # infection rate of T-cells
dA  = 1      # death rate of infected cells
pi  = 1e3    # burst rate of virus from cells
gam = 18     # virus clearance rate

n=pi/gam

#antibody parameters
cB=1000 #initial concentration in blood/plasma
cL=100   #initial number in lymph
dAb=.1  #natural+viral clearance of Ab (same in L and B?)
TAb=0.01 #transition rate for Ab from B->L

hill=1
halfmax=10

t=np.linspace(0,100,1e3)

plt.figure()

#dose response curve
dose=np.logspace(-2,3)
response = 1/(1+(dose/halfmax)**-hill)
plt.subplot(121)
plt.semilogx(dose,response)
plt.xlabel('dose (nAbs per $\mu$L)')
plt.ylabel('response (proportion inhibited)')

#inhibition over time
dose=cB*np.exp(-(dAb+TAb)*t)+cL*np.exp(-dAb*t)
prop_inhib = 1/(1+(dose/halfmax)**-hill)

plt.subplot(122)
plt.semilogy(t/7,prop_inhib)
plt.xlabel('time (weeks)')
plt.ylabel('proportion inhibited')

plt.tight_layout()

#set of diff eqs that describes immunologic dynamics simple
def model(X,t):          
    Y=np.zeros(3); #intialize ode solution vector
    
    Y[0] = aS - dS*X[0] - Bt0*X[0]*X[2];
    Y[1] = Bt0*X[0]*X[2] - dA*X[1];
    Y[2] = pi*X[1] - gam*X[2];
    
    return Y   # for odeint

#set of diff eqs that describes immunologic dynamics including bnAbs
def model_bnAbs(X,t):          
    Y=np.zeros(3); #intialize ode solution vector

    dose=cB*np.exp(-(dAb+TAb)*t)+cL*np.exp(-dAb*t)

    t_ode.append(t)
    
    bnAb_dose.append(dose)
    
    prop_inhib = 1/(1+(dose/halfmax)**-hill)

    bnAb_eff.append(prop_inhib)
    
    Bt=Bt0*(1-prop_inhib)
    
    Y[0] = aS - dS*X[0] - Bt*X[0]*X[2];
    Y[1] = Bt*X[0]*X[2] - dA*X[1];
    Y[2] = pi*X[1] - gam*X[2];
    
    return Y   # for odeint

bnAb_dose=[]
bnAb_eff=[]
t_ode=[]

t=np.linspace(0,100,1e3)

S0=aS/dS
V0=0.03
I0=n*V0

RES = spi.odeint(model_bnAbs,[S0,I0,V0],t) 

R0=aS*pi*Bt0/(dA*dS*gam) #calculate basic reproducive number (~8-11 for HIV)

print 'Basic reproductive number is', R0

plt.subplot(121)
plt.semilogy(t/7,RES[:,2]*1e3)
plt.xlabel('time (weeks)')
plt.ylabel('virus concentration (per mL)')
plt.subplot(122)
plt.semilogy(np.array(t_ode)/7,np.array(bnAb_dose))
plt.xlabel('time (weeks)')
plt.ylabel('bnAb concentration (per $\mu$L)')
plt.tight_layout()


rates = aS,dS,Bt0,dA,pi,gam

num_rules=5
num_states=3

prop = np.zeros((num_rules)) #vector of propencities
T = np.zeros((num_rules,num_states)) #vector of transition rules

change_names = ['S(+)','S(-)','I(+)','I(-)','V(-)']

def update_propencities(x,t):
    prop[0] = rates[0]           #constant production 
    prop[1] = rates[1]*x[0]      #density dependent susceptible death
    prop[2] = rates[2]*x[0]*x[2] #infection
    prop[3] = rates[3]*x[1]      #infected cell burst
    prop[4] = rates[5]*x[2]      #density dependent viral clearance
    return prop

def update_transition_mat(x,t):    
    T[0][0] =  1; T[0][1] =  0; T[0][2] =  0;
    T[1][0] = -1; T[1][1] =  0; T[1][2] =  0;
    T[2][0] = -1; T[2][1] =  1; T[2][2] = -1;
    T[3][0] =  0; T[3][1] = -1; T[3][2] =  np.random.poisson(pi);
    T[4][0] =  0; T[4][1] =  0; T[4][2] = -1;
    return T

#stochastic HIV model inspired by Keeling/Rohani+Soumpasis
change_counter = np.zeros([num_rules])

def stochastic_eqs(X,t): 
    
    Rate = update_propencities(X,t) #matrix of rates
    Change = update_transition_mat(X,t) #matrix of transitions
    
    R1=pl.rand(); R2=pl.rand(); #2 uniform random numbers
    
    dt = -np.log(R2)/(np.sum(Rate)); #next time step
    
    m=min(pl.find(pl.cumsum(Rate)>=R1*pl.sum(Rate))); #get the index of the transition
    
    X+=Change[m,:]
    
    change_counter[m]+=1
    
    return np.array(X),dt

def run_stochastic(X0,t_list):
    
    t=0; t_track=[0]; t_index=0; #initialize step, vectors for times and states
    
    X=X0; X_track=[X0]
   
    step=0; 
    while t < max_time:
        
        X,dt = stochastic_eqs(X,t) #calculate the gillespie step and new state
        
        t+=dt 
        
        #just keep track of daily counts
        if t>t_list[t_index]:
            t_track.append(t)
            X_track.append(X); #append new states
            t_index+=1
    
        step+=1
    
    return np.array(t_track),np.array(X_track)

start = timer()

max_time = 10 #day to stop simulation

tPts=100 #how many timepoints to actually tally up

t_list=np.linspace(0,max_time,tPts)
X0=np.array([np.round(aS/dS),10,0])

tv,Xv = run_stochastic(X0,t_list) #run the simulation

end = timer()
print(end - start)      

plt.figure(figsize=(7,3),dpi=600)
state_names = ['S','I','V']

RES = spi.odeint(model,X0,t_list) 

for i in range(3):
    plt.subplot(131+i)
    plt.plot(tv,Xv[:,i])
    plt.plot(t_list,RES[:,i])
    plt.title(state_names[i])
    #plt.xticks(np.arange(0,))
    if i==0:
        plt.ylabel('Concentration per $\mu$L')
    if i==1:
        plt.xlabel('Time (days)')
plt.yscale('log')
plt.tight_layout()

plt.bar(np.arange(5)-0.5,change_counter)
plt.yscale('log')
plt.xticks(range(5),change_names)
plt.tight_layout()


# In[226]:

len(tv)


# In[ ]:



