
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

num=5

#set of diff eqs that describes immunologic dynamics simple
def sde_model:          
    
    S0=aS/dS; V0=0.03; I0=n*V0
    
    t=np.linspace(0,100,1e3)
    
    X0=np.repmat([S0,I0,V0],1,num); #intialize ode solution vector
    
    for dt in t:
        
    Y[0] = aS - dS*X[0] - Bt0*X[0]*X[2];
    Y[1] = Bt0*X[0]*X[2] - dA*X[1];
    Y[2] = pi*X[1] - gam*X[2];
    
    return Y   # for odeint




