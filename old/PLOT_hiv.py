#!/usr/bin/python

# by DBR 9/2016 #
# script to read in and plot the data from the c gillespie SIV model #

import sys
import numpy as np
import matplotlib.pyplot as plt

def read_data(fn='out.txt'):
    # read the file 
    f2 = open(fn, 'r')
    # read the whole file into a single variable, a list of every row of the file.
    lines = f2.readlines()
    f2.close()

    # initialize some variable to be lists:
    t=[]; S=[]; I=[]; V=[];

    # scan the rows of the file stored in lines, and put the values into S,I,V variables
    for line in lines:
        p = line.split()
        t.append(float(p[0]))
        S.append(float(p[1]))
        I.append(float(p[2]))
        V.append(float(p[3]))

    tt=np.array(t)
    St=np.array(S)
    It=np.array(I)
    Vt=np.array(V)

    return tt,St,It,Vt

def plot_data(tt,St,It,Vt):
    plt.subplot(121)    
    plt.step(tt/7,St,color='olivedrab',lw=2,alpha=0.3)
    plt.step(tt/7,It,color='indianred',lw=2,alpha=0.3)
    plt.xlabel('Time (weeks)')
    plt.ylabel(r'Cells per $\mu$L')
    plt.legend(['Susceptible','Infected'])

    plt.subplot(122)    
    plt.semilogy(tt/7,Vt,color='mediumorchid',lw=2,alpha=0.3)
    plt.ylim([1,1e5])
    plt.xlabel('Time (weeks)')
    plt.ylabel(r'Virus per mL')

if len(sys.argv)>=2:
    fn=str(sys.argv[1])

sol=read_data(fn)
plot_data(sol[0],sol[1],sol[2],sol[3])
plt.tight_layout()
plt.show()

plt.tight_layout()
plt.show()

