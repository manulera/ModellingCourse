import random
from matplotlib import pyplot as plt
import numpy as np
import scipy
from scipy.integrate import odeint
from matplotlib import interactive
import sys
import os

k = 1
n = 100
n0=n
dt = 0.001
plt.ion()
t = np.arange(0,10,dt)
user_decision=True

def step(n0,t,k):
    return -n0*k
print(sys.path)
solution = odeint(step,n0,t,args=(k,))

ax=plt.gca()
ax.plot(t,solution)

def Guill1(n,k,dt):
    listn = list()
    i=0
    while n>0:
        i+=1
        if n*k*dt>random.random():
            listn.append(i)
            n-=1
    y=np.arange(100,0,-1)
    t_disc=np.array(listn)* dt
    ax.plot(t_disc,y)
    plt.show()
    plt.pause(0.01)
while user_decision==True:
    Guill1(n,k,dt)
    user_decision = raw_input("want more?")=="yes"