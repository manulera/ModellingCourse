import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random

k=[1,0.11,0.01,1]
init=np.array([80,40])
dt=0.001
t=np.arange(0,10,dt)

def myODE(init,t,k):
    x,y=init
    A,B,C,D=k
    dx=A*x-B*x*y
    dy=C*x*y-D*y
    return dx,dy

detsolution = odeint(myODE, init,t,args=(k,))
print 'Deterministic finished'

def StoSol(init,t,k):

    stosol=np.zeros([t.shape[0],2])
    dt = t[1] - t[0]
    for i in range(t.size):
        Ps = myODE(init,t,k)

        booly=[Ps[0]*dt>random.random(),Ps[1]*dt>random.random()]
        init = init + np.array(booly)
        stosol[i,:]=init
    return stosol

# stosolution=StoSol(init,t,k)


plt.figure()

plt.plot(t,detsolution[:,0],label='prey')
plt.plot(t,detsolution[:,1],label='predator')
# plt.step(t,stosolution[:,0],label='prey')
# plt.step(t,stosolution[:,1],label='predator')

plt.legend()
plt.show()