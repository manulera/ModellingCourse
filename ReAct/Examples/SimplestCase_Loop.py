import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random

A = 100
c = 0.01

# All steps in the reaction, show the concept of

all_z=[0]

all_A=[100]


while A > 0:
    a1 = A*c
    a0 = a1

    z = -1/a0*np.log(random.random())
    A = A-1

    all_z.append(z)
    all_A.append(A)

plt.figure()

all_z=np.array(all_z)
all_A=np.array(all_A)

plt.step(np.cumsum(all_z),all_A)

plt.show()



