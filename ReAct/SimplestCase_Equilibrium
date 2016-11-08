import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random

A = 100
B = 0

c1 = 0.01
c2 = 0.01

# 1 step

all_z=[0]

all_A=[100]

for i in range(1000):
    a1 = A*c1
    a2 = B*c2

    a0 = a1 + a2

    z = -1/a0*np.log(random.random())

    r2=random.random()

    if r2*a0 < a1:
        mu = 1
        A = A-1
        B = B+1
    else:
        mu = 2
        A = A + 1
        B = B - 1




