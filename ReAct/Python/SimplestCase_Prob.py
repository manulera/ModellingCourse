import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random

# Draw the P(z,mu)

runs=1000
all_z=np.zeros(runs)

for i in range(runs):
    A = 100
    c = 0.01

    # 1 step

    a1 = A*c
    a0 = a1

    z = -1/a0*np.log(random.random())
    all_z[i]=z

plt.figure()

plt.hist(all_z)

plt.show()


