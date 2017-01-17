# A very simple script to prove that two functions x(t) y(t), if they are periodic, they produce a closed curve when we
# plot (x,y). This was useful as a proof of concept for the plotting of the predator-pray plots.

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random

# Try other periodic functions! What happens if you include absolute values or squares?

x = np.linspace(0,2*np.pi)

A = np.absolute(np.sin(x))

B = np.cos(x)

plt.figure()

plt.plot(x,A)
plt.plot(x,B)

plt.figure()
plt.plot(A,B)
plt.axis('equal')

plt.show()