import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random

from ColorLine import *

init=np.array([80,40])
dt=0.001
t=np.arange(0,10,dt)

# Very nice values to see the oscillations, extinction happens

import numpy as np
from Gilles import *
import matplotlib.pyplot as plt


# Initial conditions
user_input = ['Pred', 40.0,
              'Prey', 160.0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k=[5,0.045,0.005,1,0.002]

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
# dx = A * x - B * x * y
# dy = C * x * y - D * y
reactions = (
    (-1,'Prey'),(1,'Prey'),k[0], # Birth of a prey, would be interesting to see if we change the stoch to 2 for A
    (-1,'Pred',1,'Prey'),(),k[1], # Prey is hunted, analogy is a catalizer for degradation
    (-1,'Pred',-1,'Prey'),(1,'Pred'),k[2], # Predators nourrish on prey
    (1,'Pred'),(),k[3], # Predators die
    (-1,'Prey',1,'Prey'),(),k[4] # If too many preys, they will consume their resources
)
# dt is used for the deterministic calculation, and the
dt=0.01
t = np.arange(0, 40, dt)

(solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t)

print 'calculated'

Gillesplot(solution,t,tgill, valsgill,rows,mode)

f, ax = plt.subplots()

colorline(valsgill[0][0,:], valsgill[0][1,:], None, cmap=plt.get_cmap('jet'), linewidth=2)

ax.plot(solution[:,0],solution[:,1])

plt.show()

