
import numpy as np
from Gilles import *
import matplotlib.pyplot as plt


# Initial conditions
user_input = ['A', 1000,
              'B', 0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (10,10)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A'),(1,'B'),k[0],
    (1,'B'),(1,'A'),k[1],
)
# dt is used for the deterministic calculation, and the
dt=0.0001
t = np.arange(0, 2, dt)

(solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t)

Gillesplot(solution,t,tgill, valsgill,rows,mode)

plt.show()