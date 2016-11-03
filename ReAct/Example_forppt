
import numpy as np
from Gilles import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pprint import pprint

# Initial conditions
user_input = ['A', 150,
              'B', 100,
              'C', 0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (0.05,0.005)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A',1,'B'),(2,'C'),k[0],
    (2, 'C'),(1,'A',1,'B'),k[1],
)
# dt is used for the deterministic calculation, and the
dt=0.0001
t = np.arange(0, 1, dt)

(solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t)

Gillesplot(solution,t,tgill, valsgill,rows,mode)

ax=plt.gca()

print ax.lines
for i in ax.lines:
    i.set_linewidth(3.0)
plt.draw()
plt.savefig('Exampleppt2.png')