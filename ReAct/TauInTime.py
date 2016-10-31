import numpy as np
from Gilles import *
import matplotlib.pyplot as plt
from DeviationAnalysis import *
#What kind of information can we extract in time out of this?

# Initial conditions
user_input = ['A', 200]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = 0.5

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A'),(),k,
)
# dt is used for the deterministic calculation, and the
dt=0.1
t = np.arange(0, 10, dt)
SolList=list()
(solution, (tgill, valsgill), rows, mode) = ReAct(user_input, reactions, t)
Gillesplot(solution,t,tgill, valsgill,rows,mode)
ax = plt.gca()
for i in range(20):
    (solution,(tgill, valsgill),rows,mode)=ReAct(user_input,reactions,t)
    diff_tau = np.diff(tgill[0])
    ax.plot(np.squeeze(tgill)[1:], diff_tau)
    print diff_tau[-1]
    print

ax=plt.gca()



plt.show()



