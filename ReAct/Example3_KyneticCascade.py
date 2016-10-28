
import numpy as np
from Gilles import *
import matplotlib.pyplot as plt

# Initial conditions
user_input = ['Rec', 10,
              '1M3', 10,
              '1M3P', 0,
              '1M2', 20,
              '1M2P', 0,
              '1M1', 30,
              '1M1P', 0]

# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (2,0.05,1,0.5,1,0.5,1)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'Rec'),(),k[0],
    (-1,'Rec',1,'1M3'),(1,'1M3P'),k[1],
    (1,'1M3P'),(1,'1M3'),k[2],
    (-1,'1M3P',1,'1M2'),(1,'1M2P'),k[3],
    (1,'1M2P'),(1,'1M2'),k[4],
    (-1, '1M2P', 1, '1M1'), (1, '1M1P'), k[5],
    (1, '1M1P'), (1, '1M1'), k[6],
)
# dt is used for the deterministic calculation, and the
dt=0.00001
t = np.arange(0, 10, dt)

(solution,(tgill, valsgill),rows,mode)=ReAct(user_input,reactions,t)


Gillesplot(solution,t,tgill, valsgill,rows,mode)

Gillesplot(solution,t,tgill, valsgill,rows,mode,['Rec','1M3P','1M2P','1M1P'])
plt.show()