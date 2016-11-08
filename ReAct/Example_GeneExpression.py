
import numpy as np
from Gilles import *
import matplotlib.pyplot as plt


# Initial conditions
user_input = ['Isite1', 1,
              'Isite1_b', 0,
              'mRNA', 0,
              'Prot', 0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (10,10)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (

    (-1,'Isite1'), (1, 'mRNA'), 0.2 ,

    (1,'mRNA'),(),0.01,

    (-1, 'mRNA'), (1, 'Prot'), 2.05,

    (1, 'Prot'), (), 4.0125,

    #(1, 'Prot2'), (), 4.0125,

    #(2, 'Prot'), (1,'Prot2'), 2,

    #(1,'Prot2'), (2,'Prot'), 0.2,

    (1, 'Prot',1,'Isite1'), (1,'Isite1_b'), 1,

    (1, 'Isite1_b'), (1, 'Prot',1,'Isite1'), 1,
)
# dt is used for the deterministic calculation, and the
dt=1
t = np.arange(0, 10000, dt)

(solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t,mode=1)

Gillesplot(solution,t,tgill, valsgill,rows,mode)

plt.show()