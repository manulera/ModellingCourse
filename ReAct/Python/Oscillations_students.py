import numpy as np
from Gilles import *
import matplotlib.pyplot as plt
from ColorLine import *

# Initial conditions
user_input = ['DNA', 1,
'DNAP2', 0,
'mRNA', 0,
'Prot', 0,
'Prot2', 0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (10,10)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (

(-1,'DNA'), (1, 'mRNA'), 100 , #This is the only thing that I changed, we almost had it!

(-1, 'DNAP2'), (1,'mRNA'), 0.005,

(1,'mRNA'),(), 50,

(-1, 'mRNA'), (1, 'Prot'), 40,

(1, 'Prot'), (), 20,

(2, 'Prot'), (1,'Prot2'), 3,

(1, 'Prot2'), (2,'Prot'), 1,

(1, 'Prot2',1,'DNA'), (1,'DNAP2'), 0.3,

(1, 'DNAP2'), (1, 'Prot2',1,'DNA'), 0.1,

)
# dt is used for the deterministic calculation, and the
dt=0.01
t = np.arange(0, 60, dt)

(solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t)

Gillesplot(solution,t,tgill, valsgill,rows,mode)
f, ax = plt.subplots()
colorline(valsgill[0][rows['Prot'],:], valsgill[0][rows['mRNA'],:], None, cmap=plt.get_cmap('jet'), linewidth=2)

ax.plot(solution[:,rows['Prot']],solution[:,rows['mRNA']])

plt.show()