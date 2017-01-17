import numpy as np
from Gilles import *
import matplotlib.pyplot as plt
from ColorLine import *
# Initial conditions
user_input = ['DNA', 1,
'DNAP2', 0,
'mRNA_n', 0,
'mRNA_c', 0,
'Prot', 0,
'Prot2_c', 0,
'Prot2_n', 0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (10,10)
# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
(-1,'DNA'), (1, 'mRNA_n'), 100,
(-1, 'DNAP2'), (1,'mRNA_n'), 0.05,
(1, 'mRNA_n'),(), 1,

(1, 'mRNA_c'), (), 10,
(-1, 'mRNA_c'), (1, 'Prot'), 400,
(1, 'Prot'), (), 10,
(2, 'Prot'), (1,'Prot2_c'), 0.1,
(1, 'Prot2_c'), (2,'Prot'), 1000,

(1, 'mRNA_n'), (1, 'mRNA_c'), 5,

#(1, 'mRNA_c'), (1, 'mRNA_n'), 5,

(1, 'Prot2_c'), (1, 'Prot2_n'), 100,

(1, 'Prot2_n'), (1, 'Prot2_c'), 10,
(1, 'Prot2_n', 1,'DNA'), (1,'DNAP2'), 10,
(1, 'DNAP2'), (1, 'Prot2_n',1,'DNA'), 1,
)
# dt is used for the deterministic calculation, and the
dt=0.01
t = np.arange(0, 20, dt)
(solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t, mode = 1)

solution[:,rows['DNA']]=solution[:,rows['DNA']]*100

Gillesplot(solution,t,tgill, valsgill,rows,mode)
f, ax = plt.subplots()
#colorline(valsgill[0][rows['Prot'],:], valsgill[0][rows['mRNA_c'],:], None, cmap=plt.get_cmap('jet'), linewidth=2)
ax.plot(solution[:,rows['Prot']],solution[:,rows['mRNA_c']])
plt.show()