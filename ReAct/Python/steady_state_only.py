
import numpy as np
from Gilles import *
import matplotlib.pyplot as plt
from DeviationAnalysis import *

# Initial conditions
user_input = ['A', 150,
              'B', 100,
              'C', 0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (0.5,0.05, 0.1, 0.3)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A',),(1,'B'),k[0],
    (1,'B',),(1,'A'),k[1],
    (1,'B',),(1,'C'),k[2],
    (1,'C',),(1,'A'),k[3],
)
print ReAct(user_input,reactions,0,mode=3)
t = np.linspace(0, 20)

(solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t)

print solution[-1,:]
Gillesplot(solution,t,tgill, valsgill,rows,mode)
plt.show()
