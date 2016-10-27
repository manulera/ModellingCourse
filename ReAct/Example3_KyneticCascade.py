
import numpy as np
from Guilles import ReAct
import matplotlib.pyplot as plt

mode=0
# Initial conditions
user_input = ['A', 500,
              'B', 500]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (10,10)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A'),(1,'B'),k[0],
    (1,'B'),(1,'A'),k[1],
)
# dt is used for the analytical calculation, and the
dt=0.00001
t = np.arange(0, 1, dt)

(solution,(tguill, valsguill),names)=ReAct(user_input,reactions,t,mode)

fig, ax = plt.subplots()
tguill=np.array(tguill)
for i in range(valsguill.shape[0]):
    ax.plot(t,solution[:,i])
    ax.step(tguill, valsguill[i, :])

# Show over time
plt.show()