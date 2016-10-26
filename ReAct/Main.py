# Trying to make everything from the reactions only
# Kinases Model II
# Time

from DeterministicFunctions import Stochy2
import numpy as np
import matplotlib.pyplot as plt
dt = 0.01
t = np.arange(0, 100, dt)

# Parameters

user_input = ['A', 100,
              'B', 0]

k = (0.2,0.1)

# Empty reaction ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A'),(1,'B'),k[0],
    (1,'B'),(1,'A'),k[1],
)

(solution, names) = Stochy2(user_input, reactions, t)

fig, ax = plt.subplots()

for i in [0,1]:
    ax.plot(t, solution[:, i], label=names[i])
legend = ax.legend(loc='center', shadow=True)
# Show over time
plt.show()