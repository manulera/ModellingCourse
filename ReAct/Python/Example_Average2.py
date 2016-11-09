
import numpy as np
from Gilles import *
import matplotlib.pyplot as plt
from DeviationAnalysis import *

# Initial conditions
user_input = ['A', 100,
              'B', 0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (0.1,0.1)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (2,'A'),(2,'B'),k[0],
    (2,'B'),(2,'A'),k[1],
)
# dt is used for the deterministic calculation, and the
dt=0.0001
t = np.arange(0, 0.5, dt)

fig, ax = plt.subplots()
for i in range(4):
    (solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t,rounds=50)

    curves = CurveAverager(tgill,valsgill,t)

    DeviPlot=GillDeviation(solution,curves)


    ax.plot(t,solution[:,1])
    ax.plot(t,solution[:,0])
    ax.plot(t,curves[-1,:,1])
    ax.plot(t,curves[-1,:,0])
    ax.plot(t,DeviPlot[:,0]*10)
    axes = plt.gca()

axes.set_ylim([0,100])
plt.show()

# fig, ax = plt.subplots()
#
# for i in range(1):
#     print i
#     (solution,(tgill, valsgill),rows,mode)=ReAct(user_input,reactions,t,rounds=100)
#     curves = CurveAverager(tgill,valsgill,t)
#     DeviPlot=GillDeviation(solution,curves)
#     ax.plot(t,DeviPlot[:,0])
#
# axes = plt.gca()
# axes.set_ylim([0,10])
# plt.show()