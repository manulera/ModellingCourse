
import numpy as np
from Gilles import *
import matplotlib.pyplot as plt
from DeviationAnalysis import *
from mpl_toolkits.mplot3d import Axes3D
# Initial conditions
user_input = ['A', 100,
              'B', 0]
# Constants (this is not necessary, they could be filled up already in the reaction tuple)
k = (10,10)

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A'),(1,'B'),k[0],
    (1,'B'),(1,'A'),k[1],
)
# dt is used for the deterministic calculation, and the
dt=0.0001
t = np.arange(0, 0.6, dt)

(solution,(tgill, valsgill,all_mus, all_taus),rows,mode) = ReAct(user_input,reactions,t,rounds=300)
fig = plt.figure()

Gillesplot(solution,t,tgill, valsgill,rows,mode)

for i in np.arange(0,0.5,0.05):

    A,X,Y = EquationMaker(reactions,tgill,all_mus, all_taus,i,i+0.05)
    print A
    Y,X=np.meshgrid(Y,X)

    fig = plt.figure()
    #ax = fig.gca(projection='3d')

    #ax.plot_surface(X,Y,A, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
    plt.imshow(A, cmap='hot')

    plt.draw()

plt.show()