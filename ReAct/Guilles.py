
import numpy as np
import math
import random
import matplotlib.pyplot as plt
# Parameters

user_input = ['A', 100,
              'B', 0]

k = (0.2,0.1)

# Empty reaction ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A'),(1,'B'),k[0],
    (1,'B'),(1,'A'),k[1],
)

#Function to calculate the possible combinations?
dt=0.001
t = np.arange(0, 10, dt)

def Stochy2(user_input, reactions, t):
    # Create elements and rows
    row = dict()
    elements = list()
    init = list()
    for i in range(0, len(user_input), 2):
        elements.append(user_input[i])
        row.update({user_input[i]: i / 2})
        init.append(user_input[i + 1])

    # Create Gamma and decay
    Gamma = np.zeros((len(elements), len(reactions) / 3), int)
    decay = tuple()
    k = tuple()
    for j in range(0, len(reactions), 3):  # Reactants/ +1 Products/ +2 Konstants
        for ind in range(0, len(reactions[j]), 2):  # In the reactants, get the value for Gamma, and maybe decay
            sto = reactions[j][ind]
            name = reactions[j][ind + 1]
            if sto > 0:
                Gamma[row[name], j / 3] = -sto;
            else:
                decay += ((row[name], j / 3, abs(sto)),)
        for ind in range(0, len(reactions[j + 1]), 2):  # In the products, get the value for Gamma
            sto = reactions[j + 1][ind]
            name = reactions[j + 1][ind + 1]
            if sto > 0:
                Gamma[row[name], j / 3] = sto;
        k += (reactions[j + 2],)

    # Call Stochy
    return (Stochy(elements, init, t, Gamma, k, decay), row)

def Stochy(elements, init, t, Gamma, k, decay):

    def Comby(quant, index):
        """
        Input is quant (number of molecules) and index (the reaction index for that particular molecule in the reaction).
        A very interesting question is why this function will always return an int, think about it!
        :param quant: int
        :param index: int
        :return:
        """
        result = 1
        for i in range(index):
            result *= (quant - i)
            result /= (i + 1)
        return result

    def calcP(init, Gamma, k):
        """
        Will calculate the cummulative sum of propensities for all reactions
        :param init:
        :param Gamma:
        :param k:
        :return:
        """
        P = list()
        # Calculate Pi, propensity of each reaction
        for i in range(0, Gamma.shape[1]):  # Along the reactions

            Pi = 1
            Ri = Gamma[:, i]
            for j in range(0, Gamma.shape[0]):  # Along the reactants
                if Ri[j] < 0:
                    # Ri[j] is the stochiometry index of the element j
                    # init[j] is the number of molecules of element j that we have
                    Pi *= Comby(init[j], -Ri[j])
                elif decay:
                    for ele in decay:
                        if (ele[0], ele[1]) == (j, i):
                            Pi *= Comby(init[j], ele[2])
            P.append(Pi * k[i])
        return np.cumsum(P)

    # ODE system (in format used for odeint)
    def GuilleStep(init, t, Gamma, k):

        P = calcP(init, Gamma, k)
        # P[-1] is what Guillespie calls a_0
        tau_step=(math.log(1)-math.log(random.random()))/P[-1]
        # Find the first non-false element in the logical array resulting from > comparison
        mu=next(i for i, x in enumerate(P > random.random() * P[-1]) if x)
        # One round of the reaction results in the substraction of the column of Gamma
        # corresponding to the reaction from the init values
        init -= Gamma[:, mu]
        t+=tau_step
        return (init, t)

    tcount=0
    tend=t[-1]
    tguill= [0]
    valsguill= np.array(init)
    while tcount < tend:
        (init, tcount)=GuilleStep(init, tcount, Gamma, k)
        valsguill=np.c_[valsguill,np.array(init)]
        tguill.append(tcount)

    return tguill,valsguill

((tguill, valsguill),names)=Stochy2(user_input,reactions,t)

print valsguill.shape
fig, ax = plt.subplots()
tguill=np.array(tguill)
for i in range(valsguill.shape[0]):

    ax.plot(tguill, valsguill[i, :])

# Show over time
plt.show()