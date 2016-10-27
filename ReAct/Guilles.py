
import numpy as np
import math
import random
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def AnaSol(elements, init, t, Gamma, k, decay):

    # Function to calculate V
    def calcV(init, Gamma, k):
        V = list()
        # Write the reaction and create V
        for i in range(0, Gamma.shape[1]):  # Along the reactions
            Vi = 1
            Ri = Gamma[:, i]
            for j in range(0, Gamma.shape[0]):  # Along the reactants
                if Ri[j] < 0:
                    Vi *= (init[j] ** abs(Ri[j]))
                elif decay:
                    for ele in decay:
                        if (ele[0], ele[1]) == (j, i):
                            Vi *= (init[j] ** ele[2])
            V.append(Vi * k[i])
        return V

    # ODE system (in format used for odeint), this function calculates and returns, for element Ri in init, the dRi/dt
    def myODE(init, t, Gamma, k):

        # t is not used inside the function, but its required as an input for the calculations
        V = calcV(init, Gamma, k)
        rates = np.dot(Gamma, V)
        return tuple(rates) # rates is d(init)/dt

    # Solve using odeint
    solution = odeint(myODE, init, t, args=(Gamma, k,))
    return solution

def Guillespy(elements, init, t, Gamma, k, decay):

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
            result *= (quant - i)  # Number of possible collisions, if we would include the line on the bottom
            # result /= (i + 1) This would be needed if the kynetic constant did not have this information inside already
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
    def GuilleStep(init, t, Gamma, k, tend):

        P = calcP(init, Gamma, k)
        if P[-1]==0:
            return (init,tend)
            # The system is dead!
        else:
            # P[-1] is what Guillespie calls a_0
            tau_step=(math.log(1)-math.log(random.random()))/(P[-1])
            # Find the first non-false element in the logical array resulting from > comparison
            mu=next(i for i, x in enumerate(P > random.random() * P[-1]) if x)
            # One round of the reaction results in the substraction of the column of Gamma
            # corresponding to the reaction from the init values
            init += Gamma[:, mu]
            t+=tau_step
            return (init, t)

    tcount=0
    tend=t[-1]
    tguill= [0]
    valsguill= np.array(init)
    while tcount < tend:
        (init, tcount)=GuilleStep(init, tcount, Gamma, k, tend)
        valsguill=np.c_[valsguill,np.array(init)]
        tguill.append(tcount)
    return np.array(tguill),valsguill

def ReAct(user_input, reactions, t, mode=0):
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
                Gamma[row[name], j / 3] = -sto
            else:
                decay += ((row[name], j / 3, abs(sto)),)
        for ind in range(0, len(reactions[j + 1]), 2):  # In the products, get the value for Gamma
            sto = reactions[j + 1][ind]
            name = reactions[j + 1][ind + 1]
            if sto > 0:
                Gamma[row[name], j / 3] = sto
        k += (reactions[j + 2],)

    # Write the reaction
    for i in range(0, Gamma.shape[1]):  # Along the reactions
        reactants = list()
        products = list()
        Ri = Gamma[:, i]
        for j in range(0, Gamma.shape[0]):  # Along the reactants
            if Ri[j] < 0:
                reactants.append(str(abs(Ri[j])) + elements[j])
            elif Ri[j] > 0:
                products.append(str(abs(Ri[j])) + elements[j])
            elif decay:
                for ele in decay:
                    if (ele[0], ele[1]) == (j, i):
                        reactants.append(str(ele[2]) + elements[j])
                        products.append(str(ele[2]) + elements[j])
        print reactants, '--k' + str(i) + '-->', products

    # Different return depending on none
    if mode==0: # mode 0: Guillespie and Analytical solution
        return (AnaSol(elements, init, t, Gamma, k, decay), Guillespy(elements, init, t, Gamma, k, decay), row, mode)
    elif mode==1:# mode 1: Analytical solution only
        return (AnaSol(elements, init, t, Gamma, k, decay), (None,None), row, mode)
    else: # mode 2: Guillespie only
        return (None, Guillespy(elements, init, t, Gamma, k, decay), row, mode)

def Guillesplot(solution,t,tguill, valsguill,rows,mode,which2plot=False):

    # If there is not a specific list, plot all of them
    if not which2plot: which2plot=rows.keys()
    fig, ax = plt.subplots()

    for i in which2plot:
        color0=(random.random()/1.1,random.random()/1.1,random.random()/1.1) #Divided by 2 not to get too light colors
        if mode != 2:
            ax.plot(t, solution[:, rows[i]], color=color0,label=i)
        if mode != 1:
            ax.step(tguill, valsguill[rows[i], :], color=color0)

    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
              fancybox=True, shadow=True, ncol=5)

    plt.draw()