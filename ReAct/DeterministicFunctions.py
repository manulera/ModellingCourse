# Calling all the libraries and creating all the functions

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def Stochy(elements, init, t, Gamma, k, decay):
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

    # ODE system (in format used for odeint)
    def myODE(init, t, Gamma, k):

        # rates
        V = calcV(init, Gamma, k)
        rates = np.dot(Gamma, V)
        return tuple(rates)

    # Solve using odeint
    solution = odeint(myODE, init, t, args=(Gamma, k,))
    fig, ax = plt.subplots()
    for i in range(0, Gamma.shape[0]):
        ax.plot(t, solution[:, i], label=elements[i])
    legend = ax.legend(loc='center', shadow=True)
    # Show over time
    plt.show()
    return solution


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