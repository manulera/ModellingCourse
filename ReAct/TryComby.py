import numpy as np

Gamma = np.array([[-2, 2,],
                  [-1, 1,],
                  [ 3,-3]], int)


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
decay=False
def calcP(init, Gamma, k):
    """
    Will calculate the propensity of all reactions
    :param init:
    :param Gamma:
    :param k:
    :return:
    """
    P = list()
    # Write the reaction and create V
    for i in range(0, Gamma.shape[1]):  # Along the reactions

        Pi = 1
        Ri = Gamma[:, i]
        for j in range(0, Gamma.shape[0]):  # Along the reactants
            if Ri[j] < 0:
                print 'hi'
                # Ri[j] is the stochiometry index of the element j
                # init[j] is the number of molecules of element j that we have
                Pi *= Comby(init[j], -Ri[j])
            elif decay:
                for ele in decay:
                    if (ele[0], ele[1]) == (j, i):
                        Pi *= Comby(init[j], ele[2])
        P.append(Pi * k[i])
    return (P,np.cumsum(P))

a = [0, 1, 2, 3, 4]
print next(item for item in a if item>1)