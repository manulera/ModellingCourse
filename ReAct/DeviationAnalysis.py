
import numpy as np
import math
import random
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def CurveAverager(tgill,valsgill,t):

    rounds = len(valsgill)
    curves=np.zeros((rounds+1,len(t),valsgill[0].shape[0]))

    for i in range(rounds): #Along the different solutions
        ind = 0
        j=0
        while j < (t.shape[0]-1): # Along the time points

            if tgill[i][ind] <= t[j]:
                #print tgill[i][ind], t[j], valsgill[i][:,ind]
                ind += 1
                continue
            #print valsgill[i][:, ind]
            curves[i,j,:] = valsgill[i][:,ind]
            j+=1
    # The last curve will be the average

    curves[rounds, :, :]=np.sum(curves,0)/rounds
    return curves

def GillDeviation(anasol,curves):

    anasol = np.repeat(np.expand_dims(anasol,0),curves.shape[0]-1,0)
    devmatrix=np.fabs(curves[:-1,:,:]-anasol)
    return np.sum(devmatrix,0)/(curves.shape[0]-1)


# C=np.expand_dims(A,2)
# print A.shape, C.shape
# D=np.repeat(C,2,2)
# print A.shape, C.shape, D.shape
