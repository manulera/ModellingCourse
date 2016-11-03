
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

def EquationMaker(reactions,tgill,all_mus, all_taus,tmin,tmax,nbins=10):

    numreactions = len(reactions)/3
    rounds = len(tgill)
    print tmin, tmax
    mu_pairs=np.array([])
    tau_pairs=np.array([])

    for i in range(rounds): #Along the different solutions

        selected1= np.logical_and(np.greater_equal(tgill[i],tmin), np.less(tgill[i],tmax ))

        mu_pairs=np.hstack((mu_pairs,all_mus[i][np.where(selected1)]))
        tau_pairs=np.hstack((tau_pairs,all_taus[i][np.where(selected1)]))


        # for j in range(numreactions): #Along the different reactions
        #     selected2= np.equal(all_mus[i][np.where(selected1)],j)
        #     print np.sum(selected2)
        #     mu_pairs = np.hstack((mu_pairs,all_mus[i][np.where(selected2)]))
        #     tau_pairs = np.hstack((tau_pairs,all_taus[i][np.where(selected2)]))

    output = np.zeros([nbins, numreactions])
    tau_step = float(np.amax(tau_pairs))/float(nbins)

    for j in range(mu_pairs.size):
        # A bit dirty solution for the maximum value, otherwise, it exceeds the maximum index
        output[int(tau_pairs[j]/tau_step)-(int(tau_pairs[j]/tau_step)>=nbins), mu_pairs[j]] += 1
        #print [int(tau_pairs[j]/tau_step)-(int(tau_pairs[j]/tau_step)>=nbins), mu_pairs[j]]
    return output,np.linspace(0, np.amax(tau_pairs),nbins),np.arange(numreactions)

def GillDeviation(anasol,curves):

    anasol = np.repeat(np.expand_dims(anasol,0),curves.shape[0]-1,0)
    devmatrix=np.fabs(curves[:-1,:,:]-anasol)

    return np.sum(devmatrix,0)/(curves.shape[0]-1)



# C=np.expand_dims(A,2)
# print A.shape, C.shape
# D=np.repeat(C,2,2)
# print A.shape, C.shape, D.shape
