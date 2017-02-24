__author__ = 'Maggie'

import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import math
import sys

def get_avg(file):

    S = np.loadtxt(file)

    set = []
    for line in S:
        set += [line[1]]

    avg = sum(set)/len(set)
    return avg



def chunks(l, n):
    b =[]
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        b += [l[i:i+n]]

    return b

def delta_G_T(file,T,L):
    Energy = np.loadtxt(file)
    E_Section = Energy[0:T]
    t_mean = np.mean(E_Section)
    Chunks = E_Section.reshape((T/L),L)
    #Chunks = chunks(Energy, L)
    Sum = 0


    for i in range (0,T/L):
        chunk = Chunks[i]
        c_mean = np.mean(chunk)
        square = (c_mean - t_mean)**2
        Sum += square

    R = (L/float(T))*np.sqrt(Sum)
    return R




if __name__ == '__main__':

    file = 'chi_1_250000steps\EnergiesDual_250000a'
    Energy = np.loadtxt(file)
    Tot = len(Energy)
    cut = 40
    L = 1000
    num = Tot/L
    Y = []
    X = []
    for T in range (L, Tot+1, L):
        y = delta_G_T(file, T, L)
        Y += [y]
        X += [float(1)/np.sqrt(T)]
        #X += [T]
        if T%1000 == 0:
            print("completed" + str(T)+ "steps")



    print len(Y)
    print num
    print X

    fit = np.polyfit(X[cut:num],Y[cut:num],1)
    fit_fn = np.poly1d(fit)

    plt.plot(X[cut:num],Y[cut:num],label='Simulation')
    plt.plot(X[cut:num], fit_fn(X[cut:num]), '--k', label='Linear Fit')
    plt.ylabel('RMSE of the Energy')
    plt.xlabel('1/sqrt(T) T=time')
    plt.legend(loc='upper left')
    plt.title('Test for Ergodic Sampling')
    plt.figure(2)
    plt.plot(np.arange(0,Tot),Energy)
    plt.show()



"""    Energy = np.loadtxt('Energies3D_90000a')

    avg_tot = np.mean(Energy)
    print(avg_tot)
    std_tot = np.std(Energy)
    print(std_tot)
    T = len(Energy)
    L = T/100
    print(T)
    print(L)

    index = L
    sum = 0
    count = 0
    stands = []
    adding_sums = []

    while index < T-6*L+1:
        chunk = Energy[index-L:index+1]
        c_mean = np.mean(chunk)
        large_mean = np.mean(Energy[0:index+1])
        stand = np.std(chunk)
        stands += [stand]
        square = (c_mean-large_mean)**2
        sum += square
        count += 1
        STD4 = np.sqrt(sum/count)
        adding_sums += [np.sum(adding_sums)+ STD4]
        index = index +L


    print(sum)
    #print(len(sum))
    print(count)
    print(stands)
    print(adding_sums)

    STD = np.sqrt(sum/100)
    STD2 = (float(L)/T)*(np.sqrt(sum))
    print(STD)
    print(STD2)

   # t = [i for i in range(0,1000)]
   # x = 1./np.sqrt(t)
    #plt.plot(x,adding_sums)
"""


