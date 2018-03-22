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
    b = []
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        b += [l[i:i+n]]

    return b

def delta_G_T(file,T,L):
    Energy = np.loadtxt(file)
    E_Section = Energy[0:T]
    t_mean = np.mean(E_Section)
    Chunks = E_Section.reshape((T/L),L)
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