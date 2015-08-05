
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import math
import collections
import itertools
import ast




def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

def last_line(file):
    with open(file) as f:
        #first = f.readline()     # Read the first line.
        #f.seek(-2, 2)            # Jump to the second last byte.
        #while f.read(1) != "\n": # Until EOL is found...
        #    f.seek(-2, 1)        # ...jump back the read byte plus one more.
        #last = f.readline()      # Read last line.
        nums = np.loadtxt(file)
        last = nums[-1]
        f.close()
        return last

def get_avg(file):

    S = np.loadtxt(file, skiprows=1)

    set = []
    for line in S:
        set += [line[1]]

    avg = sum(set)/len(set)
    return avg

def binomial_graph(file):
    bin = open(file, 'r')
    yvals = np.loadtxt(bin, delimiter=',')
    binomial = []
    for k in xrange(0,4):
        n =3
        binomial += [choose(n,k)*(0.5**k)*(0.5**(n-k))]

    t = np.array([0,1,2,3])
    plt.plot(t.T,binomial,'r')
    print yvals
    for i in range(yvals.shape[0]):
        if i == 1:
            plt.plot(t.T,yvals[i],'b')
        elif i == 2:
            plt.plot(t.T,yvals[i],'g')
        elif i == 3:
            plt.plot(t.T,yvals[i],'c')
        else:
            plt.plot(t.T, yvals[i], 'm')

    plt.title(file)
    plt.show()






if __name__ == '__main__':
    energy_10 = []
    for c in "abcd":
        #print [last_line("Energies-10000" + str(c))]
        energy_10 += [last_line("Energies10000" + str(c))]
    energy_50 = []
    for d in "abcd":
        energy_50 += [last_line("Energies50000" + str(d))]
    energy_100 = []
    for e in "abcd":
        energy_100 += [last_line("Energies100000" + str(e))]
    energy_200 = []
    for f in "abc":
        energy_200 += [last_line("Energies200000" + str(f))]
    energy_250 = []
    for g in "abc":
        energy_250 += [last_line("Energies250000" + str(g))]
    energy_300 = []
    for h in "abc":
        energy_300 += [last_line("Energies300000" + str(h))]


    avg10 = sum(energy_10)/len(energy_10)
    avg50 = sum(energy_50)/len(energy_50)
    avg100 = sum(energy_100)/len(energy_100)
    avg200 = sum(energy_200)/len(energy_200)
    avg250 = sum(energy_250)/len(energy_250)
    avg300 = sum(energy_300)/len(energy_300)

    std_10 = np.std(energy_10)
    std_50 = np.std(energy_50)
    std_100 = np.std(energy_100)
    std_200 = np.std(energy_200)
    std_250 = np.std(energy_250)
    std_300 = np.std(energy_300)

    y_vals = [avg10, avg50, avg100, avg200, avg250, avg300]
    x_vals = [10000, 50000, 100000, 200000, 250000, 300000]
    yerr = [std_10, std_50, std_100, std_200, std_250, std_300]
    #plt.errorbar(x_vals, y_vals, yerr=yerr, marker = 'o')
    #plt.xlim(5000, 105000)
    #axarr[0].title('Plot of Energies for each Number of Runs')
    #plt.scatter(x_vals, energy_50, label='Energy for 50,000')
    #plt.scatter(x_vals, energy_100, label='Energy for 100,000')

    #legend = plt.legend(loc='upper right', shadow = True, fontsize = 'x-large')
   # plt.figure(1)
    plt.show()
   # plt.figure(2)
    SSR10 = get_avg("Saved_SSRe")
    SSR50 = get_avg("Saved_SSRf")
    SSR100 = get_avg("Saved_SSRg")
    SSR200 = get_avg("Saved_SSRh")
    SSR250 = get_avg("Saved_SSRi")
    SSR300 = get_avg("Saved_SSRj")

    StandDevA = np.loadtxt("Standard_devA", skiprows=2)
    Stand10 = StandDevA[-1]
    StandDevB = np.loadtxt("Standard_devB", skiprows=1)
    Stand50 = StandDevB[-1]
    StandDevC = np.loadtxt("Standard_devC", skiprows=2)
    Stand100 = StandDevC[-1]
    StandDevD = np.loadtxt("Standard_devD", skiprows=1)
    Stand200 = StandDevD[-1]
    StandDevE = np.loadtxt("Standard_devE", skiprows=1)
    Stand250 = StandDevE[-1]
    StandDevF = np.loadtxt("Standard_devF", skiprows=1)
    Stand300 = StandDevF[-1]

    y = [SSR10, SSR50, SSR100, SSR200, SSR250, SSR300]
    err = [Stand10, Stand50, Stand100, Stand200, Stand250, Stand300]
    #plt.errorbar(x_vals, y, yerr=err, marker = 'o')
    #plt.xlim(5000, 105000)
    #axarr[1].title('Plot of SSR for each Number of Runs')


    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].errorbar(x_vals, y_vals, yerr=yerr, marker = 'o')
    axarr[1].errorbar(x_vals, y, yerr=err, marker = 'o')
    axarr[0].set_title('Varying trial moves for Bad Solvent')
    #axarr[1].set_title('Plot of SSR for each Number of Runs')
    axarr[1].set_xlim(5000, 305000)
    axarr[0].set_ylabel("Energy")
    axarr[1].set_ylabel("SSR")
    axarr[0].set_xlabel("Number of moves")
    axarr[1].set_xlabel("Number of moves")
    plt.show()

    binomial_graph("Saved_spectra300000")





