
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
    for k in xrange(0,6):
        n =5
        binomial += [choose(n,k)*(0.5**k)*(0.5**(n-k))]

    t = np.array([0,1,2,3,4,5])
    plt.plot(t.T,binomial,'r')
    print yvals
    for i in range(yvals.shape[0]):
        if i == 1:
            plt.plot(t.T,yvals[i],'b')
        elif i == 2:
            plt.plot(t.T,yvals[i],'g')
        elif i == 3:
            plt.plot(t.T,yvals[i],'c')
        elif i == 4:
            plt.plot(t.T,yvals[i],'y')
        elif i == 5:
            plt.plot(t.T,yvals[i],'r')
        else:
            plt.plot(t.T, yvals[i], 'm')

    plt.title(file)
    plt.show()






if __name__ == '__main__':
    energy_1 = []
    for c in "abc":
        #print [last_line("Energies-10000" + str(c))]
        energy_1 += [last_line("Energies3D_1050" + str(c))]
    energy_5 = []
    for d in "abc":
        energy_5 += [last_line("Energies3D_5050" + str(d))]
    energy_10 = []
    for e in "abc":
        energy_10 += [last_line("Energies3D_10050" + str(e))]
    energy_50 = []
    for f in "abc":
        energy_50 += [last_line("Energies3D_50000" + str(f))]
    energy_90 = []
    for g in "abc":
        energy_90 += [last_line("Energies3D_90000" + str(g))]
    energy_70 = []
    for h in "abc":
        energy_70 += [last_line("Energies3D_70000" + str(h))]


    avg1 = sum(energy_1)/len(energy_1)
    avg5 = sum(energy_5)/len(energy_5)
    avg10 = sum(energy_10)/len(energy_10)
    avg50 = sum(energy_50)/len(energy_50)
    avg90 = sum(energy_90)/len(energy_90)
    avg70 = sum(energy_70)/len(energy_70)

    std_1 = np.std(energy_1)
    std_5 = np.std(energy_5)
    std_10 = np.std(energy_10)
    std_50 = np.std(energy_50)
    std_90 = np.std(energy_90)
    std_70 = np.std(energy_70)

    y_vals = [avg1, avg5, avg10, avg50,avg70, avg90] # avg100, avg200, avg250, avg300]
    x_vals = [1050, 5050, 10050, 50000, 70000, 90000]# 100000, 200000, 250000, 300000]
    yerr = [std_1, std_5, std_10, std_50,std_70, std_90]# std_100, std_200, std_250, std_300]
    #plt.errorbar(x_vals, y_vals, yerr=yerr, marker = 'o')
    #plt.xlim(5000, 105000)
    #axarr[0].title('Plot of Energies for each Number of Runs')
    #plt.scatter(x_vals, energy_50, label='Energy for 50,000')
    #plt.scatter(x_vals, energy_100, label='Energy for 100,000')

    #legend = plt.legend(loc='upper right', shadow = True, fontsize = 'x-large')
   # plt.figure(1)
    plt.show()
    """ # plt.figure(2)
    SSR10 = get_avg("Saved_SSRe")
    SSR50 = get_avg("Saved_SSRf")"""
    SSR10 = get_avg("Saved_SSR3D_10050")
    SSR50 = get_avg("Saved_SSR3D_50000")
    SSR90 = get_avg("Saved_SSR3D_90000")
    SSR70 = get_avg("Saved_SSR3D_70000")

    StandDevA = np.loadtxt("Standard_dev3D_10050")
    Stand10 = StandDevA[-1]
    StandDevB = np.loadtxt("Standard_dev3D_50000")
    Stand50 = StandDevB[-1]
    StandDevC = np.loadtxt("Standard_dev3D_90000_b")
    Stand90 = StandDevC[-1]
    StandDevD = np.loadtxt("Standard_dev3D_70000_b")
    Stand70 = StandDevD
    """StandDevE = np.loadtxt("Standard_devE", skiprows=1)
    Stand250 = StandDevE[-1]
    StandDevF = np.loadtxt("Standard_devF", skiprows=1)
    Stand300 = StandDevF[-1]"""

    x = [10050, 50000,70000, 90000]
    y = [SSR10, SSR50,SSR70, SSR90]
    err = [Stand10, Stand50, Stand70, Stand90]
    #plt.errorbar(x, y, yerr=err, marker = 'o')
    #plt.xlim(0, 91000)
    #plt.title('Plot of SSR for each Number of Runs')
    #plt.xlabel('Number of Moves')



    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].errorbar(x_vals, y_vals, yerr=yerr, marker = 'o')
#    axarr[1].errorbar(x, y, yerr=err, marker = 'o')
    axarr[0].set_title('Varying trial moves for Good Solvent')
    axarr[1].set_title('Plot of SSR for each Number of Runs')
    axarr[1].set_xlim(0, 91000)
    axarr[0].set_ylabel("Energy")
    axarr[1].set_ylabel("SSR")
    axarr[0].set_xlabel("Number of moves")
    axarr[1].set_xlabel("Number of moves")
    plt.show()

    binomial_graph("Saved_spectra3D_90000")





