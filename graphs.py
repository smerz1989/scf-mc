
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
        print [last_line("Energies-10000" + str(c))]
        energy_10 += [last_line("Energies-10000" + str(c))]
    energy_50 = []
    for d in "abcd":
        energy_50 += [last_line("Energies-50000" + str(d))]
    energy_100 = []
    for e in "abcd":
        energy_100 += [last_line("Energies-100000" + str(e))]
    energy_200 = []
    for f in "abc":
        energy_200 += [last_line("Energies-200000" + str(f))]


    avg10 = sum(energy_10)/len(energy_10)
    avg50 = sum(energy_50)/len(energy_50)
    avg100 = sum(energy_100)/len(energy_100)
    avg200 = sum(energy_200)/len(energy_200)

    std_10 = np.std(energy_10)
    std_50 = np.std(energy_50)
    std_100 = np.std(energy_100)
    std_200 = np.std(energy_200)

    y_vals = [avg10, avg50, avg100, avg200]
    x_vals = [10000, 50000, 100000, 200000]
    yerr = [std_10, std_50, std_100, std_200]
    #plt.errorbar(x_vals, y_vals, yerr=yerr, marker = 'o')
    #plt.xlim(5000, 105000)
    #axarr[0].title('Plot of Energies for each Number of Runs')
    #plt.scatter(x_vals, energy_50, label='Energy for 50,000')
    #plt.scatter(x_vals, energy_100, label='Energy for 100,000')

    #legend = plt.legend(loc='upper right', shadow = True, fontsize = 'x-large')
   # plt.figure(1)
    plt.show()
   # plt.figure(2)
    SSR10 = get_avg("Saved_SSR-e")
    SSR50 = get_avg("Saved_SSR-f")
    SSR100 = get_avg("Saved_SSR-g")
    SSR200 = get_avg("Saved_SSR-h")

    StandDevA = np.loadtxt("Standard_dev-A", skiprows=2)
    Stand10 = StandDevA[-1]
    StandDevB = np.loadtxt("Standard_dev-B", skiprows=1)
    Stand50 = StandDevB[-1]
    StandDevC = np.loadtxt("Standard_dev-C", skiprows=2)
    Stand100 = StandDevC[-1]
    StandDevD = np.loadtxt("Standard_dev-D", skiprows=1)
    Stand200 = StandDevD[-1]

    y = [SSR10, SSR50, SSR100, SSR200]
    err = [Stand10, Stand50, Stand100, Stand200]
    #plt.errorbar(x_vals, y, yerr=err, marker = 'o')
    #plt.xlim(5000, 105000)
    #axarr[1].title('Plot of SSR for each Number of Runs')


    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].errorbar(x_vals, y_vals, yerr=yerr, marker = 'o')
    axarr[1].errorbar(x_vals, y, yerr=err, marker = 'o')
    axarr[0].set_title('Determining optimum number of runs based on energy and SSR')
    #axarr[1].set_title('Plot of SSR for each Number of Runs')
    axarr[1].set_xlim(5000, 205000)
    axarr[0].set_ylabel("Energy")
    axarr[1].set_ylabel("SSR")
    axarr[0].set_xlabel("Number of Runs")
    axarr[1].set_xlabel("Number of Runs")
    plt.show()

    binomial_graph("Saved_spectra10000a")





