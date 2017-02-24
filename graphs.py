
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
    energy_5 = []
    for c in "abc":
        #print [last_line("Energies-10000" + str(c))]
        energy_5 += [last_line('chi_1_5000steps\Energies3D_5000' + str(c))]
    energy_10 = []
    for d in "bc":
        energy_10 += [last_line('chi_1_10000steps\EnergiesDual_10000' + str(d))]
    energy_50 = []
    for e in "abc":
        energy_50 += [last_line('chi_1_50000steps\EnergiesDual_50000' + str(e))]
    energy_100 = []
    for f in "abc":
        energy_100 += [last_line('chi_1_100000steps\EnergiesDual_100000' + str(f))]
    energy_150 = []
    for g in "abc":
        energy_150 += [last_line('chi_1_150000steps\EnergiesDual_150000' + str(g))]
    energy_200 = []
    for h in "abc":
        energy_200 += [last_line('chi_1_200000steps\EnergiesDual_200000' + str(h))]
    energy_250 = []
    for i in "abc":
        energy_250 += [last_line('chi_1_250000steps\EnergiesDual_250000' + str(i))]



    avg5 = sum(energy_5)/len(energy_5)
    avg10 = sum(energy_10)/len(energy_10)
    avg50 = sum(energy_50)/len(energy_50)
    avg100 = sum(energy_100)/len(energy_100)
    avg150 = sum(energy_150)/len(energy_150)
    avg200 = sum(energy_200)/len(energy_200)
    avg250 = sum(energy_250)/len(energy_250)

    std_5 = np.std(energy_5)
    std_10 = np.std(energy_10)
    std_50 = np.std(energy_50)
    std_100 = np.std(energy_100)
    std_150 = np.std(energy_150)
    std_200 = np.std(energy_200)
    std_250 = np.std(energy_250)

    y_vals = [avg5, avg10, avg50,avg100, avg150,avg200, avg250 ] # avg100, avg200, avg250, avg300]
    x_vals = [5000, 10000, 50000, 100000, 150000, 200000, 250000]# 100000, 200000, 250000, 300000]
    yerr = [std_5, std_10, std_50,std_100, std_150, std_200, std_250]# std_100, std_200, std_250, std_300]
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
    #SSR5 = get_avg("Saved_SSR3D_10050")
    SSR50 = get_avg("chi_1_50000steps\Saved_SSRDual_50000")
    SSR100 = get_avg("chi_1_100000steps\Saved_SSRDual_100000")
    SSR150 = get_avg("chi_1_150000steps\Saved_SSRDual_150000")
    SSR200 = get_avg("chi_1_200000steps\Saved_SSRDual_200000")
    SSR250 = get_avg("chi_1_250000steps\Saved_SSRDual_250000")

    StandDevA = np.loadtxt("chi_1_50000steps\Standard_devDual_50000")
    Stand50 = StandDevA[-1]
    StandDevB = np.loadtxt("chi_1_100000steps\Standard_devDual_100000")
    Stand100 = StandDevB[-1]
    StandDevC = np.loadtxt("chi_1_150000steps\Standard_devDual_150000")
    Stand150 = StandDevC
    StandDevD = np.loadtxt("chi_1_200000steps\Standard_devDual_200000")
    Stand200 = StandDevD[-1]
    StandDevE = np.loadtxt("chi_1_250000steps\Standard_devDual_250000")
    Stand250 = StandDevE
    #StandDevF = np.loadtxt("Standard_devF", skiprows=1)
    #Stand300 = StandDevF[-1]

    x = [50000, 100000, 150000, 200000, 250000]
    y = [SSR50,SSR100, SSR150, SSR200, SSR250]
    err = [Stand50, Stand100, Stand150, Stand200, Stand250]
    #plt.errorbar(x, y, yerr=err, marker = 'o')
    #plt.xlim(0, 91000)
    #plt.title('Plot of SSR for each Number of Runs')
    #plt.xlabel('Number of Moves')



    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].errorbar(x_vals, y_vals, yerr=yerr, marker = 'o')
    axarr[1].errorbar(x, y, yerr=err, marker = 'o')
    axarr[0].set_title('Varying trial moves for Good Solvent')
    axarr[1].set_title('Plot of SSR for each Number of Runs')
    axarr[1].set_xlim(0, 260000)
    axarr[0].set_ylabel("Energy")
    axarr[1].set_ylabel("SSR")
    axarr[0].set_xlabel("Number of moves")
    axarr[1].set_xlabel("Number of moves")
    plt.show()

    binomial_graph("chi_1_200000steps\Saved_spectraDual_200000")





