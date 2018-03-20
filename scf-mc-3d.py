#!/usr/bin/python
import analysis as an
import initiation as init
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import random as rnd
import sys
import trial_moves as tm

# NOTES
# lattice length = 0.154 nm
# moieties: solvent = 0, substrate = 1, polymer A = 2, polymer B = 3
# x and y components of shape has to be at least 2*radius + 2*chainlength

if __name__== '__main__':

    # Initializes the lattice and chains, and defines and prints the number of graft points
    (lattice, graft_points, chains) = init.initialize_lattice((200,200,20),708,16,[2,3],20)
    print len(graft_points)

    # Initializes the step counter and saves system variables for naming purposes
    count = 0
    x = int(sys.argv[1])
    y = sys.argv[2]

    # Creates a new directory "data" and saves the initial setup to a file
    if not os.path.exists(os.getcwd()+'/data'):
        os.mkdir(os.getcwd()+'/data')
    an.chains_to_xyz(chains, 'InitDual_'+str(x)+y, lattice)

    for i in range(0,x):
        rng = rnd.uniform(0,1)
        if rng > 0.66:
            (lattice, chains, total_energy, acc) = tm.cbmc(lattice, chains)
            print "CBMC"
        elif rng > 0.33:
            (lattice, chains, total_energy, acc) = tm.take_empty(lattice, chains, graft_points)
            print "empty"
        else:
            (lattice,chains, total_energy, acc) = tm.swap(lattice, chains)
            print "SWAP"
        an.store_energies(total_energy, x, y)

        # Record the chains every 50,000 steps so that the SSR can be taken later
        if i % 50000 == 0:
            np.save(os.getcwd()+'/data/chainsSSRDual'+str(i+1)+y, chains)

        # Record the chains every 100 configurations
        if i % 100 == 0:
            an.chains_to_xyz(chains, 'LongDual_'+str(x)+y, lattice)

        count += acc

        an.acceptance_rate(i+1,count)

    an.chains_to_xyz(chains, 'ShortDual_'+str(x)+y, lattice)

    analysis = an.sep_analysis(chains)
    analysis = tuple(x/float(sum(analysis)) for x in analysis)
    print analysis
    chainsSSR = np.load(os.getcwd()+'/data/chainsSSRDual'+str(x)+y+'.npy')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for chain in chains:
        chain = chain[chain[:,0] >= 0]
        ax.plot(chain[:,0],chain[:,1],chain[:,2])

    plt.plotfile('EnergiesDual_'+str(x)+y)
    plt.savefig(os.getcwd()+'/data/energyDual_' +str(x)+y+".pdf")
    plt.close()

    binomial = []
    for k in xrange(0,6):
        m =5
        binomial += [an.choose(m,k)*(0.5**k)*(0.5**(m-k))]

    saved = np.save(os.getcwd()+'/data/Safe_SSRDual'+str(x)+y, analysis)

    #spectra = open(r'C:\Users\Maggie\Documents\GitHub\scf-mc\Saved_spectra3D_'+str(n), 'w')
    #with spectra:
    spectra = open(os.getcwd()+'/data/Saved_spectraDual_'+str(x), 'a')
    spectra.write(str(analysis) +"\n")
    spectra.flush()
    spectra.close()

    ssr = an.SSR(analysis, binomial)
    SSR = open(os.getcwd()+'/data/Saved_SSRDual_'+str(x), 'a')
    SSR.write('-1\t' + str(ssr) + '\t' + '10000\n')
    SSR.flush()
    SSR.close()

    S = np.loadtxt(os.getcwd()+'/data/Saved_SSRDual_'+str(x))

    if y != 'a':
        set = []
        for line in S:
            set += [line[1]]
        set = map(float, set)
        std = open(os.getcwd()+'/data/Standard_devDual_'+str(x),'a')
        std.write(str(np.std(set))+'\n')

    print binomial

    t = np.array([0,1,2,3,4,5])
    plt.plot(t.T,binomial,'r')
    plt.plot(t.T,analysis,'b')
    plt.savefig(os.getcwd()+'/data/binomialDual'+str(x)+y+".pdf")
