#!/usr/bin/python
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import sys
import Initiation as init
import TrialMoves as TM
import Analysis as an



if __name__=='__main__':
    (lattice,graft_points,chains) = init.initialize_lattice((200,200,20),779,16,20)
    print len(graft_points)

    count = 0
    #n = sys.argv[1]
    #alph = sys.argv[2]
    n = 7
    alph = 'a'
    an.chains_to_xyz(chains, 'Init3D_'+str(n)+alph)

    for i in range(0,int(n)):
        rando = rnd.uniform(0,1)
        if rando>.66:
            (lattice, chains,total_energy, acc) = TM.cbmc(lattice,chains)
            print "CBMC"
        elif rando > .33:
            (lattice, chains, total_energy, acc) = TM.take_empty(lattice, chains, graft_points)
            print "Empty"
        else:
            (lattice, chains, total_energy, acc) = TM.swap(lattice, chains)
            print "SWAP"
        an.store_energies(total_energy, n, alph)

        # Should add matrices for chemical moeity and identity in here

        if i % 100 == 0: # record every hundredth configuration
            an.chains_to_xyz(chains, 'Long3D_'+str(n)+alph)
        count += acc
        an.acceptance_rate(i+1,count)

    an.chains_to_xyz(chains, 'Short3D_'+str(n)+alph)

    analysis = an.sep_analysis(chains)
    analysis = tuple(x/float(sum(analysis)) for x in analysis)
    print analysis

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for chain in chains:
        chain = chain[chain[:,0] >= 0]
        ax.plot(chain[:,0],chain[:,1],chain[:,2])
    plt.show()

    plt.plotfile('Energies3D_'+str(n)+alph)
    plt.show()

    binomial = []
    for k in xrange(0,6):
        m =5
        binomial += [an.choose(m,k)*(0.5**k)*(0.5**(m-k))]

    saved = np.save('Safe_SSR'+str(n)+alph, analysis)

    #spectra = open(r'C:\Users\Maggie\Documents\GitHub\scf-mc\Saved_spectra3D_'+str(n), 'w')
    #with spectra:
    spectra = open('Saved_spectra3D_'+str(n), 'a')
    spectra.write(str(analysis) +"\n")
    spectra.flush()
    spectra.close()

    ssr = an.SSR(analysis, binomial)
    SSR = open('Saved_SSR3D_'+str(n), 'a')
    SSR.write('-1\t' + str(ssr) + '\t' + '10000\n')
    SSR.flush()
    SSR.close()

    S = np.loadtxt('Saved_SSR3D_'+str(n))

    if alph != 'a':
        set = []
        for line in S:
            set += [line[1]]

        set = map(float, set)
        std = open('Standard_dev3D_'+str(n),'a')
        std.write(str(np.std(set))+'\n')

    print binomial

    t = np.array([0,1,2,3,4,5])
    plt.plot(t.T,binomial,'r')
    plt.plot(t.T,analysis,'b')
    plt.show()