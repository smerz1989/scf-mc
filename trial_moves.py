#!/usr/bin/python
import numpy as np
import random as rnd
import math
import initiation as init
import analysis as an

# Configurational-Bias Monte Carlo
def cbmc(lattice, chains):
    count = 0
    chains_old = np.copy(chains)
    lattice_old = np.copy(lattice)
    energy_old = an.calcenergy(chains, lattice)
    chain = rnd.choice(chains)
    length = chain[chain[:,0]>=0].shape[0]

    monomer_0 = map(int,chain[0])

    # Randomly choose a monomer to start the regrowth from
    index = rnd.choice(range(0,length))

    # Look up the moiety of the zeroth monomer
    moiety = int(lattice[monomer_0[0], monomer_0[1], monomer_0[2]])

    if np.any(index) > length:
        print 'monomer is: ' + str(monomer) + 'value at index is: ' + str(chain[index,:])
        index = length - 1
        monomer = chain[index,:,:]
        print 'move index back'
    numchains = length - index
    lattice[chain[index:length,0].astype(int),chain[index:length,1].astype(int), chain[index:length,2].astype(int)] = 0
    deleted = init.hex_Saw(lattice,numchains,moiety,chain[index])
    if np.any(deleted == -1):
        print "no saw found"
        return (lattice_old, chains_old, energy_old,count)
    print chain[2],
    print "chain 2"
    chain[index:length,:] = [0,0,0]
    chain[index:length,:] = deleted
    length_new = chain[chain[:,0]>=0].shape[0]
    energy_new = an.calcenergy(chains, lattice)
    if length != length_new:
        print 'Oringinal length: ' + str(length) + ' The new length is: ' + str(length_new) + '\n there is a difference of: ' + str(length - length_new)
        print '\n The index is: ' + str(index) + ' and monomer is: ' + str(monomer)
    if energy_new < energy_old :
        count += 1
        print "accept"
        print index
        print numchains
        return (lattice, chains, energy_new, count)
    elif rnd.uniform(0,1) < math.exp(-(energy_new-energy_old)):
        count += 1
        print "accept"
        print index
        print numchains
        return (lattice, chains, energy_new, count)
    else:
        print index
        print numchains
        print "reject"
        return (lattice_old, chains_old, energy_old, count)

def take_empty(lattice, chains, graft_points):
    count = 0
    numchains = chains.shape[0]
    chains_old = np.copy(chains)
    lattice_old = np.copy(lattice)
    energy_old = an.calcenergy(chains, lattice)
    num = rnd.randint(0,numchains-1)
    chain = chains[num,:,:]
    monomer = chain[0]
    moiety = int(lattice[monomer[0], monomer[1], monomer[2]])

    chainlength = chain[chain[:,0]>=0].shape[0]

    index = num

    print index
    print chainlength

    chain[0:chainlength,:] = [0,0,0]

    graft_point = rnd.sample(graft_points,1)

    print graft_point[0]

    deleted = init.hex_Saw(lattice,chainlength,moiety, graft_point[0])
    if np.any(deleted == -1):
        print "no saw found"
        return (lattice_old, chains_old, energy_old,count)
    chains[num,0:chainlength,:] = deleted

    length_new = chain[chain[:,0]>=0].shape[0]
    energy_new = an.calcenergy(chains, lattice)
    if chainlength != length_new:
        print 'Oringinal length: ' + str(chainlength) + ' The new length is: ' + str(length_new) + '\n there is a difference of: ' + str(chainlength - length_new)
        print '\n The index is: ' + str(index)

    if energy_new < energy_old :
        count += 1
        print "accept"
        print index
        print numchains
        return (lattice, chains, energy_new, count)
    elif rnd.uniform(0,1) < math.exp(-(energy_new-energy_old)):
        count += 1
        print "accept"
        print index
        print numchains
        return (lattice, chains, energy_new, count)
    else:
        print index
        print numchains
        print "reject"
        return (lattice_old, chains_old, energy_old, count)

def swap(lattice, chains):
    count = 0
    chains_old = np.copy(chains)
    lattice_old = np.copy(lattice)
    energy_old = an.calcenergy(chains, lattice)
    chain1 = rnd.choice(chains)
    length1 = chain1[chain1[:,0]>=0].shape[0]
    chain2 = rnd.choice(chains)
    length2 = chain2[chain2[:,0]>=0].shape[0]

    # Find the moiety of the chains so we can switch the moieties of the swapped chains
    monomer1 = chain1[0]
    moiety1 = lattice[int(monomer1[0]), int(monomer1[1]), int(monomer1[2])]
    monomer2 = map(int,chain2[0])
    moiety2 = lattice[monomer2[0], monomer2[1], monomer2[2]]
    while length1 == length2:
        chain2 = rnd.choice(chains)
        length2 = chain2[chain2[:,0]>=0].shape[0]
    if length1 < length2:
        #We will be deleting some of chain2 and growing onto chain1
        #the moiety of chain2 = moiety1 and the moiety of chain1 = moiety2 and moiety of regrowth = moiety2

        # Swap the moieties
        for step1 in chain1:
            lattice[step1[0],step1[1],step1[2]] = moiety2
        for step2 in chain2:
            lattice[step2[0],step2[1],step2[2]] = moiety1

        # Find the acceptable part, or the part without zeros
        acceptable = chain1[chain1[:,0]>=0]
        # Find the number of acceptable results in acceptable
        index = acceptable.shape[0]
        # numchain is the number of monomer units to regrow on the short chain
        numchains = length2 - index
        # Change the value of the excess units in the long chain to 0 (essentially deleting them)
        lattice[chain2[index-1:length2,0].astype(int),chain2[index-1:length2,1].astype(int), chain2[index-1:length2,2].astype(int)] = 0
        # Regrow monomers onto short chain
        deleted = init.hex_Saw(lattice,numchains+1,moiety2,chain1[index-1,:])
        if np.any(deleted == -1):
            print "no saw found"
            return (lattice_old, chains_old, energy_old,count)
        # Add regrowth to the chain
        chain1[index-1:length2,:] = deleted
        # Delete monomers in long chain
        chain2[index:length2,:] = -1

    else:
        acceptable = chain2[chain2[:,0]>=0]
        index = acceptable.shape[0]
        numchains = length1 - index
        lattice[chain1[index-1:length1,0].astype(int),chain1[index-1:length1,1].astype(int), chain1[index-1:length1,2].astype(int)] = 0
        deleted = init.hex_Saw(lattice,numchains+1,moiety1,chain2[index-1,:])
        if np.any(deleted == -1):
            print "no saw found"
            return (lattice_old, chains_old, energy_old,count)
        chain2[index-1:length1,:] = deleted
        chain1[index:length1,:] = -1
    energy_new = an.calcenergy(chains, lattice)
    if energy_new < energy_old :
        count += 1
        print "accept"
        print index
        print numchains
        return (lattice, chains, energy_new, count)
    elif rnd.uniform(0,1) < math.exp(-(energy_new-energy_old)):
        count += 1
        print "accept"
        print index
        print numchains
        return (lattice, chains, energy_new, count)
    else:
        print index
        print numchains
        print "reject"
        return (lattice_old, chains_old, energy_old, count)