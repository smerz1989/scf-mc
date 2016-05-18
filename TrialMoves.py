__author__ = 'Maggie'
#!/usr/bin/python
import numpy as np
import random as rnd
import math
import Initiation as init
import Analysis as an



def cbmc(lattice, chains):  #configurational bias monte carlo
    count = 0
    chains_old = np.copy(chains)
    lattice_old = np.copy(lattice)
    energy_old = an.calcenergy(chains, lattice)
    chain = rnd.choice(chains)
    length = chain[chain[:,0]>=0].shape[0]
   # print chain
    monomer = rnd.choice(chain[chain[:,0]>=0])
    #while monomer[0] < 1 and monomer[1] < 1:
    #    monomer = rnd.choice(chain)
    #    print "new monomer"
    #chain = chain.tolist()
    index = np.intersect1d(np.where((chain[:,0] == monomer[0]))[0], np.where((chain[:,1] == monomer[1]))[0])
    #if chain[length-1, 0] < 0 and chain[length-1, 1] < 0:
    #    length = int(length*.33)+1
    if index.shape[0] > 1:
        print index
        print monomer
    if np.any(index) > length:
        print 'monomer is: ' + str(monomer) + 'value at index is: ' + str(chain[index,:])
        index = length - 1
        monomer = chain[index,:,:]
        print 'move index back'
    numchains = length - index[0]
    lattice[chain[index[0]:length,0].astype(int),chain[index[0]:length,1].astype(int), chain[index[0]:length,2].astype(int)] = 0
    deleted = init.hex_Saw(lattice,numchains,monomer)
    if np.any(deleted == -1):
        print "no saw found"
        return (lattice_old, chains_old, energy_old,count)
    #print chain
    print chain[2],
    print "chain 2"
    chain[index[0]:length,:] = [0,0,0]
    chain[index[0]:length,:] = deleted
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
    #print chain
    chainlength = chain[chain[:,0]>=0].shape[0]
    #index = chain[0][:]
    index = num
    #print chain
    print index
    print chainlength
    #print chains[num,:,:]
    chain[0:chainlength,:] = [0,0,0]
    #print chains[num,:,:]
    graft_point = rnd.sample(graft_points,1)
#    while lattice[graft_point[0]][graft_point[1]][graft_point[2]]== 1:
#        graft_point = rnd.sample(graft_points,1)
    print graft_point[0]
    deleted = init.hex_Saw(lattice,chainlength,graft_point[0])
    if np.any(deleted == -1):
        print "no saw found"
        return (lattice_old, chains_old, energy_old,count)
    chains[num,0:chainlength,:] = deleted
    #print chains[num,:,:]
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
    while length1 == length2:
        chain2 = rnd.choice(chains)
        length2 = chain2[chain2[:,0]>=0].shape[0]
    if length1 < length2:
        #print "if"
        acceptable = chain1[chain1[:,0]>=0]
        index = acceptable.shape[0]
        numchains = length2 - index
        lattice[chain2[index-1:length2,0].astype(int),chain2[index-1:length2,1].astype(int), chain2[index-1:length2,2].astype(int)] = 0
        deleted = init.hex_Saw(lattice,numchains+1,chain1[index-1,:])
        if np.any(deleted == -1):
            print "no saw found"
            return (lattice_old, chains_old, energy_old,count)
        chain1[index-1:length2,:] = deleted
        chain2[index:length2,:] = -1

    else:
        #print "else"
        acceptable = chain2[chain2[:,0]>=0]
        index = acceptable.shape[0]
        numchains = length1 - index
        lattice[chain1[index-1:length1,0].astype(int),chain1[index-1:length1,1].astype(int), chain1[index-1:length1,2].astype(int)] = 0
        deleted = init.hex_Saw(lattice,numchains+1,chain2[index-1,:])
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