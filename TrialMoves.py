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

    #monomer = rnd.choice(chain[chain[:,0]>=0])

    #find the zeroith monomer
    monomer_0 = chain[0]
    #choose a random monomer to start the regrowth from
    index = rnd.choice(range(0,length))

    #Look up moiety of the zeroith monomer
    moiety = lattice[monomer_0[0], monomer_0[1], monomer_0[2]]

    #while monomer[0] < 1 and monomer[1] < 1:
    #    monomer = rnd.choice(chain)
    #    print "new monomer"
    #chain = chain.tolist()


    #init_index = np.intersect1d(np.where((chain[:,0] == monomer[0]))[0], np.where((chain[:,1] == monomer[1]))[0])
    #index = np.intersect1d(init_index, np.where((chain[:,2] == monomer[2]))[0])


    #if chain[length-1, 0] < 0 and chain[length-1, 1] < 0:
    #    length = int(length*.33)+1

#A check to make sure there is only one index. Not needed if doing the random integer method.
    #if index.shape[0] > 1:
     #   print index
      #  print monomer
    if np.any(index) > length:
        print 'monomer is: ' + str(monomer) + 'value at index is: ' + str(chain[index,:])
        index = length - 1
        monomer = chain[index,:,:]
        print 'move index back'
    numchains = length - index #have to do index[0] if using array form
    lattice[chain[index:length,0].astype(int),chain[index:length,1].astype(int), chain[index:length,2].astype(int)] = 0
    deleted = init.hex_Saw(lattice,numchains,moiety,chain[index])
    if np.any(deleted == -1):
        print "no saw found"
        return (lattice_old, chains_old, energy_old,count)
    #print chain
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
    monomer = chain[0]
    moiety = lattice[monomer[0], monomer[1], monomer[2]]
    deleted = init.hex_Saw(lattice,chainlength,moiety, graft_point[0])
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


#We need to find the moiety of the chains so we can switch the moeities of the swaped chains.
    monomer1 = chain1[0]
    moiety1 = lattice[monomer1[0], monomer1[1], monomer1[2]]
    monomer2 = chain2[0]
    moiety2 = lattice[monomer2[0], monomer2[1], monomer2[2]]



    while length1 == length2:
        chain2 = rnd.choice(chains)
        length2 = chain2[chain2[:,0]>=0].shape[0]
    if length1 < length2:
        #We will be deleting some of chain2 and growing onto chain1
        #the moiety of chain2 = moiety1 and the moiety of chain1 = moiety2 and moiety of regrowth = moiety2
        #print "if"

        #Swap moieties
        for step1 in chain1:
            lattice[step1[0],step1[1],step1[2]] = moiety2
        for step2 in chain2:
            lattice[step2[0],step2[1],step2[2]] = moiety1


        #Part without zeros = acceptable
        acceptable = chain1[chain1[:,0]>=0]
        #finding the # in acceptable
        index = acceptable.shape[0]
        #numchain = thee number of monomer units to regrow on the short chain
        numchains = length2 - index
        #changing the value of the excess units in the long chain to 0. deleting them
        lattice[chain2[index-1:length2,0].astype(int),chain2[index-1:length2,1].astype(int), chain2[index-1:length2,2].astype(int)] = 0
        #Regrowing monomers onto short chain
        deleted = init.hex_Saw(lattice,numchains+1,moiety2,chain1[index-1,:])
        if np.any(deleted == -1):
            print "no saw found"
            return (lattice_old, chains_old, energy_old,count)
        #adding regrowth to the chain
        chain1[index-1:length2,:] = deleted
        #deleting in long chain
        chain2[index:length2,:] = -1

    else:
        #print "else"
        #lenght2<length1
        #
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