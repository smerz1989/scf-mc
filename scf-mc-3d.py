#!/usr/bin/python
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import math
import sys
import collections
import itertools
from mpl_toolkits.mplot3d import Axes3D

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

def calcenergy(chains, lattice):
    (L, M, P) = lattice.shape
    nns_odd = np.array([(-1, -1, 0), (0, -1, 0), (1, -1, 0), (-1, 0, 0), (1, 0, 0), (0, 1, 0), (0,0,-1),(0,0,1)])
    nns_even = np.array([(0, -1, 0), (-1, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (-1, 1, 0),(0,0,-1),(0,0,1)])
    total_energy = 0
    chiAS = -1
    for chain in chains:
        #print chain
        length = chain.shape[0]
        chainlength = chain[chain[:,0]>=0].shape[0]
        if chainlength < length:
            for monomer in chain[0:chainlength,:]:
               # print monomer
                if monomer[0]%2 == 0:
                    nns = nns_even
                else:
                    nns = nns_odd
                neighbors = [[(neighbor[0]+monomer[0])%(L-1) ,(neighbor[1]+monomer[1])%(M-1) , (neighbor[2]+monomer[2])%(P-1) ] for neighbor in nns]
                for neighbor in neighbors:
                    if lattice[(neighbor[0]),(neighbor[1]),(neighbor[2])] == 0:
                        total_energy += chiAS/2.0
        else:
            for monomer in chain:
                if monomer[0]%2 == 0:
                    nns = nns_even
                else:
                    nns = nns_odd
                neighbors = [[(neighbor[0]+monomer[0])%(L-1) ,(neighbor[1]+monomer[1])%(M-1) , (neighbor[2]+monomer[2])%(P-1) ] for neighbor in nns]
                for neighbor in neighbors:
                   # print neighbor[0]
                    #print str(L) + "  " + str(M) +"  "+str(P)
                    if lattice[(neighbor[0]),(neighbor[1]), (neighbor[2])] == 0:
                        total_energy += chiAS/2.0
    return total_energy


def hex_Saw(lattice, numsteps, grafted_to=None):
    (L, M, D) = lattice.shape
    # print "Lattice at beginning is "+str(lattice)
    nns_odd = np.array([(-1, -1, 0), (0, -1, 0), (1, -1, 0), (-1, 0, 0), (1, 0, 0), (0, 1, 0), (0,0,-1),(0,0,1)])
    nns_even = np.array([(0, -1, 0), (-1, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (-1, 1, 0),(0,0,-1),(0,0,1)])
    complete = False
    numrestarts = 0
    maxtries = 1
    tries = 0
    while (not complete) and (tries < maxtries):
        steps = np.zeros([numsteps, 3], dtype=int)
        if (grafted_to == None):
            (row, col, layer) = (np.random.randint(0, L - 1), np.random.randint(0, M - 1),np.random.randint(0, D - 1))
        else:
            (row, col, layer) = grafted_to
        if lattice[row,col,layer] == 1:
            return -1
        # print "Starting location is "+str(row)+","+str(col)
        steps[0, :] = [row, col,layer]
        lattice[row, col,layer] = 1
        for i in xrange(numsteps - 1):
            if (row % 2 == 0):
                currentNNs = np.array([[int(row + nn[0]) % L, int(col + nn[1]) % M,int(layer + nn[2]) % D] for nn in nns_even])  # periodic bc
            else:
                currentNNs = np.array([[int(row + nn[0]) % L, int(col + nn[1]) % M,int(layer + nn[2]) % D] for nn in nns_odd])

            possibleMoves = currentNNs[np.where((lattice[currentNNs[:, 0], currentNNs[:, 1], currentNNs[:, 2]] != 1))[0], :]
            #print "Possible moves are "+str(possibleMoves)
            if (possibleMoves.shape[0] > 0):
                move = rnd.choice(possibleMoves)
                complete = True
            else:
                numrestarts += 1
                lattice[steps[0:(i + 1), 0], steps[0:(i + 1), 1],steps[0:(i + 1), 2]] = 0
                complete = False
                break
            # print "Chosen move is "+str(move)
            steps[i + 1, :] = move
            lattice[move[0], move[1], move[2]] = 1
            (row, col, layer) = (move[0], move[1], move[2])
        tries += 1;
    if (not complete):
        return -1
    print str(numsteps) + " long Chain grown with " + str(numrestarts) + " restarts"
    return steps


def graft_chains(lattice,graft_points,chainlength,numchains):
    chosen_points = rnd.sample(graft_points,numchains) if numchains<=graft_points.shape[0] else -1
    chains = np.empty((numchains,chainlength,3))
    idx = 0
   # print chosen_points[0]
    for graft_point in (chosen_points):
        chain = -1
        point = graft_point
        while np.any(chain==-1):
            chain = hex_Saw(lattice,chainlength,point)
            if np.any(chain == -1):
                print 'moving chain to new graft point'
                point = rnd.choice(graft_points)
        chains[idx,:,:] = chain
        idx += 1
        #print idx
        #print graft_point
    return chains

def fill_cylinder(lattice,radius):
	center = [lattice.shape[0]/2,lattice.shape[1]/2,lattice.shape[2]/2]
	startoffset=0
	endoffset = 0
	num_graft_points = 4*radius+2*radius
	graft_points = np.empty((num_graft_points,3))
	for i in xrange(radius):
		if(((i+center[0])%2)==0):
			startoffset+=1
		else:
			endoffset+=1
		lattice[center[0]+i,(center[1]-radius+startoffset):(center[1]+radius-endoffset),:]=1
		lattice[center[0]-i,(center[1]-radius+startoffset):(center[1]+radius-endoffset),:]=1
		lattice[(center[0]-i),[(center[1]-radius+startoffset-1),(center[1]+radius-endoffset)],:]=-1
		lattice[(center[0]+i),[(center[1]-radius+startoffset-1),(center[1]+radius-endoffset)],:]=-1
	lattice[center[0]+radius,(center[1]-radius+startoffset):(center[1]+radius-endoffset),:]=-1
	lattice[center[0]-radius,(center[1]-radius+startoffset):(center[1]+radius-endoffset),:]=-1
	graft_points = np.transpose(np.vstack(np.where((lattice==-1))))
	print "Graft points look like this "+str(graft_points)
	lattice[graft_points[:,0],graft_points[:,1],graft_points[:,2]]=0
	return graft_points

	
def initialize_lattice(shape,numchains,chainlength,radius):
    chainsize_A = int(chainlength*.75)
    chainsize_B = chainlength - chainsize_A
    numchainsA = int(numchains/2)
    numchainsB = numchains - numchainsA

    lattice = np.zeros(shape)
    #radius = int(0.25*lattice.shape[1])
    graft_points = fill_cylinder(lattice,radius)
    #chains = graft_chains(lattice,graft_points,chainlength,numchains)
    chains_A = graft_chains(lattice, graft_points, chainsize_A, numchainsA)
    chains_Bold = graft_chains(lattice, graft_points, chainsize_B, numchainsB)
    chains_B = np.empty([numchainsB, chainsize_A, 3])
    for i in range(numchainsB):
        addition = np.full((chainsize_A-chainsize_B,3),-1, dtype=np.int)
        chains_B[i] = np.concatenate((chains_Bold[i],addition),0)
    chains = np.concatenate((chains_A,chains_B), 0)
    return (lattice,graft_points,chains)

def cbmc(lattice, chains):  #configurational bias monte carlo
    count = 0
    chains_old = np.copy(chains)
    lattice_old = np.copy(lattice)
    energy_old = calcenergy(chains, lattice)
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
    deleted = hex_Saw(lattice,numchains,monomer)
    if np.any(deleted == -1):
        print "no saw found"
        return (lattice_old, chains_old, energy_old,count)
    #print chain
    print chain[2],
    print "chain 2"
    chain[index[0]:length,:] = [0,0,0]
    chain[index[0]:length,:] = deleted
    length_new = chain[chain[:,0]>=0].shape[0]
    energy_new = calcenergy(chains, lattice)
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
    energy_old = calcenergy(chains, lattice)
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
    deleted = hex_Saw(lattice,chainlength,graft_point[0])
    if np.any(deleted == -1):
        print "no saw found"
        return (lattice_old, chains_old, energy_old,count)
    chains[num,0:chainlength,:] = deleted
    #print chains[num,:,:]
    length_new = chain[chain[:,0]>=0].shape[0]
    energy_new = calcenergy(chains, lattice)
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
    energy_old = calcenergy(chains, lattice)
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
        deleted = hex_Saw(lattice,numchains+1,chain1[index-1,:])
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
        deleted = hex_Saw(lattice,numchains+1,chain2[index-1,:])
        if np.any(deleted == -1):
            print "no saw found"
            return (lattice_old, chains_old, energy_old,count)
        chain2[index-1:length1,:] = deleted
        chain1[index:length1,:] = -1
    energy_new = calcenergy(chains, lattice)
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


def chains_to_xyz(chains,filename):
    chainlength_A = chains[0].shape[0]
    chainlength_B = int(chainlength_A*.33)+1
    chainlength_A = 12
    chainlength_B = 4
    numchains = chains.shape[0]

    xyzfile = open(filename ,'a')
    numatoms = chainlength_A*numchains/2 + chainlength_B*numchains/2;
   # numatoms = numchains*30
    xyzfile.write(str(numatoms)+'\n\n')
    for chain in chains:
        chain = chain[chain[:,0] >= 0]
       # xyzfile.write('# '+str(chain))
        for monomer in chain:
            distance = math.sqrt((monomer[0]-chain[0][0])**2 + (monomer[1] - chain[0][1])**2 + (monomer[2]-chain[0][2])**2)
            #if distance > 13:
            #    print 'distance between tether is ' + str(distance)
            if chain.shape[0] == chainlength_A:
                xyzfile.write('C '+str(monomer[0])+' '+str(monomer[1])+' '+ str(monomer[2])+'\n')
            elif chain.shape[0] == chainlength_B:
                xyzfile.write('N ' + str(monomer[0])+' '+str(monomer[1])+' '+ str(monomer[2])+'\n')
            else:
                xyzfile.write('O '+str(monomer[0])+' '+str(monomer[1])+' '+ '0'+'\n')
    xyzfile.close()

def store_energies(energy):
    energyfile = open('Energies3D_'+str(n)+alph,'a')
    energyfile.write(str(energy)+"\n")

def acceptance_rate(i, count):
    rate = count/float(i)
    print "\n\n The rate is: " + str(rate*100) + "%"
    return rate

def sep_analysis(chains):
    sort = np.array([0,0,0])
    zero_b = 0
    one_b = 0
    two_b = 0
    three_b = 0
    four_b = 0
    five_b = 0
    for chain_sel in chains:
        nn_list = [(chain, (math.sqrt((chain_sel[0][0]-chain[0][0])**2 + (chain_sel[0][1]-chain[0][1])**2 + (chain_sel[0][2]-chain[0][2])**2)),(-1 if len(chain[chain[:,0]>=0]) == 4 else 1 )) \
         for chain in chains if ((math.sqrt((chain_sel[0][0]-chain[0][0])**2 + (chain_sel[0][1]-chain[0][1])**2 + (chain_sel[0][2]-chain[0][2])**2)) <2)]
        #print len(nn_list)
        type = []
        for i in xrange(len(nn_list)-1):
            count = 0
            type += [nn_list[i][2]]
        for item in type:
            if item == 1:
                count += 1
        if count == 0:
            zero_b += 1
        elif count == 1:
            one_b += 1
        elif count == 2:
            two_b += 1
        elif count == 3:
            three_b += 1
        elif count == 4:
            four_b += 1
        elif count == 5:
            five_b += 1
    return(zero_b, one_b, two_b, three_b, four_b, five_b)


def SSR(analysis, binomial):
    ssr = 0
    for i in range(0,6):
        diff = analysis[i]-binomial[i]
        ssr += diff**2
    return ssr


if __name__=='__main__':
    (lattice,graft_points,chains) = initialize_lattice((200,200,200),7000,16,20)
    print len(graft_points)

    count = 0
    n = sys.argv[1]
    alph = sys.argv[2]
    chains_to_xyz(chains, 'Init3D_'+str(n)+alph)

    for i in range(0,10):
        rando = rnd.uniform(0,1)
        if rando>.66:
            (lattice, chains,total_energy, acc)=cbmc(lattice,chains)
            print "CBMC"
        elif rando > .33:
            (lattice, chains, total_energy, acc) = take_empty(lattice, chains, graft_points)
            print "Empty"
        else:
            (lattice, chains, total_energy, acc) = swap(lattice, chains)
            print "SWAP"
        store_energies(total_energy)
        chains_to_xyz(chains, 'Long3D_'+str(n)+alph)
        count += acc
        acceptance_rate(i+1,count)

    chains_to_xyz(chains, 'Short3D_'+str(n)+alph)

    analysis = sep_analysis(chains)
    analysis = tuple(x/float(2400) for x in analysis)
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
        n =5
        binomial += [choose(n,k)*(0.5**k)*(0.5**(n-k))]

    spectra = open('Saved_spectra3D_'+str(n),'a')
    spectra.write(str(analysis)+"\n")

    ssr = SSR(analysis, binomial)
    SSR = open('Saved_SSR3D_'+str(n), 'a')
    SSR.write('+1\t' + str(ssr) + '\t' + '50\n')
    SSR.close()

    S = np.loadtxt('Saved_SSR3D_'+str(n))

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