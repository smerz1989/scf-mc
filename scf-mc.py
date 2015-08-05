#!/usr/bin/python
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import math
import collections
import itertools


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



def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def calcenergy(chains, lattice):
    (L, M) = lattice.shape
    nns_odd = np.array([(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (0, 1)])
    nns_even = np.array([(0, -1), (-1, 0), (1, 0), (1, 1), (0, 1), (-1, 1)])
    total_energy = 0
    chiAS = 1
    for chain in chains:
        if chain[len(chain)-1,0]<1 and chain[len(chain)-1,1]<1:
            for monomer in chain[0:int(len(chain)*.333)+1,:]:
                if monomer[0]%2 == 0:
                    nns = nns_even
                else:
                    nns = nns_odd
                neighbors = [[neighbor[0]+monomer[0],neighbor[1]+monomer[1]] for neighbor in nns]
                for neighbor in neighbors:
                    if lattice[(neighbor[0]%L),(neighbor[1]%M-1)] == 0:
                        total_energy += chiAS/2.0
        else:
            for monomer in chain:
                if monomer[0]%2 == 0:
                    nns = nns_even
                else:
                    nns = nns_odd
                neighbors = [[neighbor[0]+monomer[0],neighbor[1]+monomer[1]] for neighbor in nns]
                for neighbor in neighbors:
                    if lattice[(neighbor[0]%L),(neighbor[1]%M-1)] == 0:
                        total_energy += chiAS/2.0
    return total_energy

def hex_Saw(lattice, numsteps, grafted_to=None):
    (L, M) = lattice.shape
    # print "Lattice at beginning is "+str(lattice)
    nns_odd = np.array([(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (0, 1)])
    nns_even = np.array([(0, -1), (-1, 0), (1, 0), (1, 1), (0, 1), (-1, 1)])
    complete = False
    numrestarts = 0
    maxtries = 50
    tries = 0
    while (not complete) and (tries < maxtries):
        steps = np.zeros([numsteps, 2], dtype=int)
        if (grafted_to == None):
            (row, col) = (np.random.randint(0, L - 1), np.random.randint(0, M - 1))
        else:
            (row, col) = grafted_to
        if lattice[row,col] == 1:
            return -1
        # print "Starting location is "+str(row)+","+str(col)
        steps[0, :] = [row, col]
        lattice[row, col] = 1
        for i in xrange(numsteps - 1):
            if (row % 2 == 0):
                currentNNs = np.array([[int(row + nn[0]) % L, int(col + nn[1]) % M] for nn in nns_even])  # periodic bc
            else:
                currentNNs = np.array([[int(row + nn[0]) % L, int(col + nn[1]) % M] for nn in nns_odd])
            # for neighbor in currentNNs:
            #	print str(neighbor)+" has lattice value "+str(lattice[neighbor[0],neighbor[1]])

            possibleMoves = currentNNs[np.where((lattice[currentNNs[:, 0], currentNNs[:, 1]] != 1))[0], :]
            #print "Possible moves are "+str(possibleMoves)
            if (possibleMoves.shape[0] > 0):
                move = rnd.choice(possibleMoves)
                complete = True
            else:
                # print "Reached dead end restarting chain"
                # print steps[0:(i+1),0]
                # print steps[0:(i+1),1]
                numrestarts += 1
                lattice[steps[0:(i + 1), 0], steps[0:(i + 1), 1]] = 0
                complete = False
                break
            # print "Chosen move is "+str(move)
            steps[i + 1, :] = move
            lattice[move[0], move[1]] = 1
            (row, col) = (move[0], move[1])

        tries += 1;
    if (not complete):
        return -1
    print str(numsteps) + " long Chain grown with " + str(numrestarts) + " restarts"
    return steps


def fill_nanoparticle(lattice, diameter):
    (L, M, D) = lattice.shape
    print L
    print M
    x0 = int(L/2)
    y0 = int(M/2)

  #  for x in range (0,L):
   #     for y in range (0,M):
    #     totalProb = 0
     #       if abs(y-y0)<=diameter and (abs(sqrt(3)/2*(x-x0) + 1/2*(y-y0))<=diameter) and abs(sqrt(3)/2*(x-x0) - 1/2*(y-y0)) <=diameter:
      #          lattice[x+1,y+1] = 1
    r = int(diameter/2)
    i = 0
    j = 0

    points = np.array([x0-r,y0+1])

    for y in range(1, r+1):
        lattice[x0-r:x0+r, y0] = 1
        points = np.vstack((points, [x0+r, y0+1]))
        x_val=[]
        y_val = []

        if y % 2 == 0 and y%r != 0:
            j += 1
            lattice[x0-r+i:x0+r-j, y0+i+j] = 1
            lattice[x0-r+i:x0+r-j, y0-i-j] = 1
            x_val = [x0-r+i-1, x0+r-j+1,x0-r+i-1,x0+r-j+1]
            y_val = [y0+i+j, y0+i+j, y0-i-j, y0-i-j]
            points = np.vstack((points, [x_val[0], y_val[0]]))
            points = np.vstack((points, [x_val[1], y_val[1]]))
            points = np.vstack((points, [x_val[2], y_val[2]]))
            points = np.vstack((points, [x_val[3], y_val[3]]))
        elif y % r == 0:
            if i == j:
                i += 1
            else:
                j += 1
            lattice[x0-r+i:x0+r-j, y0+i+j] = 1
            lattice[x0-r+i:x0+r-j, y0-i-j] = 1
            x_val += range(x0-r+i,x0+r-j) + range(x0-r+i, x0+r-j)
            num = int(len(range(x0-r+i,x0+r-j)))
            y_val += [y0+i+j+1] *num + [y0-i-j-1] *num
            for k in range(1,len(x_val)):
                points = np.vstack((points, [x_val[k], y_val[k]]))
        else:
            i += 1
            lattice[x0-r+i:x0+r-j, y0+i+j] = 1
            lattice[x0-r+i:x0+r-j, y0-i-j] = 1
            x_val = [x0-r+i-1, x0+r-j+1,x0-r+i-1,x0+r-j+1]
            y_val = [y0+i+j, y0+i+j, y0-i-j, y0-i-j]
            points = np.vstack((points, [x_val[0], y_val[0]]))
            points = np.vstack((points, [x_val[1], y_val[1]]))
            points = np.vstack((points, [x_val[2], y_val[2]]))
            points = np.vstack((points, [x_val[3], y_val[3]]))


    return (lattice, points)


def init_config(lattice_size, numchains, chainsize):
    chainsize_A = int(chainsize*.75)
    chainsize_B = chainsize - chainsize_A
    numchainsA = 29
    numchainsB = 19
    lattice = np.empty((lattice_size[0], lattice_size[1]))
    (np_lattice, points) = fill_nanoparticle(lattice, 40)
    (L,M) = (lattice_size[0],lattice_size[1])
    chains_A = np.empty([numchainsA, chainsize_A, 2])
    chains_B = np.empty([numchainsB, chainsize_A, 2])
    #plt.scatter(x_val, y_val)
    #plt.show
    graft_points = np.empty((numchains, 2))
    graft_step = int((int(3 * L / 4) - int(L / 4)) / numchains)
    #lattice[int(L / 4):int(3 * L / 4), int(M / 2)] = 1
    #graft_points[:, 0] = np.arange(start=int(L / 4), stop=int(3 * L / 4) - 5 * graft_step,
                                   #step=graft_step)
    #graft_points[:, 1] = int(lattice_size[1] / 2) + 1
    #possible_points = np.transpose(np.vstack((x_val, y_val)))

    points_convert = [math.atan2(y,x) for x,y in points]
    points_sort = np.array(points[np.argsort(points_convert)])
    print str(points_sort.shape[0])


    graft_points = rnd.sample(points_sort,numchains)
    #graft_points = points_sort[::int(len(points_sort)/numchains+1)]


   #graft_points[:,0] = x_val
    #graft_points[:,1] = y_val
    for i in xrange(numchainsA):
        chain = -1
        graft_index = np.intersect1d(np.where(graft_points[i][0]==points_sort[:,0])[0],np.where(graft_points[i][1]==points_sort[:,1])[0])

        while np.any(chain == -1):
            # print points_sort[graft_index]
            # print points_sort[graft_index][0][0]
            #print points_sort[graft_index][0][1]
            chain = hex_Saw(np_lattice, chainsize_A, grafted_to=(int(points_sort[graft_index%points_sort.shape[0]][0][0]), int(points_sort[graft_index%points_sort.shape[0]][0][1])))
            if np.any(chain == -1):
                print 'moving chain 2'
                graft_index += 2

        chains_A[i] = chain
        if ((i + 1) % 10 == 0):
            print "\n\n" + str(i + 1) + "/" + str(numchains) + " chains grown\n\n"
            # print "Chain "+str(i+1)+" out of "+str(numchains)+" grown"
    for i in xrange(numchainsB):
        chain = -1
        graft_index = np.intersect1d(np.where(graft_points[i][0]==points_sort[:,0])[0],np.where(graft_points[i][1]==points_sort[:,1])[0])

        while np.any(chain == -1):
            chain = hex_Saw(np_lattice, chainsize_B, grafted_to=(int(points_sort[graft_index%points_sort.shape[0]][0][0]), int(points_sort[graft_index%points_sort.shape[0]][0][1])))
            if np.any(chain == -1) :
                print 'moving chain 2'
                graft_index += 2

        addition = np.full((chainsize_A-chainsize_B,2),-1, dtype=np.int)
        chains_B[i] = np.concatenate((chain,addition),0)
        if ((i + 1) % 10 == 0):
            print "\n\n" + str(i + 1) + "/" + str(numchains) + " chains grown\n\n"
            # print "Chain "+str(i+1)+" out of "+str(numchains)+" grown"
    chains = np.concatenate((chains_A,chains_B), 0)
    #print chains[40,:,:]
#    print chains_B[10,:,:]
#    print len(chains_B[10,:,:])
    return (chains, np_lattice, points)


def cbmc(lattice, chains):  #configurational bias monte carlo
    count = 0
    chains_old = np.copy(chains)
    lattice_old = np.copy(lattice)
    energy_old = calcenergy(chains, lattice)
    chain = rnd.choice(chains)
    length = chain[chain[:,0]>=0].shape[0]
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
    if index > length:
        print 'monomer is: ' + str(monomer) + 'value at index is: ' + str(chain[index,:])
        index = length - 1
        monomer = chain[index,:]
        print 'move index back'
    numchains = length - index
    lattice[chain[index:length,0].astype(int),chain[index:length,1].astype(int)] = 0
    deleted = hex_Saw(lattice,numchains,monomer)
    if np.any(deleted == -1):
        print "no saw found"
        return (lattice_old, chains_old, energy_old,count)
    chain[index:length,:] = [0,0]
    chain[index:length,:] = deleted
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

def chains_to_xyz(chains,filename):
    #chainlength_A = chains[0].shape[0]
    #chainlength_B = int(chainlength_A*.33)+1
    chainlength_A = 12
    chainlength_B = 4
    numchains = chains.shape[0]

    xyzfile = open(filename ,'a')
    numatoms = chainlength_A*numchains/2 + chainlength_B*numchains/2;
    xyzfile.write(str(numatoms)+'\n\n')
    for chain in chains:
        chain = chain[chain[:,0] >= 0]
       # xyzfile.write('# '+str(chain))
        for monomer in chain:
            distance = math.sqrt((monomer[0]-chain[0][0])**2 + (monomer[1] - chain[0][1])**2)
            if distance > 13:
                print 'distance between tether is ' + str(distance)
            if chain.shape[0] == chainlength_A:
                xyzfile.write('C '+str(monomer[0])+' '+str(monomer[1])+' '+ '0'+'\n')
            elif chain.shape[0] == chainlength_B:
                xyzfile.write('N ' + str(monomer[0])+' '+str(monomer[1])+' '+ '0'+'\n')
            #else:
                #xyzfile.write('O '+str(monomer[0])+' '+str(monomer[1])+' '+ '0'+'\n')
    xyzfile.close()

def store_energies(energy):
    energyfile = open('EnergiesHex42b','a')
    energyfile.write(str(energy)+"\n")

def acceptance_rate(i, count):
    rate = count/float(i)
    print "\n\n The rate is: " + str(rate*100) + "%"
    return rate

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
        acceptable = chain1[chain1[:,0]>=0]
        index = acceptable.shape[0]
        numchains = length2 - index
        lattice[chain1[index-1:length2,0].astype(int),chain1[index-1:length2,1].astype(int)] = 0
        deleted = hex_Saw(lattice,numchains+1,chain1[index-1,:])
        if np.any(deleted == -1):
            print "no saw found"
            return (lattice_old, chains_old, energy_old,count)
        chain1[index-1:length2,:] = deleted
        chain2[index:length2,:] = -1

    else:
        acceptable = chain2[chain2[:,0]>=0]
        index = acceptable.shape[0]
        numchains = length1 - index
        lattice[chain2[index-1:length1,0].astype(int),chain2[index-1:length1,1].astype(int)] = 0
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

def check(lattice, chains):
    for chain in chains:
        chain = chain[chain[:,0]>=0]
    #    for monomer in chain:
    #        for j in xrange(0,chain.shape[0]-1):
    #            index = np.intersect1d(np.where((chain[:,0] == monomer[0]))[0], np.where((chain[:,1] == monomer[1]))[0])
    #    intersect = np.array([x for x in set(tuple(x) for x in chain) & set(tuple(x) for x in chain)])
        intersect = [tuple(x) for x in chain]
        dups = [item for item,count in collections.Counter(intersect).items() if count > 1]
    #            if index == j:
    #               print "\n"
        if len(dups) > 0:
            print "ERROR: Duplicate in chain"
            print dups
            index = np.intersect1d(np.where((chain[:,0] == dups[0][0])), np.where((chain[:,1] == dups[0][1])))
            print index
            #print "Monomer is: " + str(monomer) + "at index " + str(index)
    for chain1, chain2 in itertools.combinations(chains,2):
            chain1 = chain1[chain1[:,0]>=0]
            chain2 = chain2[chain2[:,0]>=0]
            array = np.concatenate((chain1, chain2), 0)
            intersect2 = [tuple(x) for x in array]
            dups2 = [item for item,count in collections.Counter(intersect2).items() if count > 1]

            if len(dups2) > 0:
                print "overlap"
                print dups2
                index2 = np.intersect1d(np.where((chain[:,0] == dups2[0][0])), np.where((chain[:,1] == dups2[0][1])))
                print index2

    return 0

def sep_analysis(chains):
    sort = np.array([0,0,0])
    zero_b = 0
    one_b = 0
    two_b = 0
    three_b = 0
    for chain in chains:
        chain = chain[chain[:,0]>=0]
        chain_len = len(chain)
        c = 0
        if chain_len  == 4:
            c = -1
        if chain_len == 12:
            c = 1
        sort = np.vstack((sort, [chain[0][0],chain[0][1], c]))
    for i in xrange(sort.shape[0]-1):
        count = 0
        if i == 0 or i == 1:
            count += 0
        else:
            type = [sort[i-1,2], sort[i,2], sort[i+1,2]]
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
    return(zero_b, one_b, two_b, three_b)


def SSR(analysis, binomial):
    ssr = 0
    for i in range(0,4):
        diff = analysis[i]-binomial[i]
        ssr += diff**2
    return ssr

if __name__ == '__main__':
    print "Starting initial configuration"
    (chains, lattice, points) = init_config((2000, 500,1), 48, 16)
    print chains
    count = 0
    chains_to_xyz(chains, 'InitHex42b')
    for i in range(0,250000):
        if rnd.uniform(0,1)<.5:
            (lattice, chains,total_energy, acc)=cbmc(lattice,chains)
            print "CBMC"
        else:
            (lattice, chains, total_energy, acc) = swap(lattice, chains)
            print "SWAP"
        store_energies(total_energy)
        chains_to_xyz(chains, 'LongHex42b')
        count += acc
        acceptance_rate(i+1,count)
<<<<<<< Updated upstream
        if i%10==0:
            np.save("Saved_ChainsHex42b", chains)
    #chains[4,4,:] = chains [20,4,:]
    #chains[4,6,:] = chains[4,10,:]
    #print chains[4,10,:]
    chains_to_xyz(chains, 'ShortHex42b')
=======
    chains_to_xyz(chains, 'Short200000c')
>>>>>>> Stashed changes
    check(lattice, chains)
    analysis = sep_analysis(chains)
    analysis = tuple(x/float(48) for x in analysis)
    print analysis
    for saw in chains:
            saw = saw[saw[:,0] >= 0]
            #saw[::2, 0] += 0.5
            plt.plot(saw[:, 0], saw[:, 1])
    plt.show()
    profile = lattice.sum(axis=0)
    plt.plot(profile[0:400])
    plt.show()
    plt.scatter(points[:,0],points[:,1])
   # plt.show()

    plt.plotfile('EnergiesHex42b')
    plt.show()

    binomial = []
    for k in xrange(0,4):
        n =3
        binomial += [choose(n,k)*(0.6**k)*(0.4**(n-k))]

    spectra = open('Saved_spectraHex40b','a')
    spectra.write(str(analysis)+"\n")

    ssr = SSR(analysis, binomial)
    SSR = open('Saved_SSRhex40b', 'a')
    SSR.write('+1\t' + str(ssr) + '\t' + '250000\n')
    SSR.close()

    S = np.loadtxt('Saved_SSRhex40b', skiprows=1)

    set = []
    for line in S:
        set += [line[1]]

    set = map(float, set)
    std = open('Standard_dev40','a')
    std.write(str(np.std(set))+'\n')

    print binomial

    t = np.array([0,1,2,3])
    plt.plot(t.T,binomial,'r')
    plt.plot(t.T,analysis,'b')
    plt.show()


