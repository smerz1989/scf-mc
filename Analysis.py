#!/usr/bin/python
import numpy as np
import math

# Creates coefficients for the binomial distribution (Andrew Dalke)
#   n = number of trials
#   k = number of successes
def choose(n, k):
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

# Calculates the free energy of a given configuration
#   chains = array of chain coordinates for the lattice
#   lattice = matrix of available lattice sites
def calcenergy(chains, lattice):
    (L, M, P) = lattice.shape
    nns_odd = np.array([(-1, -1, 0), (0, -1, 0), (1, -1, 0), (-1, 0, 0), (1, 0, 0), (0, 1, 0), (0,0,-1),(0,0,1)])
    nns_even = np.array([(0, -1, 0), (-1, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (-1, 1, 0),(0,0,-1),(0,0,1)])
    total_energy = 0
    chiAS = 0
    chiAB = 0.25
    chiAA = 0
    chiBB = 0
    chi = np.array([[0,0,chiAS,chiAS],[0,0,0,0],[chiAS,0,chiAA,(chiAB/2.0)],[chiAS,0,(chiAB/2.0),chiBB]])
    for chain in chains:
        #print chain
        length = chain.shape[0]
        chainlength = chain[chain[:,0]>=0].shape[0]
        for monomer in chain[0:chainlength,:]:
            if monomer[0]%2 == 0:
                    nns = nns_even
            else:
                    nns = nns_odd
            neighbors = [[int((neighbor[0]+monomer[0])%(L-1)) ,int((neighbor[1]+monomer[1])%(M-1)) , int((neighbor[2]+monomer[2])%(P-1)) ] for neighbor in nns]
            for neighbor in neighbors:
                total_energy += chi[int(lattice[int(monomer[0]),int(monomer[1]),int(monomer[2])]) ,int(lattice[int(neighbor[0]),int(neighbor[1]),int(neighbor[2])])]
                #total_energy += chiAS/2.0

    """ if chainlength < length:
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
                        total_energy += chiAS/2.0"""
    return total_energy

# Converts the chains array into an xyz file with carbon labeled as the long chains and nitrogen as the short chains
#   chains = array of coordiantes of each monomer
#   filename = name of file to be created
def chains_to_xyz(chains,filename,lattice):
    chainlength_A = chains[0].shape[0]
    chainlength_B = int(chainlength_A*.33)+1
    numchains = chains.shape[0]

    xyzfile = open(filename ,'a')
    numatoms = chainlength_A*int(numchains/2) + chainlength_B*(numchains-int(numchains/2));
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
                moietyA = lattice[int(monomer[0]), int(monomer[1]), int(monomer[2])]
                xyzfile.write(str(moietyA)+ " " +str(monomer[0])+' '+str(monomer[1])+' '+ str(monomer[2])+'\n')
            elif chain.shape[0] == chainlength_B:
                moietyB = lattice[int(monomer[0]), int(monomer[1]), int(monomer[2])]
                xyzfile.write(str(moietyB)+" " + str(monomer[0])+' '+str(monomer[1])+' '+ str(monomer[2])+'\n')
            else:
                xyzfile.write('O '+str(monomer[0])+' '+str(monomer[1])+' '+ '0'+'\n')
    xyzfile.close()

# Writes the energies to a file
#   energy = the value calculated using calcenergy
#   n = number of trial moves
def store_energies(energy, n, alph):
    energyfile = open('EnergiesDual_'+str(n)+alph,'a')
    energyfile.write(str(energy)+"\n")

# Number of moves that were accepted
#   i = number of trial moves taken
#   count = number of trial moves accepted
def acceptance_rate(i, count):
    rate = count/float(i)
    print "\n\n The rate is: " + str(rate*100) + "%"
    return rate

# Sorts segments into bins based on nearest neighbors (simulated MAULDI, binomial)
#   chains = array of chain coordinates
def sep_analysis(chains):
    count = 0
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

# Takes the sum of square residual between the expected binomial and analysis from sep_analysis
#   analysis = the distribution using bins
#   binomial = binomial distribution
def SSR(analysis, binomial):
    ssr = 0
    for i in range(0,6):
        diff = analysis[i]-binomial[i]
        ssr += diff**2
    return ssr