"""
Contains all functions necessary to initialize the chain and nanoparticle positions inside the lattice.
"""

import numpy as np
import random as rnd

def hex_Saw(lattice, numsteps, moiety, grafted_to=None):
    """Finds a self avoiding walk on the given lattice with an optional specified starting point and sets the corresponding points in the lattice to the given moiety.

    Parameters
    ----------
    lattice : int array
        A numpy array containing the lattice positions of all current nanoparticles and chains
    numsteps : int
        The number of steps of the requested self-avoiding walk
    moiety : int
        The moiety of the monomers used in the self-avoiding walk
    grafted_to : int tuple
        The x,y,z coordinates of the desired starting point of the self-avoiding walk if left blank a random starting point is chosen.
    
    Returns
    -------
    int array
        Steps in the found self-avoiding walk, if no self-avoiding walk is found then -1 is returned.
    """
    (L, M, D) = lattice.shape
    nns_odd = np.array([(-1, -1, 0), (0, -1, 0), (1, -1, 0), (-1, 0, 0), (1, 0, 0), (0, 1, 0), (0,0,-1),(0,0,1)])
    nns_even = np.array([(0, -1, 0), (-1, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (-1, 1, 0),(0,0,-1),(0,0,1)])
    complete = False
    numrestarts = 0
    maxtries = 1
    tries = 0
    while (not complete) and (tries < maxtries):
        steps = np.zeros([numsteps, 3], dtype=int)
        if grafted_to is None:
            (row, col, layer) = (np.random.randint(0, L - 1), np.random.randint(0, M - 1),np.random.randint(0, D - 1))
        else:
            (row, col, layer) = grafted_to
        # Checks to see if there is already something in that lattice position
        if lattice[row,col,layer] > 0:
            return -1

        steps[0, :] = [row, col,layer]
        lattice[row, col,layer] = moiety
        for i in xrange(numsteps - 1):
            if (row % 2 == 0):
                currentNNs = np.array([[(row + nn[0]) % L, (col + nn[1]) % M, (layer + nn[2]) % D] for nn in nns_even])  # periodic bc
            else:
                currentNNs = np.array([[(row + nn[0]) % L, (col + nn[1]) % M, (layer + nn[2]) % D] for nn in nns_odd])

            possibleMoves = currentNNs[np.where((lattice[currentNNs[:, 0], currentNNs[:, 1], currentNNs[:, 2]] == 0))[0], :]

            if (possibleMoves.shape[0] > 0):
                move = rnd.choice(possibleMoves)
                complete = True
            else:
                numrestarts += 1
                lattice[steps[0:(i + 1), 0], steps[0:(i + 1), 1],steps[0:(i + 1), 2]] = 0
                complete = False
                break

            steps[i + 1, :] = move
            lattice[move[0], move[1], move[2]] = moiety
            (row, col, layer) = (move[0], move[1], move[2])
        tries += 1;

    if (not complete):
        return -1

    print str(numsteps) + " long chain grown with " + str(numrestarts) + " restarts."
    return steps

def graft_chains(lattice, graft_points, chainlength, numchains, moieties):
    """Places the given number of chains with the specified chainlength in self-avoiding walks starting at randomly chosen points chosen from the given graft points.

    Parameters
    ----------
    lattice : int array
        The system lattice specifying the current position of all chains and nanoparticles currently residing in the system.
    graft_points : int list
        A list of the lattice (x,y,z) points of all possible graft points on the nanoparticle surface.
    chainlength : int
        The desired chainlength of the grafted chains.
    numchains : int
        The number of chains to graft on the nanoparticle surface.
    moieties : int
        The integer ID number of the moiety that composes the chain.

    Returns
    -------
    chains
        A list of all the chains grafted onto the surface.
    """

    chosen_points = rnd.sample(graft_points,numchains) if numchains<=graft_points.shape[0] else -1
    chains = np.empty((numchains,chainlength,3),dtype=int)
    idx = 0

    for graft_point in (chosen_points):
        chain = -1
        point = graft_point
        while np.any(chain==-1):
            chain = hex_Saw(lattice,chainlength,moieties,point)
            if np.any(chain == -1):
                print 'Moving chain to new graft point.'
                point = rnd.choice(graft_points)
        chains[idx,:,:] = chain
        idx += 1

    return chains

# Finds the graft points on a lattice by filling in valid coordinates
#   lattice = 3D lattice created using initialize_lattice
#   radius =
def fill_cylinder(lattice, radius):
    """Fills a faceted cylindrical nanoparticle in the center of the given system lattice with the specified radius.

    Parameters
    ----------
    radius : float
        The desired radius of the nanoparticle cylinder in lattice units

    Returns
    -------
    graft_points : int list
        A list of the x,y,z coordinates of all possible grafting points on the newly created cylindrical nanoparticle.
    """
    # Finds the center of the 3D lattice by dividing each dimension by 2
    center = [lattice.shape[0]/2, lattice.shape[1]/2, lattice.shape[2]/2]

    startoffset = 0
    endoffset = 0

    # Goes through all values for a defined radius to find boundaries of the cylinder in the lattice
    for i in xrange(radius):

        # startoffset and endoffset are the left and right offsets from the center of the cylinder
        if (i + center[0]) % 2 == 0:
            startoffset += 1
        else:
            endoffset += 1

        lattice[center[0]+i,(center[1]-radius+startoffset):(center[1]+radius-endoffset),:] = 1
        lattice[center[0]-i,(center[1]-radius+startoffset):(center[1]+radius-endoffset),:] = 1
        lattice[(center[0]-i),[(center[1]-radius+startoffset-1),(center[1]+radius-endoffset)],:] = -1
        lattice[(center[0]+i),[(center[1]-radius+startoffset-1),(center[1]+radius-endoffset)],:] = -1

    lattice[center[0]+radius,(center[1]-radius+startoffset):(center[1]+radius-endoffset),:]=-1
    lattice[center[0]-radius,(center[1]-radius+startoffset):(center[1]+radius-endoffset),:]=-1
    graft_points = np.transpose(np.vstack(np.where(lattice == -1)))

    print "The graft points look like this...\n"+str(graft_points)
    lattice[graft_points[:,0],graft_points[:,1],graft_points[:,2]]=0
    return graft_points

# Preallocates the lattice and chains, and defines the graft points
#   shape = list containing 3 integers to define the dimensions of the 3D lattice
#   numchains = total number of chains of A and B
#   chainlength = total chain length of A and B
#   moieties = list of 2 integers to identify A and B monomers
#   radius =
def initialize_lattice(shape, numchains, chainlength, moieties, radius):
    """Initializes the system lattice with specified number of chains, chainlengths, and chain moieties, with a faceted nanoparticle of specified radius.
    
    Parameters
    ----------
    shape : int tuple
        A tuple of ints containing the three dimensions of the system x,y,z in lattice units
    numchains : int
        The number of chains to be placed in the system.
    chainlegnths : int tuple
        A tuple of ints containing the chain length of each of the two chain types in the system.
    moieties : int tuple
        A tuple containing the integer identifier of the moieties to be used in the simulation.
    radius : int
        An integer specifying the radius of the faceted nanoparticle in lattice units.

    Returns
    -------
    (lattice, graft_points, chains)
        Returns the given lattice filled with the chains and nanoparticles.  
        Additionally returns all graft points on the nanoparticle surface in the form of a list of points (x,y,z).  Finally returns a list of all the chains grafted onto the nanoaprticle.
    """
    # chainlengthA is 75% of chainlength, and chainlengthB is the remainder
    chainlengthA = chainlength * 3 / 4
    chainlengthB = chainlength - chainlengthA
    ## test: chainlengthB = int(chainlength/2)
    ## test: chainlengthA = chainlength - chainlengthB

    # numchainsA is 50% of numchains, and numchainsB is the remainder
    numchainsA = numchains / 2
    numchainsB = numchains - numchainsA

    # Creates a lattice with dimensions defined by shape list
    lattice = np.zeros(shape, dtype=int)

    # Finds the valid graft points on the lattice
    graft_points = fill_cylinder(lattice, radius)

    chains_A = graft_chains(lattice, graft_points, chainlengthA, numchainsA, moieties[0])
    chains_Bold = graft_chains(lattice, graft_points, chainlengthB, numchainsB, moieties[1])
    chains_B = np.empty([numchainsB, chainlengthA, 3], dtype=int)
    for i in range(numchainsB):
        addition = np.full((chainlengthA - chainlengthB, 3), -1, dtype=np.int)
        chains_B[i] = np.concatenate((chains_Bold[i], addition), 0)
        chains = np.concatenate((chains_A, chains_B), 0)
    return (lattice, graft_points, chains)
