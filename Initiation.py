__author__ = 'Maggie'

import numpy as np
import random as rnd

#Don't choose any negative moieties

def hex_Saw(lattice, numsteps, moiety, grafted_to=None):
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
        if lattice[row,col,layer] > 0:          #makes sure there is nothing already there
            return -1
        # print "Starting location is "+str(row)+","+str(col)
        steps[0, :] = [row, col,layer]
        lattice[row, col,layer] = moiety
        for i in xrange(numsteps - 1):
            if (row % 2 == 0):
                currentNNs = np.array([[int(row + nn[0]) % L, int(col + nn[1]) % M,int(layer + nn[2]) % D] for nn in nns_even])  # periodic bc
            else:
                currentNNs = np.array([[int(row + nn[0]) % L, int(col + nn[1]) % M,int(layer + nn[2]) % D] for nn in nns_odd])

            possibleMoves = currentNNs[np.where((lattice[currentNNs[:, 0], currentNNs[:, 1], currentNNs[:, 2]] == 0))[0], :]
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
            lattice[move[0], move[1], move[2]] = moiety
            (row, col, layer) = (move[0], move[1], move[2])
        tries += 1;
    if (not complete):
        return -1
    print str(numsteps) + " long Chain grown with " + str(numrestarts) + " restarts"
    return steps


def graft_chains(lattice,graft_points,chainlength,numchains, moiety):
    chosen_points = rnd.sample(graft_points,numchains) if numchains<=graft_points.shape[0] else -1
    chains = np.empty((numchains,chainlength,3))
    idx = 0
  #  print chosen_points[0]
    for graft_point in (chosen_points):
        chain = -1
        point = graft_point
        while np.any(chain==-1):
            chain = hex_Saw(lattice,chainlength,moiety,point)
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
    #Check functionality without
	#num_graft_points = 4*radius+2*radius
	#graft_points = np.empty((num_graft_points,3))
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


def initialize_lattice(shape,numchains,chainlength,moieties,radius):

    chainsize_A = int(chainlength*.75)
    chainsize_B = chainlength - chainsize_A
    numchainsA = int(numchains/2)
    numchainsB = numchains - numchainsA

    lattice = np.zeros(shape)
    #radius = int(0.25*lattice.shape[1])
    graft_points = fill_cylinder(lattice,radius)
    #chains = graft_chains(lattice,graft_points,chainlength,numchains)
    chains_A = graft_chains(lattice, graft_points, chainsize_A, numchainsA,moieties[0])
    chains_Bold = graft_chains(lattice, graft_points, chainsize_B, numchainsB,moieties[1])
    chains_B = np.empty([numchainsB, chainsize_A, 3])
    for i in range(numchainsB):
        addition = np.full((chainsize_A-chainsize_B,3),-1, dtype=np.int)
        chains_B[i] = np.concatenate((chains_Bold[i],addition),0)
    chains = np.concatenate((chains_A,chains_B), 0)
    return (lattice,graft_points,chains)