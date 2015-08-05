#!/usr/bin/python
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import math
import collections
import itertools
from mpl_toolkits.mplot3d import Axes3D

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
	for idx,graft_point in enumerate(chosen_points):
		chains[idx,:,:] = hex_Saw(lattice,chainlength,graft_point)
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

	
def initialize_lattice(shape,numchains,chainlength):
	lattice = np.zeros(shape)
	radius = int(0.25*lattice.shape[1])
	graft_points = fill_cylinder(lattice,radius)
	chains = graft_chains(lattice,graft_points,chainlength,numchains)
	return (lattice,graft_points,chains)



if __name__=='__main__':
	(lattice,graft_points,chains) = initialize_lattice((200,200,200),60,30)
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for chain in chains:
		ax.plot(chain[:,0],chain[:,1],chain[:,2])
	plt.show()
