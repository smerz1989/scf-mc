#!/usr/bin/python
import numpy as np
import random as rnd
import matplotlib.pyplot as plt

L=1000
M=1000
hexagonal_lattice = np.empty((L,M))

def hex_Saw(lattice,numsteps,grafted_to=None):
	(L,M) = lattice.shape
	#print "Lattice at beginning is "+str(lattice)
	nns_odd = np.array([(-1,-1),(0,-1),(1,-1),(-1,0),(1,0),(0,1)])
	nns_even = np.array([(0,-1),(-1,0),(1,0),(1,1),(0,1),(-1,1)])
	complete=False
	numrestarts = 0
	while(not complete):	
		steps = np.zeros([numsteps,2],dtype=int)
		if(grafted_to==None):
			(row,col) = (np.random.randint(0,L-1),np.random.randint(0,M-1))
		else:
			(row,col) = grafted_to
		#print "Starting location is "+str(row)+","+str(col)
		steps[0,:] = [row,col]
		lattice[row,col] = 1
		for i in xrange(numsteps-1):
			if(row%2==0):
				currentNNs = np.array([[(row+nn[0])%L,(col+nn[1])%M] for nn in nns_even])
			else:
				currentNNs = np.array([[(row+nn[0])%L,(col+nn[1])%M] for nn in nns_odd])
			#for neighbor in currentNNs:
			#	print str(neighbor)+" has lattice value "+str(lattice[neighbor[0],neighbor[1]])
			possibleMoves = currentNNs[np.where((lattice[currentNNs[:,0],currentNNs[:,1]]!=1))[0],:]
			#print "Possible moves are "+str(possibleMoves)
			if(possibleMoves.shape[0]>0):
				move = rnd.choice(possibleMoves)
				complete=True
			else:
				#print "Reached dead end restarting chain"
				#print steps[0:(i+1),0]
				#print steps[0:(i+1),1]
				numrestarts+=1
				lattice[steps[0:(i+1),0],steps[0:(i+1),1]]=0
				complete=False
				break
			#print "Chosen move is "+str(move)
			steps[i+1,:] = move
			lattice[move[0],move[1]] = 1
			(row,col) = (move[0],move[1])
	print str(numsteps)+" long Chain grown with "+str(numrestarts)+" restarts"
	return steps


def fill_nanoparticle(lattice,diameter):
	print "placeholder"

def init_config(lattice_size,numchains,chainsize):
	lattice = np.empty((lattice_size[0],lattice_size[1]))
	chains = np.empty([numchains,chainsize,2])
	graft_points = np.empty((numchains,2))
	graft_step = int((int(3*lattice_size[0]/4)-int(lattice_size[0]/4))/numchains)
	lattice[int(lattice_size[0]/4):int(3*lattice_size[0]/4),int(lattice_size[1]/2)] = 1
	graft_points[:,0] = np.arange(start=int(lattice_size[0]/4),stop=int(3*lattice_size[0]/4)-5*graft_step,step=graft_step)
	graft_points[:,1] = int(lattice_size[1]/2)+1
	for i in xrange(numchains):
		chains[i] = hex_Saw(lattice,chainsize,grafted_to = (int(graft_points[i,0]),int(graft_points[i,1])))
		if((i+1)%10==0):
			print "\n\n"+str(i+1)+"/"+str(numchains)+" chains grown\n\n"
		#print "Chain "+str(i+1)+" out of "+str(numchains)+" grown"
	return (chains,lattice)

if __name__=='__main__':
	print "Starting initial configuration"
	(chains,lattice) = init_config((2000,500),120,10)
	for saw in chains:
		saw[::2,0]+=0.5
		plt.plot(saw[:,0],saw[:,1])
	plt.show()
	profile = lattice.sum(axis=0)
	plt.plot(profile[251:266])
	plt.show()
