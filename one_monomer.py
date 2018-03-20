import numpy as np
import json



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

    return (lattice,graft_points)



if __name__=='__main__':
    (lattice,graft_points) = initialize_lattice((200,200,20),779,16,20)

    next = np.loadtxt('Short3D_90000c', dtype=str, delimiter=" ", skiprows = 1)
    print next
    shell = []

    for monomer in next:
        #print int(float(monomer[1]))
        for point in graft_points:
            #print point[0]
            if point[0] == int(float(monomer[1])) and point[1] == int(float(monomer[2])) and point[2] == int(float(monomer[3])):
                shell += [monomer]

    print shell

    xyzfile = open('corrected_S90000c','a')
    xyzfile.write(str(1295)+'\n\n')
    for item in shell:
        print item
        xyzfile.write(str(item[0]) + " " + str(item[1])+ ' ' + str(item[2]) + ' ' + str(item[3]) + '\n')
    xyzfile.close()

    print xyzfile