__author__ = 'Maggie'
import numpy as np
import matplotlib.pyplot as plt


def get_avg(file):

    S = np.loadtxt(file, skiprows=1)

    set = []
    for line in S:
        set += [line[1]]

    avg = sum(set)/len(set)
    return avg

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


def binomial_graph(file, p):
    bin = open(file, 'r')
    yvals = np.loadtxt(bin, delimiter=',')
    binomial = []
    for k in xrange(0,4):
        n =3
        binomial += [choose(n,k)*((1-p)**k)*(p**(n-k))]

    t = np.array([0,1,2,3])
    width = .35
    fig, ax = plt.subplots()
    rects1 = ax.bar(t.T,binomial,width, color ='r' )
    rects2 = ax.bar(t.T+width, yvals, width, color='b')
    #plt.plot(t.T,binomial,'r')
    print binomial
    print yvals
    #y = [p, 0, 0, (1-p)]

    #print yvals
   # for i in range(yvals.shape[0]):
   #     if i == 1:
    #        plt.plot(t.T,yvals[i],'b')
     #   elif i == 2:
      #      plt.plot(t.T,yvals[i],'g')
      #  elif i == 3:
      #      plt.plot(t.T,yvals[i],'c')
      #  else:
      #      plt.plot(t.T, yvals[i], 'm')

   # plt.plot(t.T, yvals)
    ax.set_xticks(t.T+width)
    ax.set_xticklabels( ('0', '1', '2', '3') )

    ax.legend( (rects1[0], rects2[0]), ('Expected', 'Theoretical Spectra') )

    plt.title("Binomial Distribution for a Short Chain Concentration of " + str(p))
    plt.xlabel("Number of short chains")
    plt.ylabel("Percent of total chains")
    plt.show()

if __name__ == '__main__':
    yvals = []
    Ys = []
    for i in range(10,100,10):
        if i%20 == 0:
            S = np.loadtxt("Saved_SSRhex" + str(i)+'b', skiprows=1)
            yvals += [S[1]]
            speca = np.loadtxt("Saved_spectraHex"+str(i)+'b', delimiter=',')
            ya = [float("."+str(i)), 0, 0, (1-(float("."+str(i))))]
            ysa = 0
            for i in range(0,4):
                ysa += (speca[i] - ya[i])**2
            print ysa
            Ys += [ysa]

        elif i == 50:
            SSR200 = get_avg("Saved_SSRi")
            yvals += [SSR200]
            specb = np.loadtxt("Saved_spectra250000", delimiter=',')
            yb = [float("."+str(i)), 0, 0, (1-(float("."+str(i))))]
            ysb = 0
            for i in range(0,4):
                ysb += (specb[0,i] - yb[i])**2
            print ysb
            Ys += [ysb]
        elif i%10 == 0 and (i<20 or i>80):
            T = np.loadtxt("Saved_SSRhex" + str(i) + 'b', skiprows=1)
            print T[1]
            yvals += [T[1]]
            specc = np.loadtxt("Saved_spectraHex"+str(i) + 'b', delimiter=',')
            yc = [float("."+str(i)), 0, 0, (1-(float("."+str(i))))]
            ysc = 0
            for i in range(0,4):
                ysc += (specc[i] - yc[i])**2
            print ysc
            Ys += [ysc]

    xvals = [.1042,.208,.396,.5,.604,.792,.896]

    print yvals
    print xvals
    print Ys

   # badErr = [.0164822411]*7

    plt.plot(xvals, yvals)
    #plt.errorbar(xvals, yvals, xerr=0,yerr=badErr)
    #plt.plot(xvals, Ys)
    plt.title("Change in SSR with concentration of short chains (Bad Solvent)", fontsize= 20)
    plt.xlabel("phiB")
    plt.xlim(0,1)
    plt.ylabel("SSR")
    plt.show()

    binomial_graph("Saved_spectraHex10b", .1)

    binomial_graph("Saved_spectraHex20b", .2)
    binomial_graph("Saved_spectraHex40b", .4)
    binomial_graph("Saved_spectraHex60b", .6)
    binomial_graph("Saved_spectraHex90b", .9)
    binomial_graph("Saved_spectraHex80b", .8)