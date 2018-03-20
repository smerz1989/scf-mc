__author__ = 'Maggie'

import sys
import numpy
import matplotlib.pyplot as plt


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

def SSR(analysis, binomial):
    ssr = 0
    for i in range(0,4):
        diff = analysis[i]-binomial[i]
        ssr += diff**2
    return ssr

if __name__ == '__main__':

    p = sys.argv[1]
    p = float(p)
    print p
    binomial = []
    for k in xrange(0,4):
        n =3
        binomial += [choose(n,k)*((1-p)**k)*((p)**(n-k))]

    print binomial
    spectra = numpy.loadtxt(sys.argv[2], delimiter=",")
    print spectra

    ssr = SSR(spectra, binomial)
    t = numpy.array([0,1,2,3])
    plt.plot(t.T,binomial,'r')
    plt.plot(t.T,spectra,'b')
    plt.show()


    print ssr