#!/usr/bin/env python

from math import sqrt, atan2, pi
from numpy import array, zeros

from time import sleep
FCOST = 0
evals = [0]


def simplex_search(F, P, alpha, beta, gamma, tol=1e-6, verbose=False):
    """
    Inputs:
        F       function we want to minimize
        P       list of n+1 points defining a simplex
        alpha   reflection coefficient (alpha > 0)
        beta    contraction coefficient (0 < beta < 1)
        gamma   expansion coefficient (gamma > 1)
        tol     convergence criterion

    Outputs:
        x_min   point at which F achieves a minimum
        y_min   value of F(x_min)
    """

    # space dimension
    n = len(P) - 1

    # sanity checks
    assert n == len(P[0])
    assert tol > 0
    assert n > 0

    # start by evaluating F on the simplex vertices
    Y = [F(p) for p in P]

    def argmaxmin():
        h,l = 0, 0
        yh, yl = Y[0], Y[0]
        for i,y in enumerate(Y):
            if y > yh:
                h, yh = i, y
            if y < yl:
                l, yl = i, y
        return h,l

    def centroid(k):
        """
        calculate the centroid of the simplex, excluding the k-th point
        """
        return sum(P[i] for i in xrange(n+1) if i != k) / n

    def reflection(alpha, P_bar, P_h):
        P_star = (1 + alpha) * P_bar - alpha * P_h
        y_star = F(P_star)
        if verbose:
            print "reflection:", P_star, y_star
        return (P_star, y_star)

    def expansion(gamma, P_star, P_bar):
        P_star2 = gamma * P_star + (1 - gamma) * P_bar
        y_star2 = F(P_star2)
        if verbose:
            print "expansion:", P_star2, y_star2
        return (P_star2, y_star2)

    def contraction(beta, P_h, P_bar):
        P_star2 = beta * P_h + (1 - beta) * P_bar
        y_star2 = F(P_star2)
        if verbose:
            print "contraction:", P_star2, y_star2
        return (P_star2, y_star2)

    def shrink_simplex(k):
        if verbose:
            print "shrinking towards", k, P[k]
        for i in xrange(n+1):
            if i != k:
                P[i] = (P[i] + P[k])/2
                Y[i] = F(P[i])
        return

    def y_avg():
        return sum(Y)/len(Y)

    def y_std():
        ybar = y_avg()
        yvar = sum((y-ybar)**2 for y in Y)/n
        return sqrt(yvar)

    count = [0]
    def step():
        count[0] += 1
        if verbose:
            print "starting iteration", count[0]

        h,l = argmaxmin()
        if verbose:
            print "max:", h, P[h], Y[h]
            print "min:", l, P[l], Y[l]

        P_bar = centroid(h)
        if verbose: print "centroid:", P_bar

        P_star, y_star = reflection(alpha, P_bar, P[h])

        if y_star < Y[l]:
            P_star2, y_star2 = expansion(gamma, P_star, P_bar)
            if y_star2 < Y[l]:
                if verbose: print "P_h <- P_star2 (using expansion)"
                P[h],Y[h] = P_star2, y_star2
                return
            else:
                if verbose: print "P_h <- P_star (using reflection)"
                P[h],Y[h] = P_star, y_star
                return
        else:
            # y_star >= Y[l]
            if not all(y_star > Y[i] for i in xrange(n+1) if i != h):
                if verbose: print "P_h <- P_star (using reflection)"
                P[h],Y[h] = P_star, y_star
                return
            else:
                if not (y_star > Y[h]):
                    if verbose: print "P_h <- P_star (swapping with reflection)"
                    P[h],Y[h] = P_star, y_star
                P_star2, y_star2 = contraction(beta, P[h], P_bar)
                if y_star2 > Y[h]:
                    shrink_simplex(l)
                    return
                else:
                    if verbose: print "P_h <- P_star2 (using contraction)"
                    P[h],Y[h] = P_star2, y_star2
                    return
        return

    def fitquadric():
        from numpy import dot
        from numpy.linalg import inv

        yy = zeros((n+1,n+1), dtype=float)
        for i in xrange(n+1):
            for j in xrange(n+1):
                if i == j:
                    yy[i,i] = Y[i]
                elif i < j:
                    x = (P[i] + P[j])/2
                    yy[i,j] = F((P[i] + P[j]) / 2)
                else:
                    yy[i,j] = yy[j,i]
        #print "yy =\n%s" % yy

        a = zeros((n,1), dtype=float)
        for i in xrange(n):
            a[i,0] = 2*yy[0,i+1] - (Y[i+1] + 3*Y[0])/2
        #print "a =\n%s" % a

        B = zeros((n,n), dtype=float)
        for i in xrange(n):
            for j in xrange(n):
                if i == j:
                    B[i,i] = 2*(Y[i+1] + Y[0] - 2*yy[0,i+1])
                else:
                    B[i,j] = 2*(yy[i+1,j+1] + Y[0] - yy[0,i+1] - yy[0,j+1])
        #print "B =\n%s" % B

        Q = zeros((n,n), dtype=float)
        for i in xrange(n):
            for j in xrange(n):
                Q[i,j] = P[i+1,j] - P[0,j]
        #print "Q =\n%s" % Q

        Binv = inv(B)
        #print "B^(-1) =\n%s" % Binv

        xi_min = -dot(Binv, a)
        x_min = P[0] + dot(Q, xi_min).reshape((n,))
        y_min = Y[0] + float(dot(a.T, xi_min))
        return x_min, y_min

    def show():
        for i in xrange(n+1):
            print P[i], Y[i]
        print "evals", evals[0]
        print "std", y_std()

    if verbose:
        show()
        print "-" * 30
    while y_std() > tol:
        sleep(FCOST)
        step()
        if verbose:
            show()
            print "-" * 30

    if verbose:
        x_min, y_min = fitquadric()
        print "quadric fit minimum ->", x_min, y_min

    h,l = argmaxmin()
    return P[l], Y[l]


def wedge(x, a):
    """
    Make a wedge simplex centered at point x, with size a
    """
    n = len(x)
    P = zeros((n+1,n), dtype=float)
    for i in xrange(n+1):
        for j in xrange(n):
            if i-1 != j:
                P[i,j] = x[j] - a
            else:
                P[i,j] = x[j] + a
    return P

def rosenbrock(p):
    evals[0] += 1
    x1,x2 = p
    return 100 * (x2 - x1**2)**2 + (1 - x1)**2

def powell_quartic(p):
    evals[0] += 1
    x1,x2,x3,x4 = p
    return (x1 + 10*x2)**2 + 5*(x3 - x4)**2 + (x2 - 2*x3)**4 + 10*(x1 - x4)**4

def theta(x1,x2):
    val = atan2(x2,x1)
    if x1 < 0:
        val += pi
    return val / (2 * pi)

def powell_helical_valley(p):
    evals[0] += 1
    x1,x2,x3 = p
    return 100 * (x3 - 10 * theta(x1,x2))**2 + (sqrt(x1**2 + x2**2) - 1)**2 + x3**2

def test_rosenbrock():
    # clear number of function evals
    evals[0] = 0

    # initial point
    x,y = -1.2, 1

    # make initial simplex
    P = wedge([x,y], 1)

    # what should these be?
    alpha, beta, gamma = 1, 0.5, 2

    # run the simplex search
    return simplex_search(rosenbrock, P, alpha, beta, gamma, tol=1e-12, verbose=True)

def test_powell_quartic():
    evals[0] = 0
    x,y,z,w = 3, -1, 0, 1
    P = wedge([x,y,z,w], 1.5)
    alpha, beta, gamma = 1, 0.5, 2
    return simplex_search(powell_quartic, P, alpha, beta, gamma, tol=1e-12, verbose=True)

def test_powell_helical_valley():
    evals[0] = 0
    x,y,z = -1, 0, 0
    P = wedge([x,y,z], 0.8)
    alpha, beta, gamma = 1, 0.5, 2
    return simplex_search(powell_helical_valley, P, alpha, beta, gamma, tol=1e-12, verbose=True)


if __name__ == '__main__':
    print "starting test"
    #x_min, y_min = test_rosenbrock()
    x_min, y_min = test_powell_quartic()
    #x_min, y_min = test_powell_helical_valley()
    print "used %d function evaluations" % evals[0]
    print "minimum ->", x_min, y_min


# EOF
