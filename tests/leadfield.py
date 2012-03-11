import numpy
import random
import os
import sys
from pylab import errorbar, xlabel, ylabel, show

sys.path.append(u'/home/louis/dev/sloc/tests')
import nelder_mead

# Calculate the actual potential using the bem_solve function.
def real_potential(object, position, orientation):
    p_string = str(position[0]) + ' ' + str(position[1]) + ' ' + str(position[2])
    o_string = str(orientation[0]) + ' ' + str(orientation[1]) + ' ' + str(orientation[2])
    dipoledata = open('data/{0}.dipole'.format(object), 'w')
    dipoledata.write('1\n')
    dipoledata.write('{0} {1}'.format(p_string, o_string))
    dipoledata.close()
    os.system('bin/bem_solve {0}.prm'.format(object))
    os.system('mv phi.dat real_phi.dat')

# Calculates the m x 3 lead field matrix for a given position and list of m electrode points.
def lead_field_matrix(object, position, electrodes):
    m = len(electrodes)
    L = numpy.zeros((m, 3))
    unit_vectors = [' 1 0 0', ' 0 1 0', ' 0 0 1']
    p_string = str(position[0]) + ' ' + str(position[1]) + ' ' + str(position[2])
    for i, o_string in enumerate(unit_vectors):
        file = open('data/{0}.dipole'.format(object), 'w')
        file.write('1\n')
        file.write(p_string)
        file.write(o_string)
        file.close()
        
        os.system('bin/bem_solve {0}.prm'.format(object))

        # Potentials are a list in a file. Fill in columns.
        L[:,i] = [open('phi.dat').readlines()[x] for x in electrodes]
    return L

# Find the cost function given a point x.
def cost_function(x):
    evals[0] += 1

    L = lead_field_matrix(mesh, x, electrodes)

    p = numpy.dot(numpy.linalg.pinv(L), phi+error)
    phi_hat = numpy.dot(L, p)
    cost = numpy.dot(phi-phi_hat, phi-phi_hat)**0.5
    cost_file.write('{0} {1} {2} {3}\n'.format(x[0], x[1], x[2], cost))
    return cost

def dipole_fit(x):
    L = lead_field_matrix(mesh, x, electrodes)
    p = numpy.dot(numpy.linalg.pinv(L), phi+error)
    return p

def minimize_position():
    evals[0] = 0
    x, y, z = 0.5, 0.5, 0.5
    P = nelder_mead.wedge([x, y, z], 0.01)
    alpha, beta, gamma = 1, 0.5, 1.5
    x_min, y_min = nelder_mead.simplex_search(cost_function, P, alpha, beta, gamma, tol=1e-12, verbose=True)
    return x_min, y_min


mesh = 'hex'
evals = [0]

# Set up actual result.
real_position = [.5, .5, .5]
dipole_orientation = [5, 7, 10]
real_potential(mesh, real_position, dipole_orientation)

# Pick m electrode points.
m = 4
electrodes = [int(i) for i in numpy.linspace(0, len(open('real_phi.dat').readlines())-1, m)]

phi = numpy.array([float(open('real_phi.dat').readlines()[x]) for x in electrodes])

cost_file = open('cost.dat', 'w')

# Introduce error.
SNR = numpy.linspace(50,100,1)
sigma = (numpy.dot(phi,phi) / SNR)**0.5
n = 1
position_error = numpy.zeros((n,len(sigma)))
for i in range(n):
    for j, s in enumerate(sigma):
        error = numpy.random.normal(0, s, len(phi))
        x_min, y_min = minimize_position()
        position_error[i,j] = (numpy.dot(x_min - real_position, 
                                        x_min - real_position)/len(phi))**0.5

cost_file.close()
"""
# n measurements
n = 10
dipole_error = numpy.zeros((n, len(sigma)))
for i in range(n):
    for j, s in enumerate(sigma):
        error = numpy.random.normal(0, s, len(phi))
        diff = dipole_fit(real_position) - dipole_orientation
        # RMS error is this, I think.
        dipole_error[i, j] = (numpy.dot(diff, diff)/len(phi))**0.5
"""


#print 'Starting simplex search:'
#x_min, y_min = minimize_position()
print 'used {0} function evaluations'.format(evals[0])
print 'minimum ->', x_min, y_min

errorbar(SNR, numpy.mean(position_error, axis=0),
         yerr=numpy.std(position_error, axis=0))

xlabel('SNR')
ylabel('RMS position error')
show()    
