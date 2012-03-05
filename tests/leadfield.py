import numpy
import random
import os

# First calculate the potential for a dipole source at the origin with some
# dipole orientation.
object = 'spheres'
dipoledata = open('data/{0}.dipole'.format(object), 'w')
dipole_orientation = '5 7 3.5'

# Calculate the field for the given orientation.
dipoledata.write('1\n')
dipoledata.write('0 0 0 {0}'.format(dipole_orientation))
dipoledata.close()
os.system('bin/bem_solve {0}.prm'.format(object))
os.system('mv phi.dat real_phi.dat')

# List of positions and orientations used for building L.
dipole_positions = ['0 0 0']
dipole_orientations = [' 1 0 0', ' 0 1 0', ' 0 0 1']

# Build list of dipole positions: a sphere of radius 6, spaced on a cube.
# Use a volume mesh to build it?
#N = 4
#dipole_positions = []
#for x in numpy.linspace(-3,3,N):
#    for y in numpy.linspace(-3,3,N):
#        for z in numpy.linspace(-3,3,N):
#            if (x**2 + y**2 + z**2) < 9:
#                dipole_positions.append(str(x) + ' ' + str(y) + ' ' + str(z))

# Initialize the lead field matrix.
m = 3
n = len(dipole_positions)
lead_field_matrix = numpy.zeros((m, 3*n))

#m = len(open('real_phi.dat').readlines())
# Pick m electrode points at random.
electrodes = random.sample(range(len(open('real_phi.dat').readlines())), m)

# Build the matrix.
for i, position in enumerate(dipole_positions):
    for j, orientation in enumerate(dipole_orientations):
        file = open('data/{0}.dipole'.format(object), 'w')
        file.write('1\n')
        file.write(position)
        file.write(orientation)
        file.close()

        os.system('bin/bem_solve {0}.prm'.format(object))

        # Potentials are a list in a file. Fill in columns.
        lead_field_matrix[:,3*i+j] = [open('phi.dat').readlines()[x] for x in electrodes]

# Solve for p using linear least squares. L*p = phi -> p = pinv(L)*phi
#phi = [float(i) for i in open('real_phi.dat')]
phi = numpy.array([float(open('real_phi.dat').readlines()[x]) for x in electrodes])
p_dipoles = numpy.dot(numpy.linalg.pinv(lead_field_matrix), phi)

print 'Electrode numbers', electrodes
print 'Approximate dipoles:', p_dipoles
print 'Difference:', [float(comp) for comp in dipole_orientation.split()] - p_dipoles

output = open('lead_field.dat', 'w')
for row in lead_field_matrix:
    for element in row:
        output.write(str(element) + ' ')
    output.write('\n')
