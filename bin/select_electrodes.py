#!/usr/bin/env python
"""
(1) Read potentials from a file called phi.dat
(2) Read indices from a file called 'electrode_indices.dat'
(3) Write out a file called 'electrodes.txt'
"""

def read_potentials(filename):
    print("Reading {0}".format(filename))
    phi = []
    with open(filename, 'r') as fp:
        lines = fp.readlines()
        phi = [float(val) for val in lines]
    return phi

def read_indices(filename):
    print("Reading {0}".format(filename))
    indices = []
    with open(filename, 'r') as fp:
        lines = fp.readlines()
        indices = [int(val) for val in lines]
    return indices

def write_indices(filename, indices):
    print("Writing {0}".format(filename))
    with open(filename, 'w') as fp:
        M = len(indices)
        for index in indices:
            fp.write("{0}\n".format(index))
    return

def write_electrodes(filename, phi, indices):
    print("Writing {0}".format(filename))
    with open(filename, 'w') as fp:
        M = len(indices)
        fp.write("{0} 0\n".format(M))
        for index in indices:
            fp.write("{0} {1}\n".format(index, phi[index]))
    return

def main():

    import os.path
    from random import randint

    # read potential measurements from 'phi.dat'
    phi = read_potentials('phi.dat')
    N = len(phi)

    indices_filename = 'electrode_indices.dat'
    if not os.path.exists(indices_filename):
        # if we're not reading electrode indices from a file,
        # let's just choose the number of electrodes
        M = 30
        # make a list of M random integers in the range [0,N)
        # and save it to 'electrode_indices.dat'
        indices = [randint(0,N-1) for i in xrange(M)]
        write_indices(indices_filename, indices)
    else:
        # the indices are there. let's read them
        indices = read_indices(indices_filename)
        M = len(indices)

    # finally, write out the file 'electrodes.txt' that we can use
    # to run bem_inverse_solver
    write_electrodes('electrodes.txt', phi, indices)


if __name__ == '__main__':
    main()

# EOF
