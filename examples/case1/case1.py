#!/usr/bin/env python
"""
Generate the input files that we will be using in case1
"""

import numpy

FWD_EXE = "../../bin/bem_solve"
INV_EXE = "../../bin/bem_inverse_solver"

TEMPLATE_FWD_PRM = """\
# -----------------------------------------------------------------------------
# file sphere4_{lv}_{dp}.fwd.prm

# forward model parameters
set surface_mesh = ./data/sphere4_{lv}.surf.ucd
set material_data = ./data/sphere4.sigma
set dipole_sources = ./data/true_{dp}.dipole
set output_vtk = ./out/phi_{lv}_{dp}.vtk
set output_phi = ./out/phi_{lv}_{dp}.dat
set output_mat = ./out/mat_{lv}.dat
set debug = true
set verbose = true
"""

TEMPLATE_INV_PRM = """\
# -----------------------------------------------------------------------------
# file sphere4_{lv}_{dp}_{ec}_{snr}_{tr}.inv.prm

# inverse model parameters
set verbose = true
set electrodes = ./out/electrodes_{lv}_{dp}_{ec}_{snr}_{tr}.dat
set output_sources = ./out/source_{lv}_{dp}_{ec}_{snr}_{tr}.dipole

subsection Forward Problem Parameters
    set verbose = false
    set material_data = ./data/sphere4.sigma
    set surface_mesh = ./data/sphere4_{lv}.surf.ucd
end

subsection Simplex Search Parameters
    set verbose = true
    set debug = false
    #set initial_search_point = 0.0, 0.0, 0.0
    set initial_search_radius = 0.009
    set tolerance = 1e-8
    set max_iterations = 1000
    set reflection_coefficient = 1.0
    set contraction_coefficient = 0.5
    set expansion_coefficient = 2.0
    set reduction_coefficient = 0.5
end
"""

# -----------------------------------------------------------------------------

# make a range function we can use
fullrange = lambda n: range(1, n+1)

# number of models to run
models = ['sphere4']

# number of refinement levels for each model
levels = [("lv%d" % i) for i in [0,1]]

# number of electrode configurations
electrodes = [("ec%d" % i) for i in fullrange(2)]
total_electrodes = 20
electrodes_layout = {
    'ec1': [{'layer': 4, 'num': total_electrodes}],
    'ec2': [{'layer': 4, 'num': total_electrodes-8}, {'layer':1, 'num': 8}],
}

# number of distinct dipoles to consider
dipoles = [("dp%d" % i) for i in fullrange(3)]

# number of trials per ec/dp combinations
trials = [("tr%d" % i) for i in fullrange(3)]

# number of SNR levels (let's get 4 points per plot)
num_snr_points = 4
snr_codes = [("snr%d" % i) for i in fullrange(num_snr_points)]
snr_values = numpy.logspace(0, 5, num=num_snr_points)
snr_factor = 1 / numpy.sqrt(snr_values)

# -----------------------------------------------------------------------------

def model_fwd_prm(**kw):
    return 'data/sphere4_{lv}_{dp}.fwd.prm'.format(**kw)

def model_inv_prm(**kw):
    return 'data/sphere4_{lv}_{dp}_{ec}_{snr}_{tr}.inv.prm'.format(**kw)

def true_dipole_dat(**kw):
    return 'data/true_{dp}.dipole'.format(**kw)

def true_phi_dat(**kw):
    return 'out/phi_{lv}_{dp}.dat'.format(**kw)

def noise_dat(**kw):
    return 'out/noise_{lv}_{dp}_{ec}_{snr}_{tr}.dat'.format(**kw)

def material_dat(**kw):
    return 'out/mat_{lv}.dat'.format(**kw)

def electrode_indices_dat(**kw):
    return 'out/electrode_indices_{lv}_{ec}.dat'.format(**kw)

def electrodes_dat(**kw):
    return 'out/electrodes_{lv}_{dp}_{ec}_{snr}_{tr}.dat'.format(**kw)

def sources_dat(**kw):
    return 'out/sources_{lv}_{dp}_{ec}_{snr}_{tr}.dat'.format(**kw)

# -----------------------------------------------------------------------------

def every_combination():
    for lv in levels:
        for dp in dipoles:
            for ec in electrodes:
                for snr in snr_codes:
                    for tr in trials:
                        yield dict(lv=lv, dp=dp, ec=ec, snr=snr, tr=tr)
    return

def unpack(**kw):
    return [kw[col] for col in ('lv','dp','ec','snr','tr')]

def every_fwd_prm():
    for lv in levels:
        for dp in dipoles:
            yield model_fwd_prm(lv=lv, dp=dp)
    return

def every_inv_prm():
    for kw in every_combination():
        yield model_inv_prm(**kw)
    return

# -----------------------------------------------------------------------------

def random_dipole(radius, strength):
    """
    Generate a random dipole within a given radius and strength
    """
    from numpy.random import random

    # generate random dipole location
    # 's' is a random scale on the radial direction
    # 'u' is a randomly chosen unit vector
    s = 2 * random(3) - 1
    u = 2 * random(3) - 1
    u /= numpy.sqrt(sum(u*u))
    x,y,z = (radius * s) * u

    # generate random dipole strength
    # 's' is a random scaling factor
    # 'v' aligns the dipole components along z-axis
    s = 1 + 5 * random(1)
    v = numpy.array([0, 0, 1], dtype=float)
    px,py,pz = (strength * s) * v

    return (x,y,z,px,py,pz)

def generate_true_dipoles():
    for dp in dipoles:
        p = random_dipole(radius=70e-3, strength=1e-8)
        filename = true_dipole_dat(dp=dp)
        with open(filename, 'w') as fp:
            print('Writing ' + filename)
            fp.write('1\n')
            fp.write('{0} {1} {2} {3} {4} {5}\n'.format(*p))

def generate_fwd_params():
    for lv in levels:
        for dp in dipoles:
            kw = dict(lv=lv, dp=dp)
            filename = model_fwd_prm(**kw)
            with open(filename, 'w') as fp:
                print('Writing ' + filename)
                fp.write(TEMPLATE_FWD_PRM.format(**kw))
    pass

# -----------------------------------------------------------------------------

def generate_noise():
    """
    Pre-generate the noise we will apply to every electrode
    """
    from numpy.random import normal
    factor = dict(zip(snr_codes, snr_factor))

    # cache the potentials phi per (lv,dp) pair
    phi = dict()

    for kw in every_combination():

        lv,dp,ec,snr,tr = unpack(**kw)

        # read the true potentials
        if (lv,dp) not in phi:
            filename = true_phi_dat(lv=lv, dp=dp)
            with open(filename, 'r') as fp:
                phi[lv,dp] = numpy.array([float(x) for x in fp.readlines()])
        true_phi = phi[lv,dp]

        # generate the noise we'll be adding
        phi_rms = numpy.sqrt(sum(true_phi * true_phi) / len(true_phi))
        sigma = phi_rms * factor[snr]
        noise = normal(loc=0, scale=sigma, size=total_electrodes)

        # write it out
        filename = noise_dat(**kw)
        with open(filename, 'w') as fp:
            print('Writing ' + filename)
            for x in noise:
                fp.write('{0}\n'.format(x))
    return

def generate_electrode_indices():
    from random import sample

    for lv in levels:

        # first, read the material ids
        filename = material_dat(lv=lv)
        with open(filename, 'r') as fp:
            layer_indices = dict()
            for line in fp.readlines():
                i,m = [int(s) for s in line.split()]
                if m not in layer_indices:
                    layer_indices[m] = []
                layer_indices[m].append(i)

        for ec in electrodes:

            # obey the electrodes layout
            indices = []
            for layout in electrodes_layout[ec]:
                # select n random vertices from layer m
                m,n = layout['layer'], layout['num']
                indices += sample(layer_indices[m], n)

            # write out the indices
            filename = electrode_indices_dat(lv=lv, ec=ec)
            with open(filename, 'w') as fp:
                print('Writing ' + filename)
                for i in indices:
                    fp.write('{0}\n'.format(i))
    return

def generate_electrodes():
    """
    Generate the electrodes.dat input files needed to run the inverse solver
    """
    # cache the potentials phi per (lv,dp) pair
    phi = dict()

    # cache the indices per (lv,ec) pair
    indices = dict()

    for kw in every_combination():

        lv,dp,ec = kw['lv'], kw['dp'], kw['ec']

        # read the true potentials
        if (lv,dp) not in phi:
            filename = true_phi_dat(lv=lv, dp=dp)
            with open(filename, 'r') as fp:
                phi[lv,dp] = numpy.array([float(x) for x in fp.readlines()])
        true_phi = phi[lv,dp]

        # read the electrode indices
        if (lv,ec) not in indices:
            filename = electrode_indices_dat(lv=lv, ec=ec)
            with open(filename, 'r') as fp:
                indices[lv,ec] = [int(x) for x in fp.readlines()]
        idx = indices[lv,ec]

        # read the pre-generated noise
        filename = noise_dat(**kw)
        with open(filename, 'r') as fp:
            noise_values = numpy.array([float(x) for x in fp.readlines()])
        assert len(noise_values) == len(idx)

        # apply the noise to phi
        filename = electrodes_dat(**kw)
        with open(filename, 'w') as fp:
            print('Writing ' + filename)
            fp.write('{0} 0\n'.format(len(idx)))
            for i,noise in zip(idx, noise_values):
                fp.write('{0} {1}\n'.format(i, true_phi[i] + noise))

    return

def generate_inv_params():
    for kw in every_combination():
        filename = model_inv_prm(**kw)
        with open(filename, 'w') as fp:
            print('Writing ' + filename)
            fp.write(TEMPLATE_INV_PRM.format(**kw))
    return

# -----------------------------------------------------------------------------

def prepare_fwd_problem():
    generate_true_dipoles()
    generate_fwd_params()

def run_fwd_problem(dry=False):
    import os
    if dry:
        print("#!/bin/bash -x")
    for filename in every_fwd_prm():
        cmd = "{exe} {prm}".format(exe=FWD_EXE, prm=filename)
        print(cmd)
        if not dry:
            os.system(cmd)
    return

def prepare_inv_problem():
    generate_noise()
    generate_electrode_indices()
    generate_electrodes()
    generate_inv_params()

def run_inv_problem(dry=True):
    import os
    if dry:
        print("#!/bin/bash -x")
    for filename in every_inv_prm():
        cmd = "{exe} {prm}".format(exe=INV_EXE, prm=filename)
        print(cmd)
        if not dry:
            os.system(cmd)
    return

# -----------------------------------------------------------------------------

def main():
    import sys

    if len(sys.argv) <= 2:
        print("Usage: {0} {{prepare,run}} {{fwd,inv}}".format(sys.argv[0]))
        sys.exit(1)

    action, prob = sys.argv[1:3]
    dry = True
    if action == 'prepare':
        if prob == 'fwd':
            prepare_fwd_problem()
        elif prob == 'inv':
            prepare_inv_problem()
        else:
            print("Uknown problem {0!r}".format(prob))
            sys.exit(2)
    elif action == 'run':
        if prob == 'fwd':
            run_fwd_problem(dry=dry)
        elif prob == 'inv':
            run_inv_problem(dry=dry)
        else:
            print("Uknown problem {0!r}".format(prob))
            sys.exit(2)
    else:
        print("Unkown action {0!r}".format(action))
        sys.exit(2)

    return

if __name__ == '__main__':
    main()

# EOF
