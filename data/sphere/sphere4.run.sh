#!/bin/bash -x

./make_spherical_head_model 0

#./make_dipoles sphere4.dipoles 0.010,-0.012,0.013/5.5e-8,-2.0e-7,1.0e-7

if [ ! -f sphere4.electrodes ]; then
    ./select_electrodes \
        -m sphere4.surf.mesh -i sphere4.surf.mat \
        -v sphere4.electrodes.vtk -o sphere4.electrodes \
        1/30 4/35 3/30
fi

./bem_forward_solver sphere4.fwd.prm

if [ ! -f sphere4.electrodes.dat ]; then
    ./measure_electrodes \
        -p sphere4.phi.dat \
        -e sphere4.electrodes \
        -o sphere4.electrodes.dat \
        -n 1e12
fi

./bem_inverse_solver sphere4.inv.prm

colordiff -u sphere4.dipoles sphere4.sources

