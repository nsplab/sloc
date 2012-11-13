#!/bin/bash -x

if [ ! -f head.electrodes ]; then
    ./select_electrodes \
        -m head.surf.mesh -i head.surf.mat \
        -v head.electrodes.vtk -o head.electrodes \
        0/30 3/5 4/0
        #0/30 1/35 3/30
        #1/30 4/35 3/30
fi

./bem_forward_solver head.fwd.prm

if [ ! -f head.electrodes.dat ]; then
    ./measure_electrodes \
        -p head.phi.dat \
        -e head.electrodes \
        -o head.electrodes.dat \
        -n 1e12
fi

#./bem_inverse_solver head.inv.prm
#colordiff -u head.dipoles head.sources

