#!/bin/bash -x

#./make_dipoles tetra.dipoles 1,1,1/1,1,1
./select_electrodes -m tetra.surf.mesh -i tetra.surf.mat -v tetra.electrodes.vtk -o tetra.electrodes 1/4
./bem_forward_solver tetra.fwd.prm
./measure_electrodes -p tetra.phi.dat -e tetra.electrodes -o tetra.electrodes.dat -n 1e12
./bem_inverse_solver tetra.inv.prm
colordiff -u tetra.dipoles tetra.sources
