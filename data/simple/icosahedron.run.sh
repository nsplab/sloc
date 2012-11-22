#!/bin/bash -x

N=1
if [[ -n "$1" ]]; then
    N="$1"
fi

./make_dipoles \
    -v icosahedron.dipoles.vtk -o icosahedron.dipoles \
    0,0,0/0,0,1

if [ ! -f icosahedron.electrodes ]; then
    ./select_electrodes \
        -m icosahedron.surf.mesh -i icosahedron.surf.mat \
        -v icosahedron.electrodes.vtk -o icosahedron.electrodes \
        4/12
fi

cat >icosahedron.fwd.prm <<DOC
# icosahedron.fwd.prm
set verbose = true
set debug = true

set surface_mesh = icosahedron.surf.mesh
set surface_mesh_materials = icosahedron.surf.mat 
set material_data = icosahedron.sigma
set dipole_sources = icosahedron.dipoles

set output_vtk = icosahedron.p$N.phi.vtk
set output_phi = icosahedron.p$N.phi.dat
DOC

mpiexec -n $N ./bem_forward_solver icosahedron.fwd.prm >icosahedron.p$N.fwd.log

