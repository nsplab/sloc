#!/bin/bash -x

N=1
if [[ -n "$1" ]]; then
    N="$1"
fi

center=0,0,0

./make_dipoles \
    -v icosahedron.dipoles.vtk \
    -o icosahedron.dipoles \
    ${center}/0,0,1

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

if [ ! -f icosahedron.electrodes.dat ]; then
    ./measure_electrodes \
        -p icosahedron.p$N.phi.dat \
        -e icosahedron.electrodes \
        -o icosahedron.electrodes.dat
fi

cat >icosahedron.grid.prm <<DOC
set electrodes = icosahedron.electrodes.dat
set output_file = icosahedron.cost_at_grid_pts

subsection Forward Problem Parameters
    set surface_mesh = icosahedron.surf.mesh
    set surface_mesh_materials = icosahedron.surf.mat 
    set material_data = icosahedron.sigma
end

subsection Grid Parameters
    set grid_dims = 10, 10, 10
    set grid_lengths = 0.01, 0.01, 0.01
    set grid_center = ${center}
end
DOC

mpiexec -n $N ./bem_cost_function icosahedron.grid.prm >icosahedron.p$N.grid.log
