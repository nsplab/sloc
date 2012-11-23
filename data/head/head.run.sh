#!/bin/bash -x

N=1
if [[ -n "$1" ]]; then
    N="$1"
fi

./make_head_model

./make_dipoles \
    -v head.dipoles.vtk -o head.dipoles \
    4.81851,-129.829,-117.922/0,0,1000e-9

if [ ! -f head.electrodes ]; then
    ./select_electrodes \
        -m head.surf.mesh -i head.surf.mat \
        -v head.electrodes.vtk -o head.electrodes \
        0/30 3/5 4/0
        #0/30 1/35 3/30
        #1/30 4/35 3/30
fi

cat >head.fwd.prm <<DOC
# head.fwd.prm
set verbose = true
set debug = true

set surface_mesh = head.surf.mesh
set surface_mesh_materials = head.surf.mat
set material_data = head.surf.sigma
set dipole_sources = head.dipoles

set output_vtk = head.phi.vtk
set output_phi = head.phi.dat
DOC

mpiexec -n $N ./bem_forward_solver head.fwd.prm >head.p$N.fwd.log

if [ ! -f head.electrodes.dat ]; then
    ./measure_electrodes \
        -p head.phi.dat \
        -e head.electrodes \
        -o head.electrodes.dat \
        -n 1e12
fi

cat >head.inv.prm <<DOC
# head.inv.prm
set verbose = true
set electrodes = head.electrodes.dat
set output_sources = head.sources

subsection Forward Problem Parameters
    set verbose = false
    set debug = false
    set surface_mesh = head.surf.mesh
    set surface_mesh_materials = head.surf.mat
    set material_data = head.surf.sigma
end

subsection Simplex Search Parameters
    set verbose = true
    set debug = true
    set initial_search_point = 3.5, -127, -119
    set tolerance = 1e-8
    set max_iterations = 1000
    set reflection_coefficient = 1.0
    set contraction_coefficient = 0.5
    set expansion_coefficient = 2.0
    set reduction_coefficient = 0.5
end
DOC

#./bem_inverse_solver head.inv.prm
#colordiff -u head.dipoles head.sources

