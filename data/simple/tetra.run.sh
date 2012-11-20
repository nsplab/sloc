#!/bin/bash -x

N=1
if [[ -n "$1" ]]; then
    N="$1"
fi

./make_dipoles \
    -v tetra.dipoles.vtk -o tetra.dipoles \
    0.3,0.3,0.3/0,0,1

./select_electrodes \
    -m tetra.surf.mesh -i tetra.surf.mat \
    -v tetra.electrodes.vtk -o tetra.electrodes \
    1/4

cat >tetra.fwd.prm <<DOC
# tetra.fwd.prm
set verbose = true
set debug = true

set surface_mesh = tetra.surf.mesh
set surface_mesh_materials = tetra.surf.mat
set material_data = tetra.sigma
set dipole_sources = tetra.dipoles

set output_vtk = tetra.vtk
set output_phi = tetra.phi.dat
DOC

mpiexec -n $N ./bem_forward_solver tetra.fwd.prm | tee tetra.fwd.log

./measure_electrodes -p tetra.phi.dat -e tetra.electrodes -o tetra.electrodes.dat -n 1e12

cat >tetra.inv.prm <<DOC
# tetra.inv.prm
set verbose = true
set electrodes = tetra.electrodes.dat
set output_sources = tetra.sources

subsection Forward Problem Parameters
    set verbose = false
    set material_data = tetra.sigma
    set surface_mesh = tetra.surf.mesh
    set surface_mesh_materials = tetra.surf.mat
end

subsection Simplex Search Parameters
    set verbose = true
    set debug = true
    #set initial_search_point = 1.1, 0.9, 1.2
    set initial_search_point = 0.5, 0.2, 0.3
    set initial_search_radius = 0.5
    set tolerance = 1e-16
    set max_iterations = 1000
    set reflection_coefficient = 1.0
    set contraction_coefficient = 0.5
    set expansion_coefficient = 2.0
    set reduction_coefficient = 0.5
end
DOC

#./bem_inverse_solver tetra.inv.prm | tee tetra.inv.log
#colordiff -u tetra.dipoles tetra.sources
