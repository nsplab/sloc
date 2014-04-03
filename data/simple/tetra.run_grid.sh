#!/bin/bash -x

N=2
if [[ -n "$1" ]]; then
    N="$1"
fi

# Create file "tetra.dipoles"
make_dipoles \
    -v tetra.dipoles.vtk -o tetra.dipoles \
    0.3,0.3,0.3/0,0,1.0

if [ ! -f tetra.electrodes ]; then
    select_electrodes \
        -m tetra.surf.mesh -i tetra.surf.mat \
        -v tetra.electrodes.vtk -o tetra.electrodes \
        1/4
fi

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

mpiexec -n $N bem_forward_solver tetra.fwd.prm >tetra.p$N.fwd.log
#bem_forward_solver tetra.fwd.prm >tetra.fwd.log

#measure_electrodes -p tetra.phi.dat -e tetra.electrodes -o tetra.electrodes.dat -n 1e12

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

#bem_inverse_solver tetra.inv.prm | tee tetra.inv.log
#colordiff -u tetra.dipoles tetra.sources

cat >tetra.grid.prm <<DOC
set electrodes = tetra.electrodes.dat
set output_file = tetra.cost_at_grid_pts

subsection Forward Problem Parameters
    set verbose = true
    set debug = true
    set surface_mesh = tetra.surf.mesh
    set surface_mesh_materials = tetra.surf.mat
    set material_data = tetra.sigma
    set logfile = tetra.cost.log
end

subsection Grid Parameters
    set grid_dims = 1, 1, 1
    set grid_lengths = 0.01, 0.01, 0.01
    set grid_center = 0.3, 0.3, 0.3
end
DOC

#mpiexec -n $N bem_cost_function head.grid.prm #>head.p$N.grid.log
#mpiexec -n 8 bem_cost_function head.grid.prm #>head.p$N.grid.log
#bem_cost_function tetra.grid.prm #>tetra.grid.log
