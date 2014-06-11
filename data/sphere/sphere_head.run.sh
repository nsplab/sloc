#!/bin/bash -x

N=1
if [[ -n "$1" ]]; then
    N="$1"
fi

# Create file "sphere.dipoles"
make_dipoles \
    -v sphere_head.dipoles.vtk -o sphere_head.dipoles \
    0.0,0.0,0.0/0,0,1

if [ ! -f sphere.electrodes ]; then
    select_electrodes \
        -m sphere.mesh -i sphere.mat \
        -v sphere_head.electrodes.vtk -o sphere_head.electrodes \
        1/4
fi

cat >sphere.fwd.prm <<DOC
# tetra.fwd.prm
set verbose = true
set debug = true

set surface_mesh = sphere_head.surf.mesh
set surface_mesh_materials = sphere_head.surf.mat
set material_data = sphere_head.surf.sigma
set dipole_sources = sphere_head.dipoles

set output_vtk = sphere_head.vtk
set output_phi = sphere_head.phi.dat
DOC

mpiexec -n 8 bem_forward_solver sphere.fwd.prm #>tetra.p$N.fwd.log
#bem_forward_solver sphere.fwd.prm #>tetra.fwd.log

#measure_electrodes -p sphere.phi.dat -e sphere.electrodes -o sphere.electrodes.dat -n 1e12

cat >sphere.inv.prm <<DOC
# tetra.inv.prm
set verbose = true
set electrodes = sphere.electrodes.dat
set output_sources = sphere.sources

subsection Forward Problem Parameters
    set verbose = false
    set material_data = sphere.sigma
    set surface_mesh = sphere.mesh
    set surface_mesh_materials = sphere.mat
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
