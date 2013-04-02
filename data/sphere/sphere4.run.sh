#!/bin/bash -x

N=1
if [[ -n "$1" ]]; then
    N="$1"
fi

./make_spherical_head_model 2

# Other values to try:
#   0,0,0/1000e-9,0,0
#   0.010,-0.012,0.013/5.5e-8,-2.0e-7,1.0e-7
./make_dipoles \
    -v sphere4.dipoles.vtk -o sphere4.dipoles \
    0,0,0/1,0,0
    #0,0,0/100e-9,0,0

if [ ! -f sphere4.electrodes ]; then
    ./select_electrodes \
        -m sphere4.surf.mesh -i sphere4.surf.mat \
        -v sphere4.electrodes.vtk -o sphere4.electrodes \
        1/30 4/35 3/30
fi

cat >sphere4.fwd.prm <<DOC
# sphere4.fwd.prm
set verbose = true
set debug = true

set surface_mesh = sphere4.surf.mesh
set surface_mesh_materials = sphere4.surf.mat
set material_data = sphere4.sigma
set dipole_sources = sphere4.dipoles

set output_vtk = sphere4.phi.vtk
set output_phi = sphere4.phi.dat
DOC

mpiexec -n $N ./bem_forward_solver sphere4.fwd.prm >sphere4.p$N.fwd.log

if [ ! -f sphere4.electrodes.dat ]; then
    ./measure_electrodes \
        -p sphere4.phi.dat \
        -e sphere4.electrodes \
        -o sphere4.electrodes.dat \
        -n 1e12
fi

cat >sphere4.inv.prm <<DOC
# sphere4.inv.prm
set verbose = true
set electrodes = sphere4.electrodes.dat
set output_sources = sphere4.sources

subsection Forward Problem Parameters
    set verbose = false
    set debug = false
    set surface_mesh = sphere4.surf.mesh
    set surface_mesh_materials = sphere4.surf.mat
    set material_data = sphere4.sigma
end

subsection Simplex Search Parameters
    set verbose = true
    set debug = true
    set initial_search_point = 0.20, 0.20, 0.20
    set tolerance = 1e-8
    set max_iterations = 1000
    set reflection_coefficient = 1.0
    set contraction_coefficient = 0.5
    set expansion_coefficient = 2.0
    set reduction_coefficient = 0.5
end

DOC

#./bem_inverse_solver sphere4.inv.prm
#colordiff -u sphere4.dipoles sphere4.sources

