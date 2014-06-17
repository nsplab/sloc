#!/bin/bash -x

N=8

if [[ -n "$1" ]]; then
    N="$1"
fi

#center=4.81851,-129.829,-117.922
#center=28.82274,-136.70355,-159.76381
center=41.89535,-129.70389,-148.20721

make_head_model

#loop over grid points in temporal lobe
#for X in 0 
#do
#for Y in 0
#do
#for Z in 0 
#do

#for IJK in 0 
#do

#direct="0,0,0"
#if [ $IJK -eq 0 ]
#then
direct="0,0,1"
#fi

#if [ $IJK -eq 1 ]
#then
#     direct="0,1,0"
#fi

#if [ $IJK -eq 2 ]
#then
#     direct="0,0,1"
#fi

# X: 
# Y:
# Z: 

X=5
Y=5
Z=5

nz=$(echo "-148.20721 + $Z * -2.94812" | bc)
ny=$(echo "-129.70389 + $Y * -2.94812" | bc)
nx=$(echo "41.89535 + $X * -2.94812" | bc)

make_dipoles \
    -v head.dipoles.vtk -o head.dipoles \
     $nx,$ny,$nz/$direct
#     28.82274,-136.70355,-159.76381/1,0,0
#    4.81851,-129.829,-117.922/1,0,0
#    4.81851,-129.829,-117.922/0,0,1000e-9

if [ ! -f head.electrodes ]; then
    # select 10-20 electrodes
    select_electrodes_given_3d_pos \
        -m head.surf.mesh -i head.surf.mat \
        -v head.electrodes.vtk -o head.electrodes \
	22.66097/-215.17047/-106.02161 \
	-14.51707/-215.99734/-107.96710 \
	64.50138/-180.06104/-115.39294 \
	44.53565/-183.85083/-77.47854 \
	1.87596/-196.02194/-68.69604 \
	-37.78682/-183.39616/-77.88834 \
	-60.20803/-176.53879/-119.37214 \
	77.96370/-132.19577/-123.96312 \
	48.23156/-141.30139/-53.65706 \
	-1.38297/-144.61928/-40.58036 \
	-51.34485/-132.28114/-63.19135 \
	-73.09639/-129.91707/-126.25740 \
	73.95229/-94.09358/-108.64140 \
	44.53437/-85.88811/-52.85554 \
	4.47472/-87.27579/-40.55307 \
	-36.97533/-86.45993/-55.13177 \
	-60.07651/-86.27083/-120.12472 \
	42.53815/-49.62876/-98.95287 \
	-11.83266/-44.56505/-99.30063 
fi

echo $X.$Y.$Z.$direct
echo $nx.$ny.$nz


cat >head.fwd.prm <<DOC
# head.fwd.prm
set verbose = true
set debug = false

set surface_mesh = head.surf.mesh
set surface_mesh_materials = head.surf.mat
set material_data = head.surf.sigma
set dipole_sources = head.dipoles

set output_vtk = head.phi.$X.$Y.$Z.$direct.vtk
set output_phi = head.phi.$X.$Y.$Z.$direct.dat
DOC

time mpiexec -n $N bem_forward_solver head.fwd.prm #>head.p$N.fwd.log

#done
#done
#done
#done

if [ ! -f head.electrodes.dat ]; then
    measure_electrodes \
        -p head.phi.dat \
        -e head.electrodes \
        -o head.electrodes.dat \
        -n 1e12
fi

<<DOC
cat >head.inv.prm
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

#bem_inverse_solver head.inv.prm
#colordiff -u head.dipoles head.sources

cat >head.grid.prm <<DOC
set electrodes = head.electrodes.dat
set output_file = head.cost_at_grid_pts

subsection Forward Problem Parameters
    set verbose = true
    set debug = false
    set surface_mesh = head.surf.mesh
    set surface_mesh_materials = head.surf.mat
    set material_data = head.surf.sigma
end

subsection Grid Parameters
    set grid_dims = 1, 1, 1
    set grid_lengths = 0.01, 0.01, 0.01
    set grid_center = ${center}
end
DOC

#mpiexec -n $N bem_cost_function head.grid.prm #>head.p$N.grid.log
#mpiexec -n 8 bem_cost_function head.grid.prm #>head.p$N.grid.log
#bem_cost_function_omp head.grid.prm #>head.p$N.grid.log

