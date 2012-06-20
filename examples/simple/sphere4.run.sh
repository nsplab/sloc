#!/bin/bash -x

# prepare the mesh
./tests/spherical_head_model \
    --debug=0 \
    --refinement-level=0 \
    --output-mesh=./out/sphere4.surf.ucd

# solve the forward model using a specific dipole
./bin/bem_solve sphere4.fwd.prm

# select electrodes
if [ ! -f ./sphere4.ec1 ]; then
    ./bin/select_electrodes \
        -m sphere4.surf.ucd \
        -o sphere4.ec1 \
        1/10 4/30
fi

# add noise to the electrodes
#./bin/measure_electrodes -p sphere4_phi.dat -e sphere4.ec1 -o sphere4.em1
./bin/measure_electrodes \
    -p sphere4_phi.dat \
    -e sphere4.ec1 \
    -o sphere4.em1

#./bin/bem_inverse_solver sphere4.inv.prm | tee sphere4.inv.log
