#!/bin/bash -x
./tests/spherical_head_model 0
./bin/bem_solve sphere4.fwd.prm
./bin/select_electrodes -m sphere4.surf.ucd -p sphere4_phi.dat -o sphere4_electrodes.dat 1/10 4/30
./bin/bem_inverse_solver sphere4.inv.prm | tee sphere4.inv.log
