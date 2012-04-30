#!/bin/bash -x
./tests/spherical_head_model 0
./bin/bem_solve sphere4.fwd.prm
./bin/select_electrodes -m sphere4.surf.ucd -o sphere4.ec1 1/10 4/30
./bin/measure_electrodes -p sphere4_phi.dat -e sphere4.ec1 -o sphere4.em1
./bin/bem_inverse_solver sphere4.inv.prm | tee sphere4.inv.log
