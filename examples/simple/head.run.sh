#!/bin/bash -x
./bin/bem_solve head.fwd.prm
./bin/select_electrodes -m data/head.surf.ucd -o head.ec1 0/18
./bin/measure_electrodes -p head_phi.dat -e head.ec1 -o head.em1
./bin/bem_inverse_solver head.inv.prm | tee head.inv.log
