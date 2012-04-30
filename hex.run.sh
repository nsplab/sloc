#!/bin/bash -x
./bin/bem_solve hex.fwd.prm
./bin/select_electrodes -m data/hex.surf.ucd -o hex.ec1 0/8
./bin/measure_electrodes -p hex_phi.dat -e hex.ec1 -o hex.em1
./bin/bem_inverse_solver hex.inv.prm | tee hex.inv.log
