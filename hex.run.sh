#!/bin/bash -x
./bin/bem_solve hex.fwd.prm
./bin/select_electrodes -m data/hex.surf.ucd -p hex_phi.dat -o hex_electrodes.dat 0/8
./bin/bem_inverse_solver hex.inv.prm | tee hex.inv.log
