#!/bin/bash -x
./bin/bem_solve doublesphere.fwd.prm
./bin/select_electrodes -m data/doublesphere.surf.ucd -p doublesphere_phi.dat -o doublesphere_electrodes.dat 1/20 2/10
./bin/bem_inverse_solver doublesphere.inv.prm | tee doublesphere.inv.log
