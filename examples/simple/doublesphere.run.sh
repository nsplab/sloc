#!/bin/bash -x
./bin/bem_solve doublesphere.fwd.prm
./bin/select_electrodes -m data/doublesphere/surf.ucd -o doublesphere.ec1 1/20 2/10
./bin/measure_electrodes -p doublesphere_phi.dat -e doublesphere.ec1 -o doublesphere.em1
./bin/bem_inverse_solver doublesphere.inv.prm | tee doublesphere.inv.log
