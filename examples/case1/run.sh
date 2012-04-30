#!/bin/bash -x

# first, create the surface meshes we need
./surf.sh

# next, build the true solutions using the true_dpN.dipole files
./case1.py prep fwd
./case1.py run fwd > ./fwd.sh
chmod +x ./fwd.sh
time ./fwd.sh

# now try to do source localization on the noisy data
./case1.py prep inv
./case1.py run inv > ./inv.sh
chmod +x ./inv.sh
time ./inv.sh 2>&1 | tee inv.log

# EOF
