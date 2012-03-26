#!/bin/bash -x
./case1.py prep fwd
./case1.py run fwd > ./fwd.sh
./fwd.sh 
./case1.py prep inv
./case1.py run inv > ./inv.sh
time ./inv.sh 2>&1 | tee inv.log
