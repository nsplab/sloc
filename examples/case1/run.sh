#!/bin/bash -x
./case1.py prepare fwd
./case1.py run fwd > ./fwd.sh
./fwd.sh 
./case1.py prepare inv
./case1.py run inv > ./inv.sh
time ./inv.sh 2>&1 | tee inv.log
