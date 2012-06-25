#!/bin/bash -x
./01-prepare-model.sh
./02-prepare-dipoles.sh
./03-prepare-fwd-prm.sh
./04-solve-fwd-problem.sh
./05-select-electrodes.sh
./06-measure-electrodes.sh
./07-prepare-inv-prm.sh
#./08-solve-inv-problem.sh

