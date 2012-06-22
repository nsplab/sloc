#!/bin/bash

if [ ! -f ./sphere4_dp1_ec1.em ]; then
    measure_electrodes \
        -p ./sphere4_dp1.phi \
        -e ./sphere4.ec1 \
        -n 3.2 \
        -o ./sphere4_dp1_ec1.em
fi

