#!/bin/bash

if [ ! -f ./head_ec1.em ]; then
    measure_electrodes \
        --phi=./out/head.phi \
        --electrodes=./head_ec1.ec \
        --snr=3.2
        --output=./head_ec1.em
fi

