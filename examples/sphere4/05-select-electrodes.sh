#!/bin/bash

if [ ! -f sphere4.ec1 ]; then
    select_electrodes \
        -m sphere4.surf.ucd \
        -v sphere4.ec1.vtk \
        -o sphere4.ec1 \
        1/30 4/35 3/30
fi

