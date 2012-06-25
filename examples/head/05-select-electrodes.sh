#!/bin/bash

if [ ! -f head.ec1 ]; then
    select_electrodes \
        --input-mesh=head.surf.ucd \
        --output-vtk=head.ec1.vtk \
        --output=head.ec1 \
        1/30 4/35 3/30
fi

