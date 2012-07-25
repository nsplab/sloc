#!/bin/bash

create_forward_params.py \
    debug=true \
    verbose=true \
    surface_mesh=head.surf.ucd \
    material_data=head.sigma \
    dipoles=head.dipoles \
    output_vtk=head.vtk \
    output_phi=head_phi.dat \
    output=head_dp1.fwd.prm

