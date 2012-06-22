#!/bin/bash

create_forward_params.py \
    debug=true \
    verbose=true \
    surface_mesh=sphere4.surf.ucd \
    material_data=sphere4.sigma \
    dipoles=sphere4_dp1.dipole \
    output_vtk=sphere4_dp1.vtk \
    output_phi=sphere4_dp1.phi \
    output=sphere4_dp1.fwd.prm

