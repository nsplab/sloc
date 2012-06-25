#!/bin/bash

create_inverse_params.py \
    debug=true \
    verbose=true \
    electrode_measurements=sphere4_dp1_ec1.em \
    output_sources=sphere4_dp1_ec1.sources \
    fwd_surface_mesh=sphere4.surf.ucd \
    fwd_material_data=sphere4.sigma \
    fwd_verbose=false \
    fwd_debug=false \
    search_verbose=true \
    search_debug=false \
    search_initial_point='0, 0, 0' \
    search_initial_radius=10e-3 \
    search_tolerance=1e-8 \
    search_max_iterations=1000 \
    search_reflection_coefficient=1.0 \
    search_contraction_coefficient=0.5 \
    search_expansion_coefficient=2.0 \
    search_reduction_coefficient=0.5 \
    output=sphere4_dp1_ec1.inv.prm
