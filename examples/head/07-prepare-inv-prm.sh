#!/bin/bash

create_inverse_params.py \
    debug=true \
    verbose=false \
    electrode_measurements=head.em1 \
    output_sources=head.source.dipoles \
    fwd_surface_mesh=head.surf.ucd \
    fwd_material_data=head.sigma \
    fwd_verbose=true \
    search_verbose=true \
    search_debug=true \
    search_initial_point='0, 0, 0' \
    search_initial_radius=10e-3 \
    search_tolerance=1e-8 \
    search_max_iterations=1000 \
    search_reflection_coefficient=1.0 \
    search_contraction_coefficient=0.5 \
    search_expansion_coefficient=2.0 \
    search_reduction_coefficient=0.5 \
    output=head.inv.prm
