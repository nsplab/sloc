#!/usr/bin/env python

"""
Usage: create_inverse_params.py \\
        debug=____ \\
        verbose=____ \\
        electrode_measurements=____ \\
        output_sources=____ \\
        fwd_surface_mesh=____ \\
        fwd_material_data=____ \\
        fwd_verbose=____ \\
        fwd_debug=____ \\
        search_verbose=____ \\
        search_debug=____ \\
        search_initial_point=__,__,__ \\
        search_initial_radius=____ \\
        search_tolerance=____ \\
        search_reflection_coefficient=____ \\
        search_contraction_coefficient=____ \\
        search_expansion_coefficient=____ \\
        search_reduction_coefficient=____ \\
        output=____
"""

TEMPLATE = """\
# -----------------------------------------------------------------------------
# file {output}

# inverse model parameters
set debug = {debug}
set verbose = {verbose}
set electrodes = {electrode_measurements}
set output_sources = {output_sources}

subsection Forward Problem Parameters
    set debug = {fwd_debug}
    set verbose = {fwd_verbose}
    set material_data = {fwd_material_data}
    set surface_mesh = {fwd_surface_mesh}
end

subsection Simplex Search Parameters
    set debug = {search_debug}
    set verbose = {search_verbose}
    set initial_search_point = {search_initial_point}
    set initial_search_radius = {search_initial_radius}
    set tolerance = {search_tolerance}
    set max_iterations = {search_max_iterations}
    set reflection_coefficient = {search_reflection_coefficient}
    set contraction_coefficient = {search_contraction_coefficient}
    set expansion_coefficient = {search_expansion_coefficient}
    set reduction_coefficient = {search_reduction_coefficient}
end
"""

param_defaults = dict(
    debug='false',
    verbose='false',
    fwd_debug='false',
    fwd_verbose='false',
    search_debug='false',
    search_verbose='true',
    search_initial_point='0,0,0',
    search_initial_radius='0.009',
    search_tolerance='1e-8',
    search_max_iterations='1000',
    search_reflection_coefficient='1.0',
    search_contraction_coefficient='0.5',
    search_expansion_coefficient='2.0',
    search_reduction_coefficient='0.5',
)

def main():

    import sys
    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(0)

    args = param_defaults.copy()
    args.update(dict(arg.split('=') for arg in sys.argv[1:]))

    print("Writing {0!r}".format(args['output']))
    with open(args['output'],'w') as fp:
        fp.write(TEMPLATE.format(**args))

if __name__ == '__main__':
    main()

# EOF
