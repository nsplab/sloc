#!/usr/bin/env python
"""
Usage: create_forward_params.py \\
        debug=____ \\
        verbose=____ \\
        surface_mesh=____ \\
        material_data=____ \\
        dipoles=____ \\
        output_vtk=____ \\
        output_phi=____ \\
        output=____
"""

TEMPLATE = """\
# -----------------------------------------------------------------------------
# file {output}

# forward model parameters
set surface_mesh = {surface_mesh}
set material_data = {material_data}
set dipole_sources = {dipoles}
set output_vtk = {output_vtk}
set output_phi = {output_phi}
set debug = {debug}
set verbose = {verbose}
"""

def main():
    import sys

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(0)

    args = dict(verbose='false', debug='false')
    args.update(dict(arg.split('=') for arg in sys.argv[1:]))

    print("Writing {output!r}".format(**args))
    with open(args['output'],'w') as fp:
        fp.write(TEMPLATE.format(**args))

if __name__ == '__main__':
    main()

# EOF
