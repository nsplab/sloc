#!/usr/bin/env python
"""
Usage: create_material_layers.py
        material_A=__ \\
        material_B=__ \\
        material_C=__ \\
        layer_1=__/__ \\
        layer_2=__/__ \\
        layer_3=__/__ \\
        output_file=____

"""
def main():

    import sys
    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(0)

    from decimal import Decimal
    args = dict(arg.split('=') for arg in sys.argv[1:])
    mats = dict((x.split('_')[1], Decimal(args[x])) for x in args if x.startswith('material'))
    lays = dict((int(x.split('_')[1]), args[x].split('/')) for x in args if x.startswith('layer'))

    print("Writing {output!r}".format(**args))
    with open(args['output'], 'w') as fp:
        fp.write("{num_layers}\n".format(num_layers=len(lays)))
        for l in lays:
            k = dict(layer=l, sigma_int=mats[lays[l][1]], sigma_ext=mats[lays[l][0]])
            fp.write("{layer} {sigma_int:g} {sigma_ext:g}\n".format(**k))


if __name__ == '__main__':
    main()

# EOF
