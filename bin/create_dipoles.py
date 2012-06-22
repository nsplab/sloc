#!/usr/bin/env python

"""
Usage: create_dipoles.py \\
        <dipole-output-file> \\
        <x1,y1,z1/px1,py1,pz1> \\
        <x2,y2,z2/px2,py2,pz2> \\
        ...
"""

def main():

    import sys
    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(0)

    from decimal import Decimal
    output = sys.argv[1]
    dipole_args = sys.argv[2:]

    dipoles = []
    for dp in dipole_args:
        xs,ps = dp.split('/')
        x = [Decimal(xcoord) for xcoord in xs.split(',')]
        p = [Decimal(pcoord) for pcoord in ps.split(',')]
        assert len(x) == 3
        assert len(p) == 3
        dipoles.append((x,p))

    print("Writing {0!r}".format(output))
    with open(output,"w") as fp:
        fp.write("{num_dipoles}\n".format(num_dipoles=len(dipoles)))
        for x,p in dipoles:
            fp.write("{x[0]:g} {x[1]:g} {x[2]:g} {p[0]:g} {p[1]:g} {p[2]:g}\n".format(x=x,p=p))

if __name__ == '__main__':
    main()

# EOF
