#!/bin/bash -x

# refinement level 0
../../tests/spherical_head_model 0
rm -f *.stl *.sigma
mv sphere4.surf.ucd data/sphere4_lv0.surf.ucd

# refinement level 1
../../tests/spherical_head_model 1
rm -f *.stl *.sigma
mv sphere4.surf.ucd data/sphere4_lv1.surf.ucd

# EOF
