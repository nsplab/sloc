#!/bin/bash -x

sep
time stl_to_ucd brain48.stl brain48.tri.ucd
time tri_to_quad brain48.tri.ucd brain16.ucd
ln -sf brain16.ucd brain.ucd

sep
time stl_to_ucd skull48.stl skull48.tri.ucd
time tri_to_quad skull48.tri.ucd skull16.ucd
ln -sf skull16.ucd skull.ucd

sep
time stl_to_ucd skin96.stl skin96.tri.ucd
time tri_to_quad skin96.tri.ucd skin32.ucd
ln -sf skin32.ucd scalp.ucd

sep
rm -f *.tri.ucd
