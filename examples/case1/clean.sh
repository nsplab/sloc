#!/bin/bash -x

# clean up data/ directory
rm -f data/true_*.dipole
rm -f data/sphere4_*.inv.prm
rm -f data/sphere4_*.fwd.prm

# clean up out/ directory
rm -f out/phi_*.dat
rm -f out/phi_*.vtk
rm -f out/mat_*.dat
rm -f out/electrode_indices_*.dat
rm -f out/noise_*.dat
rm -f out/electrodes_*.dat

# EOF
