===============
README for sloc
===============

For installation instructions refer to ``doc/hoffman2.rst``.

===============
README for sloc
===============

5/5/13 (Ram)

Documentation is scanty.

An example overview of the analysis pipeline using an 

icosahedron mesh is provided in data/simple/icosahedron.run.sh

prepare the dipoles
select electrodes into a file
run the forward forward solution using bem_forward_solver. 

this produces output_vtk (for visualization) and output_phi 

(the raw electrodepotentials)
run measure_electrodes to add noise to the simulated 

measurement - this takes an argument that specifies SNR
  icosahedron.electrodes (the electrode locations)
  icosahedron.electrodes.dat (potential measurements at those 

electrode locations

run the bem_cost_function using as input the following files:
  icosahedron.electrodes.dat
  icosahedron.surf.mesh (surface mesh specification)
  icosahedron.sigma (conductivity values)


Details on the file 'bin/bem_cost_function.cc'.  This file 

iterates through points in the simulated brain to determine 

the cost of asserting that those points are the seizure 

location.

Other parts of this project include 
(a) the meshes and file formats that determine the various 

surfaces (scalp, skull-outer, brain-outer, ventricles, 

vessels).  meshes are visualized using meshlab.  the e-field 

projected onto the mesh is visualized using paraview, which 

reads the *.vtk file produced by 'bin/bem_cost_function.cc'.

---

===============
Units
===============

===============
Multiplying dipole magnitude by a constant 
===============
Let |phi| be the solution of the forward problem with dipole p at location r. 
|rtrue| and |ptrue| are the location and the magnitude of the dipole use in the
 forward problem to simulate the electrode measurments.

.. |phi| image:: https://github.com/nsplab/sloc/blob/master/doc/phi.png?raw=true 
.. |rtrue| image:: https://github.com/nsplab/sloc/blob/master/doc/rtrue.png?raw=true 
.. |ptrue| image:: https://github.com/nsplab/sloc/blob/master/doc/ptrue.png?raw=true 

