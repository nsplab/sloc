README for sloc
===============

For installation instructions refer to ``doc/hoffman2.rst``.

** Manuscript Drafts **

### Endovascular Source Localization (simulation study) #

https://www.writelatex.com/784824tywgtg#/1721178/

### Derivation of boundary element method (BEM) equation used in forward model, based on Luis' notes #

https://www.writelatex.com/784817tfvqwp#/1721170/

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

reads the \*.vtk file produced by 'bin/bem_cost_function.cc'.

---

===============
Units
===============
To verify the units of the equation match let consider only the first term on the right hand side:

|unit_phi_of_r|

In the SI:

|unit_phi_of_r_si|

.. |unit_phi_of_r| image:: https://github.com/nsplab/sloc/blob/master/doc/unit_phi_of_r.png?raw=true 
.. |unit_phi_of_r_si| image:: https://github.com/nsplab/sloc/blob/master/doc/unit_phi_of_r_si.png?raw=true 

===============
Multiplying dipole magnitude by a constant 
===============
Let |phi| be the solution of the forward problem with dipole p at location r. 

|rtrue| and |ptrue| are the location and the magnitude of the dipole used in the
forward problem to simulate the potential measurements, |phitrue|.

You can estimate the magnitude of the dipole for the given set of true potential 
measurements and the true location of the dipole by |ptrueasterisk|.

If you multiply the magnitude of the dipole by a constant scalar value, c, 
|pprime|, you get a new set of potential measurements, |phiprime|. Then, you
can estimate the magnitude of the dipole for the given potential measurements,

|pasterisk|.

.. |phi| image:: https://github.com/nsplab/sloc/blob/master/doc/phi.png?raw=true 
.. |rtrue| image:: https://github.com/nsplab/sloc/blob/master/doc/rtrue.png?raw=true 
.. |ptrue| image:: https://github.com/nsplab/sloc/blob/master/doc/ptrue.png?raw=true 
.. |phitrue| image:: https://github.com/nsplab/sloc/blob/master/doc/phitrue.png?raw=true 
.. |ptrueasterisk| image:: https://github.com/nsplab/sloc/blob/master/doc/ptrueasterisk.png?raw=true 
.. |pprime| image:: https://github.com/nsplab/sloc/blob/master/doc/pprime.png?raw=true 
.. |phiprime| image:: https://github.com/nsplab/sloc/blob/master/doc/phiprime.png?raw=true 
.. |pasterisk| image:: https://github.com/nsplab/sloc/blob/master/doc/pasterisk.png?raw=true 

