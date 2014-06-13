==========================
Source Localization (sloc)
==========================

Installation
============

Dependencies
------------

Boost
  (Version 1.54.0.1ubuntu1 from apt-get)
  libboost-dev

muparser
  http://muparser.beltoforion.de/ Version 2.2.3

GetFEM++
  http://download.gna.org/getfem/stable/getfem-4.2.tar.gz Version 4.2
  May need to run
  >> export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}/usr/local/lib
  
Deal.II Version 8.1.0
  http://www.dealii.org/
  http://www.ces.clemson.edu/dealii/deal.II-8.1.0.tar.gz
  cmake -DDEAL_II_WITH_MPI=ON .

Thread Building Blocks (tbb) [check if this generates an error in cmake when missing]
  libtbb-dev (4.2~20130725-1.1ubuntu1)

MPI
  apt-get install libopenmpi-dev
  libmpich-dev seems to be another option





Compiling
---------

g++ Version 4.8.3

Running
-------

Running code
add `pwd` to path
run icosahedron.run.sh from data/simple



Structure
=========

|block_diagram|

.. |block_diagram| image:: https://github.com/nsplab/sloc/blob/master/doc/block_diag.png?raw=true 


Programs
========


File Types
==========

.mesh file
  this is a standard file format, not custom made by us.
  for geometric meshes (tetrahedron) created by hand; for sphere, stl file was created in blender, and then stl_to_mesh was used to convert from .stl to .mesh.
  Not sure exactly how icosahedron was created (mosalam, 6/12/14 @ 7:38pm)
  for anatomical meshes, .stl file created in Vitrea was converted to .mesh using make_head_model

.sigma file
  custom file format made by us; 
  lookup table for the material index
  each “material index” is actually an index for an interface between two materials with same or different conductivities.
  one row for each material index - eg. with a 7 layer model, there are 7 material indices
  first column: material index
  second column: inner conductivity
  third column: outer conductivity

.mat file
  custom file format made by us; 
  each mesh has a separate .mat file
  .mat stands for “material”
  the class sloc::MaterialData defined in material_data.cc creates the .mat and .sigma files.
  first line: # of triangles in mesh
  second line onwards (one row for each triangles)
	first column: index of vertex
	second column: material index (also called material number)

.dat file
  custom file format made by Luis
  potentials on electrodes
  first line: number of electrodes N
  next N lines: electrode index (integer) and potential on electrode (floating-point)

.stl file
  format from Vitrea

.vtk file
  Format from the Visualization Toolkit
  contains information equivalent to .dat files, but in a format supported by the Visualization Toolkit

.prm file
  format from deal.II Parameter Handler
  specifies parameters for a function as a text file

.log file
  custom format
  debugging output
  not essential for later use

.cost_at_grid_pts
  created by mosalam with this program:  bem_cost_function , which reads potentials from two .dat files
  called within this shell script : head.run_grid.sh
  content: the cost computed by putting the candidate dipole source at the grid points
  stores the cost of best dipole fit at each candidate location in a grid of candidate points.
  the cost is the sum of squared errors between predicted and “measured” potentials (?)
  7 columns: (x, y, z, angleX, angleY, angleZ, cost)
  one row for each candidate point


.electrodes
  created by select_electrodes_given_3d_pos (see head.run_grid.sh that calls this)
  the vertex indices of the electrodes (10-20 system)
  takes .stl vertices for scalp and points for 10-20 electrode configuration and gives the vertices closest to the true locations

example
  head.mesh - contains nearly 16,000 triangles.
  head.mat - contains the material information for each triangle, including the material index for the inside and the outside of each triangle
  head.sigma - contains a lookup table that relates the material index to the inner and outer conductivity

Alternative Methods
===================
Need to compare results with

eeglab - NFT (directory: mfiles) version 2.3 (Mosalam)
          includes forward problem solutions
Field Trip

README for sloc
===============

For installation instructions refer to ``doc/hoffman2.rst``.

**Manuscript Drafts**

###Endovascular Source Localization (simulation study)

https://www.writelatex.com/784824tywgtg#/1721178/

###Derivation of boundary element method (BEM) equation used in forward model, based on Luis' notes

https://www.writelatex.com/784817tfvqwp#/1721170/

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

Units
===============
To verify the units of the equation match let consider only the first term on the right hand side:

|unit_phi_of_r|

In the SI:

|unit_phi_of_r_si|

.. |unit_phi_of_r| image:: https://github.com/nsplab/sloc/blob/master/doc/unit_phi_of_r.png?raw=true 
.. |unit_phi_of_r_si| image:: https://github.com/nsplab/sloc/blob/master/doc/unit_phi_of_r_si.png?raw=true 

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

