
Getting Started
===============

To login into the Hoffman2 UCLA cluster::
    $ ssh hoffman2.idre.ucla.edu

Once you're logged in, you shouldn't run computationally intensive jobs in the
login node on which you've logged in. Instead, obtain an interactive session
with the `qrsh` command::
    $ qrsh -l i,mem=1G,time=4:00:00

For more details on interactive sessions, visit
http://www.ats.ucla.edu/clusters/hoffman2/computing/sge_qrsh.htm


Setting Up Compiler And Other Software
======================================

To upgrade the version of the intel compiler to 12.1, and install other dependencies,
add these lines to the end of your `~/.bashrc` file, and log out of the cluster::
    module unload openmpi 2>/dev/null
    module unload intel 2>/dev/null
    module load intel/12.1
    module load openmpi/1.6
    module load python/2.7
    module load vtk
    module load cmake

After logging back in into the cluster, you can verify the above list of loaded modules
by issuing the following command::
    $ module list


Software Installation
=====================

Installing Boost
----------------

On a login node::
    $ cd ~/dev
    $ wget http://url/to/boost_1_49_0.tar.gz

On an interactive session::
    $ cd ~/dev
    $ tar xvfz boost_1_49_0.tar.gz
    $ cd boost_1_49_0
    $ ./bootstrap.sh 2>&1 | tee bootstrap.log
    $ ./b2 --help | less
    $ ./b2 --show-libraries
    $ vim user-config.jam  # add line "using mpi ;" to this file
    $ time ./b2 -j5 \
        --prefix=${HOME}/opt/local \
        --libdir=${HOME}/opt/local/lib \
        --layout=tagged \
        --user-config=user-config.jam \
        threading=multi \
        install \
        2>&1 | tee install.log

Installing deal.II (svn)
------------------------

On a login node::
    $ svn co http://www.dealii.org/svn/dealii/trunk/deal.II ~/dev/deal.II

Minor work-around for problem with deal.II finding the boost library directory::
    $ cd ~/opt/local/include
    $ ln -sf ../lib

On an interactive session::
    $ cd ~/opt/local/include
    $ cd ~/dev/deal.II
    $ time ./configure --enable-mpi --with-boost=${HOME}/opt/local/include 2>&1 | tee configure.log
    $ time make -j5 optimized 2>&1 | tee make.log

Installing muParser
-------------------

On a login node::
    $ svn co https://muparser.svn.sourceforge.net/svnroot/muparser/trunk/ muparser

On an interactive session::
    $ ./configure --disable-debug --disable-dependency-tracking --prefix=$HOME/opt/local 2>&1 | tee configure.log
    $ make 2>&1 | tee make.log
    $ make install 2>&1 | tee install.log
    $ pkg-config muparser --cflags --libs

Installing GetFEM++ (svn)
-------------------------

On a login node::
    $ svn co svn://svn.gna.org/svn/getfem/trunk/getfem getfem

On an interactive session::
    $ bash autogen.sh
    $ ./configure --with-pic 2>&1 | tee configure.log
    $ make 2>&1 | tee make.log

Installing sloc
---------------

On a login node::
    $ cd ~/dev
    $ git clone git@github.com:nsplab/sloc.git

On an interactive session::
    $ cd ~/dev/sloc
    $ make


Running Programs
================

An explanation on how to run MPI programs on the Hoffman2 cluster can
be found at http://www.ats.ucla.edu/clusters/common/computing/parallel/using_mpi.htm


.. vim: ft=rst
