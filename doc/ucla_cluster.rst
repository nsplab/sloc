
Getting Started
===============

To login into the Hoffman2 UCLA cluster::
    $ ssh hoffman2.idre.ucla.edu

Once you're logged in, you shouldn't run computationally intensive jobs
on the login node. Instead, obtain an interactive session with the `qrsh` command::
    $ qrsh -l i,mem=1G,time=4:00:00

For more details on interactive sessions, visit
http://www.ats.ucla.edu/clusters/hoffman2/computing/sge_qrsh.htm

Setting Up Compilers
====================

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
    $ svn co http://www.dealii.org/svn/dealii/trunk/deal.II ~/dev/deal.II-trunk

Minor work-around for problem with deal.II finding the boost library directory::
    $ cd ~/opt/local/include
    $ ln -sf ../lib

On an interactive session::
    $ cd ~/opt/local/include
    $ cd ~/dev/deal.II-trunk
    $ time ./configure --enable-mpi --with-boost=${HOME}/opt/local/include 2>&1 | tee configure.log
    $ time make -j5 optimized 2>&1 | tee make.log

Installing deal.II
------------------

On a login node::
    $ cd ~/dev
    $ wget http://www.dealii.org/download/deal.II-7.1.0.tar.gz

On an interactive session (note the `--enable-mpi` flag)::
    $ cd ~/dev
    $ tar xvfz deal.II-7.1.0.tar.gz
    $ cd deal.II
    $ time ./configure --enable-mpi 2>&1 | tee configure.log
    $ time make -j2 optimized 2>&1 | tee make.log

Installing VTK
--------------

VTK 5.8.0 is already installed on the cluster. It can be loaded
using the modules system::
    $ module avail
    $ module show vtk
    $ module load vtk
    $ module list

Installing sloc
---------------

On a login node::
    $ cd ~/dev
    $ git clone gitosis@petaturtle.com:sloc.git

On an interactive session::
    $ cd ~/dev/sloc
    $ make all debug-mode=off

Installing zeromq
-----------------

On a login node::
    $ cd ~/dev
    $ wget http://download.zeromq.org/zeromq-3.2.0-rc1.tar.gz

On an interactive session::
    $ cd ~/dev
    $ tar xvfz zeromq-3.2.0-rc1.tar.gz
    $ cd zeromq-3.2.0
    $ ./configure --prefix=$HOME/opt/local 2>&1 | tee configure.log
    $ make 2>&1 | tee make.log
    $ make install 2>&1 | tee install.log

To install the C++ bindings, do this on a login node::
    $ cd ~/dev
    $ git clone git://github.com/zeromq/cppzmq.git
    $ cd cppzmq
    $ ln -s $PWD/zmq.hpp ~/opt/local/include

Installing msgpack
------------------

On a login node::
    $ cd ~/dev
    $ wget http://msgpack.org/releases/cpp/msgpack-0.5.7.tar.gz

On an interactive session::
    $ cd ~/dev
    $ tar xvfz msgpack-0.5.7.tar.gz
    $ cd msgpack-0.5.7
    $ ./configure --prefix=$HOME/opt/local 2>&1 | tee configure.log
    $ make 2>&1 | tee make.log
    $ make install 2>&1 | tee install.log

Running Programs
================

An explanation on how to run MPI programs on the Hoffman2 cluster can
be found at http://www.ats.ucla.edu/clusters/common/computing/parallel/using_mpi.htm


.. vim: ft=rst
