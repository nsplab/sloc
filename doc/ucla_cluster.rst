
Getting Started
===============

To login into the Hoffman2 UCLA cluster::
    $ ssh hoffman2.idre.ucla.edu

Once you're logged in, you shouldn't run computationally intensive jobs
on the login node. Instead, obtain an interactive session with the `qrsh` command::
    $ qrsh -l i,mem=1G,time=4:00:00

For more details on interactive sessions, visit
http://www.ats.ucla.edu/clusters/hoffman2/computing/sge_qrsh.htm


Software Installation
=====================

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

Installing sloc
---------------

On a login node::
    $ cd ~/dev
    $ git clone gitosis@petaturtle.com:sloc.git

On an interactive session::
    $ cd ~/dev/sloc
    $ make all debug-mode=off


Running Programs
================

An explanation on how to run MPI programs on the Hoffman2 cluster can
be found at http://www.ats.ucla.edu/clusters/common/computing/parallel/using_mpi.htm


.. vim: ft=rst
