.. _Chap:Prerequisites:

Prerequisites
-------------

From the ground up, PlasmaC is built on top of the `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_ framework. To compile PlasmaC, you must therefore have the following in place:

* A Fortran compiler, typically gfortran or Intel Fortran
* A C++ compiler, typically g++ or Intel C++
* An MPI installation
* A parallel HDF5 installation
* A Chombo library

Typically, local clients (laptops and desktops) already have appropriate Fortran and C++ compilers installed, as well as a version of MPI. On clusters, HDF5 is also preinstalled (usually), and in this case, it will be sufficient to modify the Chombo build files in order to compile PlasmaC. If you already have HDF5 installed, you may skip directly to :ref:`Chap:Environment`.

.. _Chap:HDF5:

Installing HDF5
_______________

If you do not have HDF5 installed, you may do the following:

1. Compile and install zlib, which is a compression library used by HDF5. zlib can be installed by
   
   .. code-block:: c++
		
		sudo apt-get install zlib1g-dev

2. Download HDF5 (version 1.8 or newer) and install it for parallel execution

      .. code-block:: c++
		
		      ./configure --prefix=/usr/local/hdf5 --enable-production --enable-fortran --enable-parallel
		      make
		      make install

   This will install HDF5 in /usr/local/hdf5. 
