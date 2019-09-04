.. _Chap:NewSimulations:

Setting up simulations
======================

In this section we consider the setup of new mini-applications. All these examples use the Python interface for setting up problems. The complexity of these examples progressively increase. Some of the examples at the end of this chapter also include some actual C++ coding for implementation of new physics and geometries.

Advecting a scalar
__________________

This example considers the advection of a scalar in a rod geometry. We consider the example in :ref:`Chap:MyFirstCompilation`, with the addition of a geometry. We use a pre-defined geometry called `rod_sphere` which is implemented in :file:`geometries_prebuilt/rod_sphere`. In 3D, this geometry consists of rod electrode and a dielectric sphere whereas (in 2D the geometry is interpreted as a blade and a cylinder). To set up the problem, we use the Python interface

.. code-block:: bash

  ./setup.py -base_dir=./mini_apps -app_name=advection2d -plasma_kinetics=advection_kinetics -geometry=rod_sphere

The Python parser will now generate the compilation file and also fetch all the required geometry files. The parameters in :file:`geometries_prebuilt/rod_sphere.options` are

.. code-block:: bash

   # ====================================================================================================
   # ROD_SPHERE CLASS OPTIONS
   # ====================================================================================================
   rod_sphere.eps0                      = 1                                # Background permittivity
   rod_sphere.turn_off_electrode        = false                            # Turn on/off electrode
   rod_sphere.turn_off_dielectric       = false                            # Turn on/off dielectric
   rod_sphere.electrode_live            = true                             # Live electrode or not
   rod_sphere.electrode_radius          = 5.E-3                            # Electrode inner radius
   rod_sphere.electrode_center1         = 0.0 0.0 1E-2                     # Center 1
   rod_sphere.electrode_center2         = 0.0 0.0 0E-2                     # Center 2
   rod_sphere.dielectric_permittivity   = 4.0                              # Dielectric permittivity
   rod_sphere.dielectric_center         = 0.0 0.0 0.0                      # Dielectric center
   rod_sphere.dielectric_radius         = 1.0                              # Dielectric radius

The explanation of the parameters should be fairly obvious. All control parameters should *always* be lower-case, and boolean values are passed in through ``false`` or ``true``, Parameters that contain more than one entry are usually vector entries. In this case, the parameter ``rod_sphere.electrode_center1`` indicates the position of one of the endpoints of the rod in Cartesian coordinates. Likewise, the parameter ``rod_sphere.electrode_radius`` is the radius of that rod.

   We will compile this example in 3D:

.. code-block:: bash

   make -s -j8 DIM=3 main

Next, we will create an electrode with a radius of 0.1 which is centered in the middle of the domain. We change the following parameters

.. code-block:: bash

   rod_sphere.turn_off_dielectric       = true                             # Turn on/off dielectric
   rod_sphere.electrode_radius          = 0.1                              # Electrode inner radius
   rod_sphere.electrode_center1         = 0.0 0.0 10.0                     # Center 1
   rod_sphere.electrode_center2         = 0.0 0.0 0.0                      # Center 2

You may change these either in the input file, or on the command line (see :ref:`Chap:Control`). We then run the example as usual

.. code-block:: bash

   mpirun -np8 main3d.<bunch_of_options>.ex template.inputs
