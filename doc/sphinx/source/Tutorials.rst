.. _Chap:Tutorials:

Tutorials
=========

In this section we consider the setup of new mini-applications. All these examples use the Python interface for setting up problems. The complexity of these examples progressively increase. Some of the examples at the end of this chapter also include some actual C++ coding for implementation of new physics and geometries.

Scalar advection
----------------

Getting a geometry
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

We will compile this example in 2D:

.. code-block:: bash

   make -s -j8 DIM=2 main

Next, we will create an electrode with a radius of 0.1 which is centered in the middle of the domain. We change the following parameters

.. code-block:: bash

   rod_sphere.turn_off_dielectric       = true                             # Turn on/off dielectric
   rod_sphere.electrode_radius          = 0.1                              # Electrode inner radius
   rod_sphere.electrode_center1         = 0.0 10.0                         # Center 1
   rod_sphere.electrode_center2         = 0.0  0.0                         # Center 2

You may change these either in the input file, or on the command line (see :ref:`Chap:Control`). We then run the example as usual

.. code-block:: bash

   mpirun -np 8 main2d.<bunch_of_options>.ex template.inputs

If you want to run this as a 3D example, you can compile with DIM=3 and ``rod_sphere.electrode_center1 = 0.0 0.0 10.0``.

Refining the geometry
_____________________

Next, we show to to refine the geometry. In the :ref:`Chap:amr_mesh` class, the flag ``amr_mesh.max_amr_depth`` controls the maximum AMR depth that can be used for a simulation. The default in `PlasmaC` is to refine geometries down to the finest AMR level, but this *can* be changed through through :ref:`Chap:plasma_engine`. We will run the example again with ``amr_mesh.max_amr_depth = 2``, which provides us with two levels of mesh refinement:

.. code-block:: bash

   mpirun -np 8 main2d.<bunch_of_options>.ex template.inputs amr_mesh.max_amr_depth = 2. 

We see that the entire geometry was refined. If we only want to refine the tip of the geometry, we could use the refinement class :ref:`Chap:cell_tagger`, or we could use a lighter option provided through :ref:`Chap:geo_coarsener`. We will use the latter in order to remove some of the tags. In :file:`template.inputs` We set

.. code-block:: bash

   geo_coarsener.num_boxes   = 1            # Number of coarsening boxes (0 = don't coarsen)
   geo_coarsener.box1_lo     = -1.0 0.2     # Remove irregular cell tags 
   geo_coarsener.box1_hi     =  1.0 1.0     # between these two corners
   geo_coarsener.box1_lvl    = 1            # up to this level

This will remove all *geometric* refinement flags in the region defined by ``box1_lo`` and ``box1_hi``. If you want to remove other tags, you may add more boxes with the syntax ``box2_lo``, ``box2_hi`` and so on. 
