.. _Chap:RunningSimulations:

Running simulations
===================

To run simulations, the user must first compile his application. Once the application has been defined (either from :ref:`Chap:WorkedExample1` or by using the :ref:`Chap:PythonInterface`, he may compile is by

.. code-block:: bash

   make -s -j 32 DIM=N <application_name>

where *N* may be 2 or 3, and <application_name> is the name of the mini-application. This will compile an executable whose name depends on your application name and compiler settings. Please refer to the Chombo manual for explanation of the executable name. You may, of course, rename your application.

Next, applications are run by

.. code-block:: bash

		mpirun -np 32 <application_executable> <input_file>

where <input_file> is your input file. On clusters, this is a little bit different and usually requires passing the above command through a batch system.

You may also pass input parameters through the command line. For example, running

.. code-block:: bash

   mpirun -np 32 <application_executable> <input_file> plasma_engine.max_steps=10

set the ``plasma_engine.max_steps`` parameter to 10. Note that this overrides whatever is defined in the <input_file>. 

.. _Chap:RestartingSimulations:

Restarting simulations
----------------------

Restarting simulations is done in exactly the same way as running simulations, although the user must set the ``plasma_engine.restart`` flag and the ``plasma_engine.restart_step`` input variable. For example:

.. code-block:: bash

   mpirun -np 32 <application_executable> <input_file> plasma_engine.restart=true plasma_engine.restart_step=10

When a simulation is restarted, PlasmaC will look for a checkpoint file with the ``plasma_engine.output_names`` variable and the specified restart step. If this file is not found, restarting will not work. You must therefore ensure that your executable can locate this file. This also implies that you cannot change the ``plasma_engine.output_names`` variable during restarts, unless you also change the name of your checkpoint file.

Changing your physics
_____________________

During the restart step, PlasmaC will load the initial grids and checkpointed data into memory. This data resides in an HDF5 file with where appropriate headers are used to identify where the data belongs. Amongst other things, the names of these headers are taken from :ref:`Chap:plasma_kinetics`, so you cannot change the species during during restarts. Currently, PlasmaC requires the exact same number of species during restarts, as well as consistent names for these. However, you *may* change the :ref:`Chap:plasma_kinetics` core functions, allowing you to change your plasma chemistry during restarts.

Changing your integrator
________________________

Users may also change temporal integrators during restarts. However, this requires you to 1) in advance set up two different integrators and provide input variables to these or 2) recompile your executable with a different integrator.

Changing spatial discretization
_______________________________

Spatial discretization may be changed during restarts. **However, you are *not* allowed to change the geometry or physical domain.** Furthermore, the following :ref:`Chap:amr_mesh` input variables are off-limits:

* ``amr.coarsest_domain``
* ``amr.max_amr_depth``
* ``amr.ref_rat``

If you change these variables, the checkpointed data cannot be imported into memory. In principle, we *can* extend PlasmaC so that this will be allowed. If you really, really want this feature, please :ref:`Chap:contact`.

Note that whatever changes you otherwise apply to :ref:`Chap:amr_mesh` become active only after the first regrid. 

Changing other settings
_______________________

Apart from the above variables, most changes are allowed during restarts. For example, you are allowed to use different tagging criteria (or even entirely different tagging classes); you can change the solver settings or applied potential; alter the output routines, and so on.

For example, here is a code snippet (see :ref:`Chap:MiniApplications` for the full code) that allows you to change your cell tagger during restarts

.. code-block:: c++
	  
   ParmParse pp("my_application");
   bool use_my_tagger = false;
   pp.query("change_tagger", use_my_tagger);

   RefCountedPtr<cell_tagger> tagger;
   if(use_my_tagger){
      tagger = RefCountedPtr<cell_tagger> (new my_tagger());
   }
   else{
      tagger = RefCountedPtr<cell_tagger> (new field_tagger());
   }

   RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
   RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<amr_mesh> (new geo_coarsener());
   RefCountedPtr<plasma_engine> engine            = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
		                                                                                    compgeom,
												    plaskin,
												    timestepper,
												    amr,
												    tagger,
												    geocoarsen));

In the above, we assume that *my_tagger* and *field_tagger* are separate implementations of :ref:`Chap:cell_tagger`, and we have created an input variable ``my_application.change_tagger`` which allows for specification of the cell tagger at run time. 
