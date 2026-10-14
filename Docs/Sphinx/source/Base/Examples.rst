.. _Chap:Examples:

Examples
========

In ``chombo-discharge``, applications are set up so that they use the ``chombo-discharge`` source code and one ``chombo-discharge`` physics module.
These are normally set up through Python interfaces accompanying each module. 
Several example applications are given in :file:`$DISCHARGE_HOME/Exec/Examples`, which are organized by example type (e.g., plasma simulation, electrostatics, radiative transfer, etc).
If ``chombo-discharge`` built successfully, it will usually be sufficient to compile the example by navigating to the folder containing the program file (:file:`main.cpp`) and compiling it:

.. code-block:: text

   make -s -j4 main

To see how these programs are run, see :ref:`Chap:Control`.   

Positive streamer in air
------------------------

To run one of the applications that use a particular ``chombo-discharge`` physics module, we will run a simulation of a positive streamer (in air). 

The application code is located in ``$DISCHARGE_HOME/Exec/Examples/CdrPlasma/DeterministicAir`` and it uses the convection-diffusion-reaction plasma module (located in ``$DISCHARGE_HOME/Physics/CdrPlasma``).

First, compile the application by

.. code-block:: text

   cd $DISCHARGE_HOME/Exec/Examples/CdrPlasma/DeterministicAir
   make -s -j4 DIM=2 main

This will provide an executable named ``main2d.<bunch_of_options>.ex``.
If one compiles for 3D, use ``DIM=3`` either on the command-line or in the configuration file.
The executable will be named ``main3d.<bunch_of_options>.ex``.

To run the application do:

**Serial build**

.. code-block:: text

   ./main2d.<bunch_of_options>.ex positive2d.inputs

**Parallel build**
  
.. code-block:: text

   mpirun -np 8 main2d.<bunch_of_options>.ex positive2d.inputs   

If the user also compiled with HDF5, plot files will appear in the subfolder ``plt``.

.. tip::

   One can track the simulation progress through the :file:`pout.*` files, see :ref:`Chap:pout`.  

