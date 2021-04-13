.. _Chap:PythonInterface:

Python interface
----------------

To simplify the setup of simulation cases, we've included a Python script that performs a templated setup of your case based on your selected modules. The Python script resides in source directory :file:`./` and is named :file:`build.py`. To use it, you must pass the following variables through the command-line:

* ``CHOMBO_HOME`` (**optional**, defaults to ``$(CHOMBO_HOME)``. The path to your Chombo library, see :ref:`Chap:Environment` for details on how to set up your environment variables. 
* ``PLASMAC_HOME`` (**optional**, defaults to ``$(PLASMAC_HOME)``. The path to your PlasmaC library
* ``DIM`` (**optional**, defaults to 2). The problem dimensionality, which can be 2 or 3. 
* ``base_dir``. The directory in which your application will be placed
* ``app_name``. The name of your mini app. Your code will be placed in :file:`base_dir/app_name`.
* ``file_name`` (**optional**, defaults ``main``).
* ``plasma_kinetics``. Your :ref:`Chap:plasma_kinetics` implementation. PlasmaC will look for this (and an option file) in :file:`./plasma_models/<your_kinetics>`. See :ref:`Chap:Directories` for details. 
* ``geometry`` (**optional**, defaults to ``regular_geometry``). Your geometry. PlasmaC will look for this (and an option file) in :file:`./geometries_prebuilt`. See :ref:`Chap:Directories` for details.
* ``time_stepper`` (**optional**, defaults ``rk2``). Your integrator - this defaults to a second order Runge-Kutta method. If you provide your own, it should reside in the :file:`./src/time_steppers/<your_time_stepper>` directory. See :ref:`Chap:Directories` for details. 
* ``cell_tagger`` (**optional**, defaults ``NULL``). Your :ref:`Chap:cell_tagger` implementation. The tagger you provide should reside in the :file:`./src/cell_taggers/<my_tagger>` directory.

The Python interface automates the setup of a main-file through which you can compile your application, and also provides a makefile for compilation. The makefile expects that the source code for your modules reside in the folders listed above.

To get help with the Python interface, you can do

.. code-block:: bash

   ./build.py -h

This will list the input arguments that you must provide.

Using the Python script is very simple:

.. code-block:: bash

   ./build.py -base_dir=mini_applications -app_name=my_application -plasma_kinetics=my_kinetics


There are also options for direct building of your application. To do this, you must pass additionally pass ``-build=true``. You may also select the number of processes used for building and turn off compiler outputs. For example:
   
.. code-block:: bash

   ./build.py -base_dir=mini_applications -app_name=my_application -plasma_kinetics=my_kinetics -build=true -silent=true -procs=10
