.. _Chap:WorkedExample:

Worked Example
--------------

Once you have set up PlasmaC, it is much easier to use. The coding style in PlasmaC allows user to create simulation applications from existing modules. In this section we show to to build mini-apps using a Python interface that autogenerates the missing implementation pieces, sets up a simulation directory, and a makefile.

We will simulate a positive streamer discharge from a needle in a :math:`(1\textrm{cm})^2` domain. In this particular case, we will use the following

Setting up the simulation case
______________________________

1. The Morrow-Lowke plasma model for air, which resides in $(PLASMAC_HOME)/plasma_models
2. A prebuilt geometry describing a cylindrical rod for electrode, which resides in $(PLASMAC_HOME)/geometries_prebuilt
3. A second order Runge-Kutta time stepping scheme, which resides in $(PLASMAC_HOME)/time_steppers
4. A cell tagger that refines mesh based on field magnitudes and curvatures, which resides in $(PLASMAC_HOME)/cell_taggers

All these pieces have already been compiled, and reside in various directories in the source code. See :ref:`Chap:CodeStructure` for details.

To build this app, we assume that the user has already set the CHOMBO_HOME and PLASMAC_HOME environment variables. Next, navigate to PLASMAC_HOME and execute

Some body text explaining what's coming up::

    python build.py -dim=2 -base_dir=my_applications -app_name=worked_example -plasma_kinetics=morrow_lowke -geometry=rod_sphere -time_stepper=rk2 -cell_tagger=field_tagger

   
This little Python script is responsible for setting up mini-app. The flags that are passed into the script specify the dimension, the directory of the app, the application name, the plasma kinetics, the geometry, the temporal integrator, and the name of the cell tagger. For more information on the Python script, please see :ref:`Chap:PythonInterface`. 


Compiling the mini-app
______________________

The compilation of the mini-app is normally straightforward: You simply navigate to the location of wherever you installed your mini-app and compile::

  > make -j 10 DIM=2 main

If everything is set up correctly, this will compile your application. Here, we've specified the dimensionality of the application, which is two-dimensional.

If your application did not compile, something has gone wrong. Typically, there are a ton of things that can go wrong during the compilation stage; for example missing compilers, no HDF5 library (or wrong path specified). To troubleshoot miscompilations, you should make sure that

1. Chombo has a well-defined Make.defs.local file. See :ref:`Chap:GettingStarted`. 
2. HDF5 is installed with MPI
3. Your PLASMAC_HOME and CHOMBO_HOME variables are set up correctly.

Usually, once these three things are asserted, compilation is normally successful. 


Defining the simulation case
____________________________

When one uses the Python script for setting up simulation cases, all input parameters for every class that is a part of the simulation case are copied to a template file called template.inputs. This file contains all the input variables that the developer has made available for your implementation classes. Usually, there are hundreds of parameters that can be tuned. 

Note that the parameters in template.inputs are default parameters, and you will need to verify them for your simulation case. Every banner that looks like this::

  # =======================
  # SOME_TEXT CLASS OPTIONS
  # =======================

indicates that the parameters that follow belong the the SOME_TEXT class. In our case, we will only modify a few of these parameters. For a detailed explanation of all of them, see the :ref:`Chap:ImportantClasses` chapter.

For our simulation, we will only modify the computational domain, the geometry, the AMR settings, and the number of time steps that we will run. First, we specify the geometry::

  # ========================
  # ROD_SPHERE CLASS OPTIONS
  # ========================
  rod_sphere.eps0                      = 1                                # Background permittivity
  rod_sphere.turn_off_electrode        = false                            # Turn on/off electrode
  rod_sphere.turn_off_dielectric       = true                             # Turn on/off dielectric
  rod_sphere.electrode_live            = true                             # Live electrode or not
  rod_sphere.electrode_radius          = 1.E-3                            # Electrode inner radius
  rod_sphere.electrode_center1         = 0.0 1.02                         # Center 1
  rod_sphere.electrode_center2         = 0.0 0.0    2                     # Center 2
  rod_sphere.dielectric_permittivity   = 4.0                              # Dielectric permittivity
  rod_sphere.dielectric_center         = 0.0 0.0 0.0                      # Dielectric center
  rod_sphere.dielectric_radius         = 1.0                              # Dielectric radius


  
Running the simulation
______________________

Not written yet
