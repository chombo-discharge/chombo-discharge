.. _Chap:Design:

`PlasmaC` Software Design
===========================

This chapter discusses the `PlasmaC` software design. The flexibility of `PlasmaC` is enabled by segregating different parts into different modules, e.g. a time integrator class is responsible for advancing the equation set but has nothing to do with the actual physics that the equations describe. Likewise, the physics interface has no knowledge about the spatial discretization. In this chapter we discuss these modules, and how they can be combined.

.. _Chap:PlasmaCModules:

`PlasmaC` modules
-------------------

`PlasmaC` is a modular framework that allows reuse or reimplementation of many classes without affecting the underlying code. In `PlasmaC`, an entire simulation is run by an object called :ref:`Chap:plasma_engine` . The responsibility of this object is essentially to set up a simulation, run a simulation, and write output files. Although there will only be one instance of this object, there are no software restrictions in place that prevent users from creating more of them. If the user wants to run a full-blown plasma simulation, the mini-application system in `PlasmaC` requires the user to construct the :ref:`Chap:plasma_engine` object by supplying the remaining supplementary modules through its constructor. These modules (or classes), each describe geometry, physics, time integrator and so on. They are discussed below.

To facilitate code reuse, `PlasmaC` is set up such that a class :ref:`Chap:plasma_engine` takes as its input a number of modules. In this way, the user can re-use code (for example geometries or physics) over a range of applications. 

The various modules that go into :ref:`Chap:plasma_engine` are:

* :ref:`Chap:plasma_kinetics` An abstract class that defines the physics. This includes specifying the number of convected species, the number of RTE solvers, and well as specifiying how all these solvers are coupled (through e.g. source terms, boundary conditions and so on). 
* :ref:`Chap:computational_geometry` An implementation of the geometry that will be simulated. Various implementation of this class exist, which you may immediately use. New geometries are created by implementing a new :ref:`Chap:computational_geometry` class: This is a virtual class in `PlasmaC` that stores the various level set functions that describes the geometry. 
* :ref:`Chap:physical_domain` The physical domain to be simulated. This is a very lightweight class that only describes the axis-aligned box that holds your simulation domain. 
* :ref:`Chap:time_stepper` The temporal integrator; the class includes an implementation of the time stepping scheme. The base class is abstract and contains a number of useful routines that choreographs the coupling between solvers, for example by filling source terms in convection-diffusion-reaction solvers. Implementation of new integration schemes is time consuming, and additional implementations of this class is usually a developer task. `PlasmaC` comes with several high-performance time steppers that are useful. 
* :ref:`Chap:amr_mesh` The AMR mesh engine. This class includes grid generators, coarsening operators, ghost cell interpolation operators and so on. This class is one of the most important base class, and is responsible for orchestratic all spatial operations. Modifications to this class are (very) rarely made. 
* :ref:`Chap:cell_tagger` (Optional) Class that is responsible for refinement and coarsening decision. The base class is abstract, and users may implement their own classes if they like. Instance of this class will tell :ref:`Chap:plasma_engine` where to refine and coarsen cells. 
* :ref:`Chap:geo_coarsener` (Optional) A geometric coarsening class. This class is essentially a band-aid that allows users to remove mesh from uninteresting regions in the simulation domain which might otherwise have been tagged. 

Most users will only find the need to implement :ref:`Chap:plasma_kinetics`, :ref:`Chap:computational_geometry`, and possibly also :ref:`Chap:cell_tagger`. 

You will find a much more thorough explanation of these classes in the :ref:`Chap:ImportantClasses` chapter.

.. _Chap:MiniApplications:

Mini-applications
-----------------

In `PlasmaC`, simulation cases are created through a mini-application system. The user is responsible for compiling the executable (or mini-app), whose execution is controlled through an input script or through variables passed through the command line. In `PlasmaC`, the input script is read by using a Chombo class called ``ParmParse`` which read inputs from files or the command line. In `PlasmaC`, all input parameters are read in through the default constructor. In this way, all parameters are passed to their respective classes before the simulation begins. There is (currently) no support for changing input parameters during run-time. 


The mini-app executable is built by following the Chombo makefile system that tracks the dimensionality, compiler information etc. throughout your system. In reality, the C++ main file from which you will compile your executable is virtually identical across mini-applications: Users usually just replace geometries, integrators, kinetic schemes etc. Because of this, there is a python script supplied with the code that the user will find beneficial for setting up templated mini-apps. A generic setup for a mini-app looks something like this:

.. code-block:: c++

      #include "plasma_engine.H"   // Load plasma_engine class
      #include "rk2.H"             // Load temporal integrator instance (derived from time_stepper)
      #include "field_tagger.H"    // Load the cell tagger (derived from cell_tagger)
      #include "morrow_lowke.H"    // Load the plasma kinetics (derived from plasma_kinetics)
      #include "rod_sphere.H"      // Load the geometry (derived from computational_geometry)
      #include "geo_coarsener.H"   // Load the geometry grid coarsener

      #include <ParmParse.H>       // Input parameters parsing class. 

      Real g_potential;    
      Real potential_curve(const Real a_time){ // Potential curve to be simulated. This returns
        return g_potential;                    // a single value which is obtained through the input script. 
      }

      int main(int argc, char* argv[]){

      #ifdef CH_MPI
        MPI_Init(&argc,&argv);  // Initialize MPI
      #endif

        // Build argument list from input file and command line
        char* inputFile = argv[1];
        ParmParse PP(argc-2,argv+2,NULL,inputFile);
      
        { // Get the potential curve to be simulated (constant in this case)
          ParmParse pp("rod_sphere2d");
          pp.get("potential", g_potential);
        }

      	// Load the classes discussed above
        RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics> (new morrow_lowke());
        RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new rod_sphere());
        RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());
        RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(new rk2());
        RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new field_tagger());	
        RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
	RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<amr_mesh> (new geo_coarsener());
        RefCountedPtr<plasma_engine> engine            = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
		                                                                                         compgeom,
													 plaskin,
													 timestepper,
													 amr,
													 tagger,
													 geocoarsen));

      	// Give the potential curve the plasma_engine and run simulation
        engine->set_potential(potential_curve); // Provide potential curve to plasma_engine
        engine->setup_and_run();                // Run simulation
      
      
      #ifdef CH_MPI 
        MPI_Finalize(); // Finalize MPI
      #endif
      }

This is not much code. The first section of ``#include`` statements load the various `PlasmaC` modules, and the function that is defined outside ``main()`` defines the potential curve for the Poisson equation. It will be passed into ``plasma_engine`` which further distributes this function to other classes that might need it. The remaining pieces of code initializes MPI, reads the input script, and instantiates all the modules that are passed into ``plasma_engine``. Finally, ``plasma_engine`` is set up and run and MPI is finalized. In order to compile this code, you will also need a makefile that specifies how this will be compiled and linked against your Chombo library. To run the code, you will need an input script that contains all of the tunable parameters that controls your simulation. For most simulation cases, this script will contain several hundreds of parameters. Because the above steps are so similar across mini-applications, we have a Python script that automatically generates the setup of the above code, the required makefile, and a templated input file. This is discussed in the next section. 


.. _Chap:PythonInterface:

Python setup
------------

To simplify the setup of simulation cases, we've included a Python script that performs a templated setup of your case based on your selected modules. The Python script resides in source directory :file:`./` and is named :file:`setup.py`. To use it, you must pass the following variables through the command-line:

* ``CHOMBO_HOME`` (**optional**, defaults to ``$(CHOMBO_HOME)``. The path to your Chombo library, see :ref:`Chap:Environment` for details on how to set up your environment variables. 
* ``PLASMAC_HOME`` (**optional**, defaults to ``$(PLASMAC_HOME)``. The path to your `PlasmaC` library
* ``DIM`` (**optional**, defaults to 2). The problem dimensionality, which can be 2 or 3. 
* ``base_dir``. The directory in which your application will be placed
* ``app_name``. The name of your mini app. Your code will be placed in :file:`base_dir/app_name`.
* ``file_name`` (**optional**, defaults ``main``).
* ``plasma_kinetics``. Your :ref:`Chap:plasma_kinetics` implementation. `PlasmaC` will look for this (and an option file) in :file:`./plasma_models/<your_kinetics>`. See :ref:`Chap:Directories` for details. 
* ``geometry`` (**optional**, defaults to ``regular_geometry``). Your geometry. `PlasmaC` will look for this (and an option file) in :file:`./geometries_prebuilt`. See :ref:`Chap:Directories` for details.
* ``time_stepper`` The temporal integrator. If you write your own, it should reside in the :file:`./src/time_steppers/<your_time_stepper>` directory. See :ref:`Chap:Directories` for details. 
* ``cell_tagger`` (**optional**, defaults ``NULL``). Your :ref:`Chap:cell_tagger` implementation. The tagger you provide should reside in the :file:`./src/cell_taggers/<my_tagger>` directory.

The Python interface automates the setup of a main-file through which you can compile your application, and also provides a makefile for compilation. The makefile expects that the source code for your modules reside in the folders listed above. In addition to this, the Python interface will expect a file which holds the all the tunable input variables associated with a class. For example, ``amr_mesh`` contains a large number of variables that control grid generation, all of which are stored in :file:`/src/amr_mesh.options`. 

To get help with the Python interface, you can do

.. code-block:: bash

   ./setup.py -h

This will list the input arguments that you must provide.

Using the Python script is very simple:

.. code-block:: bash

   ./setup.py -base_dir=mini_applications -app_name=my_application -plasma_kinetics=my_kinetics


There are also options for direct building of your application. To do this, you must pass additionally pass ``-build=true``. You may also select the number of processes used for building and turn off compiler outputs. For example:
   
.. code-block:: bash

   ./setup.py -base_dir=mini_apps -app_name=my_app -plasma_kinetics=my_kinetics -build=true -silent=true -procs=10

.. _Chap:CodeStructure:

Code Structure
--------------

Here, we provide an overview of the `PlasmaC` directories and coding styles.

.. _Chap:Directories:

Directories
___________

The following directories in `PlasmaC` are worth noting:

* :file:`/src` contains the `PlasmaC` source code discussed in :ref:`Chap:ImportantClasses`. 
 
  * :file:`/src/amr_mesh` contains :ref:`Chap:amr_mesh` related code
  * :file:`/src/cdr_solver` contains code for the CDR solvers
  * :file:`/src/elliptic` contains operators for elliptic equations (mostly multifluid Poisson stuff)
  * :file:`/src/geometry` contains code related to the geometric interface
  * :file:`/src/global` contains some globally useful code, such as data structures, stencil types and so on.
  * :file:`/src/plasma_solver` contains the plasma framework, i.e. :ref:`Chap:plasma_kinetics`, :ref:`Chap:plasma_engine` and some related code.
  * :file:`/src/poisson_solver` contains the abstract Poisson solver class and it's geometric multigrid implementation.
  * :file:`/src/rte_solver` contains the RTE solvers
  * :file:`/src/sigma_solver` contains the surface charge solver
* :file:`/geometries_prebuilt` contains some predefined geometries.
* :file:`/plasma_models` and its subdirectories contains various implementation of :ref:`Chap:plasma_kinetics`. 
* :file:`/cell_taggers` and its subdirectories contains various implementation of :ref:`Chap:cell_tagger`.
* :file:`/time_steppers` and its subdirectories contains various implementation of :ref:`Chap:time_stepper`.
* :file:`/base_tests` contains some base tests of `PlasmaC`
* :file:`/doc` contains the documentation of `PlasmaC`
    
  * :file:`/doc/sphinx` contains the Sphinx documentation
  * :file:`/doc/doxygen` contains some markup used for the :doxy:`Doxygen API <index>`.
  * :file:`/doc/figures` contains some figures used throughout the documentation. 
* :file:`/app_builder` contains the Python interface for setting up mini-applications.


If you want to extend the `PlasmaC` code, you *may* write your own mini-apps outside of the `PlasmaC` framework. However, for maximum reuseability you might want to ensure that your changes are available in the future as well. We recommend that you place your geometries, plasma kinetics, and cell taggers in the appropriate directories listed above. This will also ensure that your work can be reached through our :ref:`Chap:PythonInterface`.

.. _Chap:InputVariables:

Input variables
_______________

Generally, the coding style for input variables is to use the class name as a prefix (where :ref:`Chap:amr_mesh` is an exception) and the variable as a suffix. All letters are lower-case. For example::

   plasma_engine.max_steps = 10

To pass input variables into `PlasmaC`, we generally refrain from hard-coding variables that should be accessible to the user. Instead, we use Chombo's ParmParse class, which is used in the following way:

.. code-block:: c++

   Real my_variable;
   ParmParse pp("prefix");
   pp.get("suffix", my_variable);

The above code segment will try to fetch an input line ``prefix.suffix`` and place it in *my_variable*. Note that the specification of ``prefix.suffix`` should be of the same type as ``my_variable`` (float in this case). For this example, passing

.. code-block:: bash

		mpirun -np 32 <my_application> <my_input_file> prefix.suffix = foo

will throw an error. There are, of course, many input parameteres that the user will want to tune when he runs a simulation. You will find a compiled list of all tunable parameters in the detailed discussion of the implementation classes in the :ref:`Chap:ImportantClasses` chapter. 

.. _Chap:Chombo:

Chombo coding guide
___________________

`PlasmaC` is mostly a large `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_ application. Chombo uses dimension-independent data structures. Since these structures are used in the physics interface, the user should familiarize himself with them. The most important structures are

* :file:`Real` - a replacement for float or double (depending on your compiler settings)
* :file:`RealVect` - a vector in space.
* :file:`Vector` - a wrapper for :file:`std::vector`.
* :file:`RefCountedPtr<T>` - a pointer class with reference counting and auto-deallocation.

The useage of these classes is straightforward. For example, a :file:`Real` is declared

.. code-block:: c++

		Real foo = 1.0;

:file:`RealVect` is a spatial vector that contains two or three entries in `PlasmaC`. To use :file:`RealVect`, one may do

.. code-block:: c++

		RealVect foo = RealVect(1.0, 0.0);

in two dimensions and

.. code-block:: c++

		RealVect foo = RealVect(1.0, 0.0, 0.0);

in three dimensions. The dimensionless way of doing this is to use Chombo macros; 

.. code-block:: c++

		RealVect foo = RealVect(D_DECL(1.0, 0.0, 0.0));

where :file:`D_DECL` is macro that returns the first two variables in 2D, and all three variables in 3D.

The :file:`Vector` class is used just as :file:`std::vector`: 		

.. code-block:: c++

   Vector<Real> foo(2);
   foo[0] = 1.0;
   foo[1] = 0.0;
		
The same goes with the smart pointer :file:`RefCountedPtr<T>`:
   
.. code-block:: c++

   RefCountedPtr<Real> ptr = RefCountedPtr<Real> (new Real(0.0));

For the full Chombo API, please see the `Chombo doxygen guide <http://davis.lbl.gov/Manuals/CHOMBO-RELEASE-3.2/classes.html>`_. 
