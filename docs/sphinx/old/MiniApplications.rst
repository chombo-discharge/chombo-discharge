.. _Chap:MiniApplications:

Mini-applications
-----------------

In PlasmaC, simulation cases are created through a mini-application system. The user is responsible for compiling the executable that controls the behavior of the mini-app, which is controlled through an input script. In PlasmaC, the input script is read by using a Chombo class called ParmParse which read inputs from files or the command line. An important feature of PlasmaC is that all input parameters are read in through the default constructor. In this way, all parameters are passed to their respective classes before the simulation even begins. We strongly encourage that the user follows this system: see :ref:`Chap:CodeStructure` for coding styles and tips. 

When parameters are passed to the PlasmaC mini-app, the nomenclature that we use is to use the class name as a prefix and the tunable input parameter as a suffix. The two are separated by exactly one period. For example, there is an option in plasma_engine that controls the maximum number of integration steps during a simulation, which is defined like this:

.. code-block:: c++
		
		plasma_engine.max_steps = 10

This will specify that the plasma_engine class, which is the class that essentially runs a plasma simulation, will perform a maximum of 10 time steps. Chombo's ParmParse can read a number of different natives (floats, integers, strings etc.), but they are always in the form above. 

There are, of course, many input parameteres that the user will want to tune when he runs a simulation. You will find a compiled list of all tunable parameters in the detailed discussion of the implementation classes in the :ref:`Chap:ImportantClasses` chapter. Following the guidelines in the :ref:`Chap:CodeStructure` chapter of this documentation, you will also know that there should be an options file accompanying every class with tunable input parameters. 

The mini-app executable is built by following the Chombo makefile system that tracks the dimensionality, compiler information etc. throughout your system. In reality, the C++ main file from which you will compile your executable is essentially identical across mini-apps: Users essentially just replace geometries, integrators, kinetic schemes etc. Because of this, there is a python script supplied with the code that the user will find beneficial for setting up templated mini-apps. In fact, the generic setup for a mini-app is like this:

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

As far as code goes, this is not much. However, you will also need a makefile that specifies how this executable will be compiled and linked against your Chombo library, as well as an input script that contains all of the tunable parameters that controls your simulation. For most simulation cases, this script will contain several hundred lines of parameters, and it is impossible to keep track of these. Because these things are so similar across mini-applications, we have a Python script that automates this setup for you. See :ref:`Chap:PythonInterface` for details. 
