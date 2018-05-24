Program description {#program-structure}
==============================

chombo-streamer is a program for performing Cartesian 2D and 3D transient simulations of fluid plasmas in complex geometries. The code is flexible with regards to geometries and plasma kinetics, and is thus applicable to a comparatively broad range of research fields. chombo-streamer is built on top of the Chombo platform, and therefore has AMR capabilities. However, subcycling in time is not supported. Because chombo-streamer is a large source code, most of the code documentation resides in the source code itself, which cannot be discussed in full. Instead, we provide an outline to the structure of chombo-streamer below. In broad strokes, we require that the user

* Builds an executable
* Provides options to this executable when it is run. 

We do have a @ref worked-example, but we strongly recommend that you finish reading this page before working your way through that example. 

The typical use of chombo-streamer consists of building and running mini-applications (i.e. executables) which obtain simulation options through an input script or the command line. These applications consist of a main file that instantiates various C++ implementations of important base classes, and the must therefore be built by the user. Because the main-file is almost identical across mini-apps, an automated setup of mini-apps will likely be supported in the future. In each case, the mini-app consists of six implementation classes:

* plasma_kinetics An abstract class that defines the overall kinetics; this includes description of source terms, velocities, diffusion coefficients, surface kinetics (e.g. secondary emission) and so on. Currently, we only have kinetic models for \f$\textrm{N}_2-\textrm{O}_2\f$ mixtures. 
* computational_geometry An implementation of the geometry that will be simulated. Various implementation of this class exist, which you may immediately use. Descriptions of new geometries must be done by the user by either implementing a new computational_geometry class, or use the scripted geometry class. 
* physical_domain The physical domain to be simulated. This is a very lightweight class that only describes the axis-aligned box that you wish to simulate. There is no need to reimplement this class. 
* time_stepper The temporal integrator; the class includes an implementation of the time stepping scheme. Implementation of new integration schemes is time consuming, and additional implementations of this class is a developer task. Currently, we support some implicit-explicit schemes and Runge-Kutta schemes. Most users will find the second order Runge-Kutta scheme to be sufficient. This class owns all the individual solvers, and has access to amr_mesh.
* cell_tagger Class that is responsible for refinement and coarsening decision. 
* amr_mesh The AMR mesh engine. This class includes grid generators, coarsening operators, ghost cell interpolation operators and so on. This class the most important base class, and is responsible for orchestratic all spatial operations.

The above classes are used to instantiate another class

* plasma_engine Class that runs the entire simulation by advancing the solvers through time_stepper. plasma_engine is also responsible for obtaining refinement and coarsening flags from cell-tagger when amr_mesh calls for a regrid operation. Additional responsibilities of this class include writing plot and checkpoint files. 

For a more thorough description of these classes, see @ref base-classes. Most users will find it sufficient to modify only two of the above classes: plasma_kinetics and computational_geometry. Because some of the above classes will be left unmodified by the user (e.g. physical_domain and amr_mesh), they might be taken out of the user interface in later design. Setting up a mini-app from existing pieces of code (i.e. pre-defined geometries and plasma kinetics) is very simple. The code snippet below is a full mini-app for simulation
of a discharge from a needle electrode with a dielectric sphere placed in the middle of the gap. 

      #include "plasma_engine.H"   // Load plasma_engine class
      #include "rk2.H"             // Load temporal integrator instance (derived from time_stepper)
      #include "field_tagger.H"    // Load the cell tagger (derived from cell_tagger)
      #include "morrow_lowke.H"    // Load the plasma kinetics (derived from plasma_kinetics)
      #include "rod_sphere.H"      // Load the geometry (derived from computational_geometry)

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
        RefCountedPtr<plasma_engine> engine            = RefCountedPtr<plasma_engine> (new plasma_engine(physdom, compgeom, plaskin, timestepper, amr, tagger));

      	// Give the potential curve the plasma_engine and run simulation
        engine->set_potential(potential_curve);
        engine->setup_and_run();
      
      
      #ifdef CH_MPI
        CH_TIMER_REPORT();
        MPI_Finalize();
      #endif
      }

In the above example, all tunable parameters for the simulation are hidden from the view of the mini-app (with the exception of the potential curve). To avoid cluttering of mini-app main files, all tunable parameters are implemented as class options which are loaded on class construction. The options are usually provided by means of an input script (although the command line is also supported) which is loaded by Chombo's ParmParse class (which is instantiated in the example above). Upon construction, this class reads an input (and possible also command-line parameters) and build an internal table of available options. This table is then read in the default constructors of each base class. Each option contains the name of the class, and then the name of the variable. For example:

      amr.coarsest_domain = 128 128 128

defines a coarsest domain of \f$(128)^3\f$ cells. Furthermore, amr_mesh requires certain information about the maximum number of AMR levels that we will use, blocking factors, refinement ratios and so on. These are all passed through an input script, which may e.g. contain

      amr.max_amr_depth   = 4           # Maximum amr depth
      amr.max_sim_depth   = 2           # Maximum simulation depth
      amr.blocking_factor = 8           # Blocking factor
      amr.max_box_size    = 16          # Maximum allowed box size
      amr.ref_rat         = 2 2 2 2 2   # Refinement ratios

The first line above restricts that maximum number of AMR levels to 4 (i.e. 5 levels in total), and the second line places another restriction: We will only use 2 of these levels during the simulation. The rationale for this design is that certain temporal parts of a simulation may safely simulated using a coarse resolution, for example when a streamer crosses a pure gas gap. Other parts may require a finer spatial resolution, for example when the streamer strikes a dielectric surface. Thus, argument amr.max_sim_depth allows you to change the number of AMR levels when you checkpoint/restart one of your simulations. The third line contains the blocking factor; this is essentially the smallest box size that will be allowed by the grid generator. Likewise, the fourth line amr.max_box_size provides information on the largest box size. Finally, amr.ref_rat defines the refinement ratios between levels; you may use refinement ratios of 2 or 4 (and they may even be mixed). Of course, a simulation case may consist of several hundres of options. Currently, one needs to obtain these either by recycling input scripts, or directly in the class header file. In the future, we will most likely implement a mini-app automated builder that obtain all class options for you. 