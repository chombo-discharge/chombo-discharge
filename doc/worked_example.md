Worked example {#worked-example}
================================

chombo-streamer is built from a mini-app principle where the user builds his own application from existing modules, or implements his own modules. In this section we show how to build mini-apps using a simple python interface that autogenerates the missing implementation pieces.

We will attempt to simulate a positive streamer discharge from a needle in a 2cm x 2cm domain; we do not consider the presence of dielectrics. As mentioned in @ref program-structure, applications primarily differ in the choice of plasma kinetic schemes and the geometries that are simulated. In this particular case we will use the Morrow-Lowke model for air discharges and a geometry that defines a needle slab geometry (although in this particular case we will eventually turn off the dielectric slab).

To create the mini-app, we do the following:

      python build.py -dim=2 -chombo_home=$(CHOMBO_HOME) -streamer_home=$(STREAMER_HOME) -base_dir=test_apps -app_name=rod_sphere2d -plasma_kinetics=morrow_lowke -geometry=rod_sphere -time_stepper=rk2 -cell_tagger=field_tagger

The python script will then autogenerate the required C++ and makefiles that you need for compiling the application. Note that, in the above, $(CHOMBO_HOME) should be the full path to the Chombo library, and $(STREAMER_HOME) should be the full path to the chombo-streamer code. If, for example, the Chombo library resides in your home folder /home/my_username/mf-chombo, you must specify

      -chombo_home=/home/my_username/mf-chombo

However, your may also use environmental variables like $(CHOMBO_HOME), if you have that specified. The remaining arguments specify where the mini-app will be stored. In this particular case, it will be placed under ./test_apps/rod2d where you will find your main.cpp file, your makefile, and a templated inputs file which contains an aggregation of all the input options that you can use to fine-tune your application (e.g. the number of AMR levels that you will use). Note that the templated file is almost never suitable for production runs; extensive modifications to this file is usually required in order to get an application up and running.

First, we will build the executable for 2D execution. Navigate to $(STREAMER_HOME)/test_apps/rod2d and type

      make -s -j 10 main

which will compile your executable. It will be named something like main.2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex, depending on your Chombo configuration. If you have more cores available for computation, you may increase the process count (e.g. -j 10 uses 10 processes, or cores, for computation).

To run the example application, the typical usage is

      mpirun -np 10 main.2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex simulation.inputs

This will run your mini-app using 10 cores, using whatever is available in 'simulation.inputs' as input parameters to your program. Our next step is to modify our template inputs file to something more coherent. For this application, you can just copy the contents below to an input script of your choice (for example simulation.inputs). 

      # ====================================================================================================
      # POTENTIAL CURVE
      # ====================================================================================================
      rod_sphere2d.potential = 1.2E4
      
      # ====================================================================================================
      # PHYSICAL_DOMAIN CLASS OPTIONS
      # ====================================================================================================
      physical_domain.lo_corner = -2E-3 -2E-3    # Low corner of problem domain
      physical_domain.hi_corner =  2E-3  2E-3    # High corner of problem domain
      
      # ====================================================================================================
      # ROD_SLAB CLASS OPTIONS
      # ====================================================================================================
      rod_sphere.eps0                      = 1.0                              # Background permittivity
      rod_sphere.turn_off_electrode        = false                            # Turn on/off electrode
      rod_sphere.turn_off_dielectric       = true                             # Turn on/off dielectric
      rod_sphere.electrode_live            = true                             # Live electrode or not
      rod_sphere.electrode_radius          = 200.E-6                          # Electrode radius
      rod_sphere.electrode_center1         = 0.0 0.0 0E-2                     # Center 1
      rod_sphere.electrode_center2         = 0.0 1.0 0E-2                     # Center 2
      rod_sphere.dielectric_permittivity   = 5.0                              # Dielectric permittivity
      rod_sphere.dielectric_center         = 0.0 -500.E-6 0.0                 # Dielectric center
      rod_sphere.dielectric_radius         = 200.E-6                         # Dielectric radius

      # ====================================================================================================
      # AMR_MESH OPTIONS
      # ====================================================================================================
      amr.verbosity       = -1          # Controls verbosity. 
      amr.coarsest_domain = 128 128     # Number of cells on coarsest domain
      amr.max_amr_depth   = 4           # Maximum amr depth
      amr.max_sim_depth   = 4           # Maximum simulation depth
      amr.fill_ratio      = 1.0         # Fill ratio for grid generation
      amr.irreg_growth    = 2           # How much to grow irregular tagged cells
      amr.buffer_size     = 1           # Number of cells between grid levels
      amr.blocking_factor = 8           # Default blocking factor (16 in 3D)
      amr.max_box_size    = 16          # Maximum allowed box size
      amr.max_ebis_box    = 64          # Maximum allowed box size
      amr.ref_rat         = 2 2 2 2 2   # Refinement ratios
      amr.num_ghost       = 3           # Number of ghost cells. Default is 3
      amr.eb_ghost        = 4           # Set number of of ghost cells for EB stuff
      amr.redist_radius   = 1           # Redistribution radius for hyperbolic conservation laws
      amr.stencil_order   = 1           # Order for interpolation/extrapolation stencils. Must be 1 or 2
      amr.stencil_radius  = 1           # Radius for interpolation/extrapolation stencils. Must be 1 or 2
      amr.stencil_type    = linear      # Default stencil type. Valid options are 'linear', 'taylor', 'lsq'
      amr.ghost_interp    = pwl         # Ghost cell interpolation type. Valid options are 'pwl' or 'quad'
      amr.load_balance    = knapsack    # Load balancing algorithm. Valid options are 'knapsack' or 'elliptic'
      amr.ebcf            = false       # Tell amr to ignore EBCF-related code. 

      # ====================================================================================================
      # PLASMA_ENGINE OPTIONS
      # ====================================================================================================
      plasma_engine.verbosity                       = 2             # Engine verbosity
      plasma_engine.plot_interval                   = 10            # Plot interval
      plasma_engine.regrid_interval                 = 10            # Regrid interval
      plasma_engine.checkpoint_interval             = 10            # Checkpoint interval
      plasma_engine.init_regrids                    = 0             # Number of initial regrids
      plasma_engine.start_time                      = 0             # Start time (fresh simulations only)
      plasma_engine.stop_time                       = 1.0           # Stop time
      plasma_engine.max_steps                       = 1000          # Maximum number of steps
      plasma_engine.geometry_only                   = false         # Special option that ONLY plots the geometry
      plasma_engine.ebis_memory_load_balance        = false         # Use memory as loads for EBIS generation
      plasma_engine.write_ebis                      = false         # Write geometry to an HDF5 file (buggy, leave as false)
      plasma_engine.read_ebis                       = false         # Read EBIS when restarting a simulation
      plasma_engine.output_mode                     = full          # Output mode
      plasma_engine.output_directory                = ./            # Output directory
      plasma_engine.output_names                    = simulation    # Simulation output names
      plasma_engine.restart                         = false         # Do restart or not
      plasma_engine.restart_step                    = 0             # Restart step
      plasma_engine.restart_mode                    = full          # Restart mode. 
      plasma_engine.refine_geometry                 = -1            # Refine geometry, -1 => Refine all the way down
      plasma_engine.refine_electrodes               = -1            # Refine electrode surfaces. -1 => equal to refine_geometry
      plasma_engine.refine_dielectrics              = -1            # Refine dielectric surfaces. -1 => equal to refine_geometry
      plasma_engine.refine_electrode_gas_interface  = -1            # Refine electrode-gas interfaces. -1 => ----"-----
      plasma_engine.refine_dielectric_gas_interface = -1            # Refine dielectric-gas interfaces. -1 => ----"-----
      plasma_engine.refine_solid_gas_interface      = -1            # Refine solid-gas interfaces. -1 => ----"-----
      plasma_engine.refine_solid_solid_interface    = -1            # Refine solid-solid interfaces. -1 => ----"-----
      
      # ====================================================================================================
      # TIME_STEPPER OPTIONS
      # ====================================================================================================
      time_stepper.verbosity        = -1    # Class verbosity
      time_stepper.solver_verbosity = -1    # Individual solver verbosities
      time_stepper.fast_rte         = 5     # Solve RTE every this time steps
      time_stepper.fast_poisson     = 1     # Solve Poisson every this time steps
      time_stepper.min_dt           = 0.    # Minimum permitted time step
      time_stepper.max_dt           = 1.E99 # Maximum permitted time step
      time_stepper.cfl              = 0.1   # CFL number
      time_stepper.relax_time       = 10.0  # Relaxation time constant
      time_stepper.source_growth    = 1000.0 # Relaxation time constant
      
      # ====================================================================================================
      # POISSON_SOLVER CLASS OPTIONS
      # ====================================================================================================
      poisson_solver.bc_x_low  = neumann               # BC type. Valid options are "neumann", "dirichlet_ground", and "dirichlet_live"
      poisson_solver.bc_x_high = neumann               # BC type. Valid options are "neumann", "dirichlet_ground", and "dirichlet_live"
      poisson_solver.bc_y_low  = dirichlet_ground      # BC type. Valid options are "neumann", "dirichlet_ground", and "dirichlet_live"
      poisson_solver.bc_y_high = dirichlet_live        # BC type. Valid options are "neumann", "dirichlet_ground", and "dirichlet_live"

      # ====================================================================================================
      # POISSON_MULTIFLUID_GMG CLASS OPTIONS (MULTIFLUID GMG SOLVER SETTINGS)
      # ====================================================================================================
      poisson_multifluid.gmg_verbosity     = -10       # GMG verbosity
      poisson_multifluid.gmg_pre_smooth    = 16         # Number of relaxations in downsweep
      poisson_multifluid.gmg_post_smooth   = 16         # Number of relaxations in upsweep
      poisson_multifluid.gmg_bott_smooth   = 16         # NUmber of relaxations before dropping to bottom solver
      poisson_multifluid.gmg_min_iter      = 5         # Minimum number of iterations
      poisson_multifluid.gmg_max_iter      = 32        # Maximum number of iterations
      poisson_multifluid.gmg_tolerance     = 1.E-6     # Residue tolerance
      poisson_multifluid.gmg_hang          = 0.2       # Solver hang
      poisson_multifluid.gmg_bottom_drop   = 4         # Bottom drop
      poisson_multifluid.gmg_bottom_solver = bicgstab  # Bottom solver type. Valid options are 'simple' and 'bicgstab'
      poisson_multifluid.gmg_bottom_relax  = 32        # Number of relaxations in bottom solve ('simple' solver only)

      # ====================================================================================================
      # CDR_LAYOUT SOLVER SETTINGS
      # ----------------------------------------------------------------------------------------------------
      cdr_layout.which_solver = godunov                # Solver type, available options are "scharfetter-gummel" and "godunov"
      cdr_gdnv.divF_nc        = covered_face           # Valid options are "covered_face" and "conservative_average"
      cdr_gdnv.limit_slopes   = true                   # Valid options are "covered_face" and "conservative_average"

      # ====================================================================================================
      # RTE_LAYOUT CLASS OPTIONS
      # ====================================================================================================
      rte_layout.which_solver = eddington_sp1 # Solver type. Available option is "eddington_sp1"
      rte_layout.stationary   = true          # Stationary solver ot no.

      # ====================================================================================================
      # EDDINGTON_SP1 CLASS OPTIONS
      # ====================================================================================================
      eddington_sp1.reflectivity      = 0.        # Reflectivity
      eddington_sp1.gmg_verbosity     = -10       # GMG verbosity
      eddington_sp1.gmg_pre_smooth    = 6         # Number of relaxations in downsweep
      eddington_sp1.gmg_post_smooth   = 6         # Number of relaxations in upsweep
      eddington_sp1.gmg_bott_smooth   = 6         # NUmber of relaxations before dropping to bottom solver
      eddington_sp1.gmg_min_iter      = 5         # Minimum number of iterations
      eddington_sp1.gmg_max_iter      = 32        # Maximum number of iterations
      eddington_sp1.gmg_tolerance     = 1.E-6     # Residue tolerance
      eddington_sp1.gmg_hang          = 0.2       # Solver hang
      eddington_sp1.gmg_bottom_drop   = 4         # Bottom drop
      eddington_sp1.gmg_bottom_solver = bicgstab  # Bottom solver type. Valid options are 'simple' and 'bicgstab'
      eddington_sp1.gmg_bottom_relax  = 32        # Number of relaxations in bottom solve ('simple' solver only)
      
      # ====================================================================================================
      # RK2 CLASS OPTIONS
      # ====================================================================================================
      rk2.alpha = 0.5 # Set alpha. 0.5 = Heuns method, 1.0 = Midpoint method

      # ====================================================================================================
      # FIELD_TAGGER CLASS OPTIONS
      # ====================================================================================================
      field_tagger.coarsen_curvature = 0.1    # Sets eps_curv for coarsening
      field_tagger.coarsen_magnitude = 0.1    # Sets eps_mag for coarsening 
      field_tagger.refine_curvature  = 0.2    # Sets eps_curv for refinement
      field_tagger.refine_magnitude  = 0.75   # Sets eps_mag for refinement

      # ====================================================================================================
      # MORROW_LOWKE CLASS OPTIONS
      # ====================================================================================================
      morrow_lowke.gas_temperature               = 300                   # Gas temperature (K)
      morrow_lowke.gas_N2_frac                   = 0.8                   # Mixing fraction of nitrogen
      morrow_lowke.gas_O2_frac                   = 0.2                   # Mixing fraction of oxygen
      morrow_lowke.gas_pressure                  = 1.0                   # Gas pressure (atm)
      morrow_lowke.gas_quenching_pressure        = 0.03947               # Quenching pressure for photo-emission (atm)
      morrow_lowke.positive_species_mobility     = 2.E-4                 # Positive species mobility
      morrow_lowke.negative_species_mobility     = 2.E-4                 # Negative species mobility
      morrow_lowke.excitation_efficiency         = 0.6                   # Impact excitation efficiency
      morrow_lowke.photoionization_efficiency    = 0.1                   # Photo-ionization efficiency
      morrow_lowke.electrode_townsend2           = 1.E-6                 # Townsend coefficient on electrodes
      morrow_lowke.electrode_quantum_efficiency  = 1.E-6                 # Quantum efficiency on electrodes
      morrow_lowke.dielectric_townsend2          = 1.E-6                 # Townsend coefficient on dielectrics
      morrow_lowke.dielectric_quantum_efficiency = 1.E-6                 # Quantum efficiency on dielectrics
      morrow_lowke.dielectric_work_function      = 1.E-6                 # Dielectric work function (ev)
      morrow_lowke.uniform_density               = 1.E10                 # Initial background ionization
      morrow_lowke.seed_density                  = 0.E16                 # Initial seed ionization
      morrow_lowke.seed_radius                   = 1E-3                  # Initial seed radius
      morrow_lowke.seed_position                 = 0. 0. 0.              # Initial seed position
      morrow_lowke.noise_amplitude               = 0.0                   # Initial noise amplitude
      morrow_lowke.noise_octaves                 = 1                     # Initial noise octaves
      morrow_lowke.noise_persistence             = 0.5                   # Reduction factor between each noise octave
      morrow_lowke.noise_frequency               = 1.0 1.0 1.0           # Spatial noise frequency
      morrow_lowke.photon1_A_coeff               = 1.12E-4               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon1_lambda_coeff          = 4.15E-2               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon2_A_coeff               = 2.88E-2               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon2_lambda_coeff          = 1.09E-1               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon3_A_coeff               = 2.76E-1               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon3_lambda_coeff          = 6.69E-1               # Parameters from Bourdon et. al photoionization model


The first step to running new mini-apps is usually to inspect the geometry and initial mesh. There is a special option in plasma_engine that allows us to dump the geometry to an HDF5 file. Our first command is therefore

      mpirun -np 10 main2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex simulation.inputs plasma_engine.geometry_only=true

where the option plasma_engine.geometry_only overrides whatever this option is set to in the input script. When you run the application, it will write a file called simulation.geometry.2d.hdf5. Note that you can call your simulation whatever you want, and have the data dump to a folder of your choice. This is controlled through plasma_engine.output_names and plasma_engine.output_directory. Modify them as you see fit. In addition, you will get a number of pout.# files. These files are the output of every process through which you can monitor mesh size, resolution, time step size, elapsed time and so on. Let us investigate the geometry file. Open the file in VisIt and plot the Filled Boundary and Mesh. You should see something like this (exact plots vary with VisIt versions):


<img src="./simulation_geometry.png" alt="Initial geometry and mesh" style="width: 1024px;"/>

That is our initial mesh and geometry, which looks what we want. Next, we will perform some screening simulations of this geometry. The current input script uses a resolution of about 2 microns and a CFL number of 0.1. Before we look for grid convergence, we attempt to resolve the simulation using larger spatial and temporal steps, at least for a short time. 


      mpirun -np 10 main2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex simulation.inputs amr.max_sim_depth=3 time_stepper.cfl=0.8 time_stepper.max_dt = 0.5E-9

The option amr.max_sim_depth=3 restricts to a maximum number of 3 AMR levels; the other two species a CFL number of 0.8 and a maximum simulation time of 0.5ns. To make this simulation run even faster, it might even be justifiable to only solve the RTE equations at certain intervals. To try this, you may add an option time_stepper.fast_rte=5, which will solve the RTE equations every 5 time steps. Depending on your hardware, the execution time for each time can be anywhere from a few hundred milliseconds to many seconds.

When the program has finished running, it will have executed 86 time steps. You will find plotfiles and checkpoint files in your specified directory. In our input script, we've specified that these should be written every 10 timesteps. This is controlled through plasma_engine.plot_interval and plasma_engine.checkpoint_interval. The first and final timesteps are always written, as long as plasma_engine.plot_interval > 0. Note that there is also an option plasma_engine.regrid_interval which specifies how often one should regrid. After 86 time steps, we examine the plasma (electron density):

<img src="./plasma_density_500ps.png" alt="plasma density" style="width: 1024px;"/>

Here, we've plotted the data on a logarithmic scale. Realizing that we might need better resolution when the streamer starts propagating, we will simulate the next 0.5ns by using a finer mesh with 2 micron resolution and a CFL number of 0.4. We will restart from where we left off, time step 86:

      mpirun -np 10 main2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex simulation.inputs plasma_engine.restart=true plasma_engine.restart_step=86

The options above will restart our simulation from time step 86. Note that you may change a bunch of other options as well, including the way the grids are handled. However, you may NOT change the AMR levels since this interferes with interpretatin of where the data lives in space. On other hand, you can restrict amr.max_sim_depth as you see fit.
