Worked example {#worked-example}
================================

chombo-streamer is built from a mini-app principle where the user builds his own application from existing modules, or implements his own modules. In this section we show how to build mini-apps using a simple python interface that autogenerates the missing implementation pieces.

We will attempt to simulate a positive streamer discharge from a needle in a 2cm x 2cm domain; we do not consider the presence of dielectrics. As mentioned in @ref program-structure, applications primarily differ in the choice of plasma kinetic schemes and the geometries that are simulated. In this particular case we will use the Morrow-Lowke model for air discharges and a geometry that defines a needle slab geometry (although in this particular case we will eventually turn off the dielectric slab).

To create the mini-app, we do the following:

      python build.py -dim=2 -chombo_home=$(CHOMBO_HOME) -streamer_home=$(STREAMER_HOME) -base_dir=test_apps -app_name=rod2d -plasma_kinetics=morrow_lowke -geometry=rod_slab -time_stepper=rk2 -cell_tagger=field_tagger

The python script will then autogenerate the required C++ and makefiles that you need for compiling the application. Note that, in the above, $(CHOMBO_HOME) should be the full path to the Chombo library, and $(STREAMER_HOME) should be the full path to the chombo-streamer code. If, for example, the Chombo library resides in your home folder /home/my_username/mf-chombo, you must specify

      -chombo_home=/home/my_username/mf-chombo

However, your may also use environmental variables like $(CHOMBO_HOME), if you have that specified. The remaining arguments specify where the mini-app will be stored. In this particular case, it will be placed under ./test_apps/rod2d where you will find your main.cpp file, your makefile, and a templated inputs file which contains an aggregation of all the input options that you can use to fine-tune your application (e.g. the number of AMR levels that you will use). Note that the templated file is almost never suitable for production runs; extensive modifications to this file is usually required in order to get an application up and running.

First, we will build the executable for 2D execution. Navigate to $(STREAMER_HOME)/test_apps/rod2d and type

      make -s -j 10 main

which will compile your executable. It will be named something like main.2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex, depending on your Chombo configuration. If you have more cores available for computation, you may increase the process count (e.g. -j 10 uses 10 processes, or cores, for computation).

To run the example application, the typical usage is

      mpirun -np 10 main.2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex simulation.inputs

This will run your mini-app using 10 cores, using whatever is available in 'simulation.inputs' as input parameters to your program. Our next step is to modify our template inputs file to something more coherent. For this application, you can just copy the contents below to an input script of your choice (for example simulation.inputs). Most 


      # ====================================================================================================
      # POTENTIAL CURVE
      # ====================================================================================================
      my_first_application.potential = 1.E4
      
      # ====================================================================================================
      # PHYSICAL_DOMAIN CLASS OPTIONS
      # ====================================================================================================
      physical_domain.lo_corner = -1E-2 -1E-2    # Low corner of problem domain
      physical_domain.hi_corner =  1E-2  1E-2    # High corner of problem domain
      
      # ====================================================================================================
      # ROD_SLAB CLASS OPTIONS
      # ====================================================================================================
      rod_slab.eps0                      = 1                                # Background permittivity
      rod_slab.corner_curvatures         = 250E-6                           # Corner curvatures on slab
      rod_slab.turn_off_electrode        = false                            # Turn on/off electrode
      rod_slab.turn_off_dielectric       = true                             # Turn on/off dielectric
      rod_slab.electrode_live            = true                             # Live electrode or not
      rod_slab.electrode_radius          = 1.E-3                            # Electrode inner radius
      rod_slab.electrode_center1         = 0.0 0.0                          # Center 1
      rod_slab.electrode_center2         = 0.0 2.0                          # Center 2
      rod_slab.dielectric_permittivity   = 4.0                              # Dielectric permittivity
      rod_slab.dielectric_corner_lo      = -1.0123E-2 -2.0123E-2            # Low corner
      rod_slab.dielectric_corner_hi      =  1.0        2.0123E-2            # High corner
      
      
      # ====================================================================================================
      # AMR_MESH OPTIONS
      # ====================================================================================================
      amr.verbosity       = -1          # Controls verbosity. 
      amr.coarsest_domain = 128 128 128 # Number of cells on coarsest domain
      amr.max_amr_depth   = 3           # Maximum amr depth
      amr.max_sim_depth   = 3           # Maximum simulation depth
      amr.fill_ratio      = 1.0         # Fill ratio for grid generation
      amr.irreg_growth    = 2           # How much to grow irregular tagged cells
      amr.buffer_size     = 1           # Number of cells between grid levels
      amr.blocking_factor = 8           # Default blocking factor (16 in 3D)
      amr.max_box_size    = 16          # Maximum allowed box size
      amr.max_ebis_box    = 32          # Maximum allowed box size
      amr.ref_rat         = 4 4 2 2     # Refinement ratios
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
      plasma_engine.max_steps                       = 100           # Maximum number of steps
      plasma_engine.geometry_only                   = false         # Special option that ONLY plots the geometry
      plasma_engine.ebis_memory_load_balance        = false         # Use memory as loads for EBIS generation
      plasma_engine.write_ebis                      = false         # Write geometry to an HDF5 file (buggy, leave as false)
      plasma_engine.read_ebis                       = false         # Read EBIS when restarting a simulation
      plasma_engine.output_mode                     = full          # Output mode
      plasma_engine.output_directory                = ./            # Output directory
      plasma_engine.output_names                    = simulation    # Simulation output names
      plasma_engine.restart                         = false         # Do restart or not
      plasma_engine.restart_file                    = default_file  # Restart file
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
      time_stepper.fast_rte         = 1     # Solve RTE every this time steps
      time_stepper.fast_poisson     = 1     # Solve Poisson every this time steps
      time_stepper.min_dt           = 0.    # Minimum permitted time step
      time_stepper.max_dt           = 1.E99 # Maximum permitted time step
      time_stepper.cfl              = 0.5   # CFL number
      time_stepper.relax_time       = 10.0  # Relaxation time constant
      time_stepper.source_growth    = 100.0 # Relaxation time constant
      
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
      poisson_multifluid.gmg_verbosity     = -1        # GMG verbosity
      poisson_multifluid.gmg_pre_smooth    = 8         # Number of relaxations in downsweep
      poisson_multifluid.gmg_post_smooth   = 8         # Number of relaxations in upsweep
      poisson_multifluid.gmg_bott_smooth   = 8         # NUmber of relaxations before dropping to bottom solver
      poisson_multifluid.gmg_min_iter      = 5         # Minimum number of iterations
      poisson_multifluid.gmg_max_iter      = 32        # Maximum number of iterations
      poisson_multifluid.gmg_tolerance     = 1.E-6     # Residue tolerance
      poisson_multifluid.gmg_hang          = 0.2       # Solver hang
      poisson_multifluid.gmg_bottom_drop   = 8         # Bottom drop
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
      eddington_sp1.gmg_verbosity     = -1        # GMG verbosity
      eddington_sp1.gmg_pre_smooth    = 12        # Number of relaxations in downsweep
      eddington_sp1.gmg_post_smooth   = 12        # Number of relaxations in upsweep
      eddington_sp1.gmg_bott_smooth   = 12        # NUmber of relaxations before dropping to bottom solver
      eddington_sp1.gmg_min_iter      = 5         # Minimum number of iterations
      eddington_sp1.gmg_max_iter      = 32        # Maximum number of iterations
      eddington_sp1.gmg_tolerance     = 1.E-6     # Residue tolerance
      eddington_sp1.gmg_hang          = 0.2       # Solver hang
      eddington_sp1.gmg_bottom_drop   = 8         # Bottom drop
      eddington_sp1.gmg_bottom_solver = simple    # Bottom solver type. Valid options are 'simple' and 'bicgstab'
      eddington_sp1.gmg_bottom_relax  = 32        # Number of relaxations in bottom solve ('simple' solver only)
      
      # ====================================================================================================
      # RK2 CLASS OPTIONS
      # ====================================================================================================
      rk2.alpha = 0.5 # Set alpha. 0.5 = Heuns method, 1.0 = Midpoint method		

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

