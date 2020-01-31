#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>
using std::ifstream;
using std::ios;

#include "FABView.H"
// this lets us use dumpIVS, other dump functions
#include "DebugDump.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"

#include "AMR.H"
#include "AMRLevel.H"
#include "AMRLevelMappedConsFactory.H"
#include "AMRLevelMappedCons.H"
#include "LevelGridMetrics.H"
#include "EulerEquationsMappedStabilityStrategy.H"
#include "DensityGradientMappedTaggingStrategy.H"
#include "GaussianAnalyticMappedTaggingStrategy.H"

#include "MOLPolytropicPhysics.H"

//#include "RampIBC.H"
//#include "ExplosionIBC.H"
//#include "Explosion1dIBC.H"
//#include "GaussianIBC.H"
//#include "Gaussian1dIBC.H"
//#include "Gaussian1dvIBC.H"
//#include "GaussianPBC.H"
//#include "GaussianSmoothBC.H"
//#include "CosinePowerIBC.H"
//#include "CosineEachPowerIBC.H"
//#include "VortexIBC.H"
//#include "ShearIBC.H"
//#include "TrigIBC.H"
//#include "SineIBC.H"
//#include "WaveIBC.H"
//#include "SourceIBC.H"
//#include "ShockTubeIBC.H"
//#include "ChannelShockIBC.H"
//#include "ChannelModianoIBC.H"

//--Switching to mapped IBCs

#include "RampMappedIBC.H"
#include "GaussianMappedIBC.H"
#include "TrigMappedIBC.H"

#include "CartesianCS.H"
#include "RThetaZCS.H"
#include "RThetaPhiCS.H"
#include "TwistedCS.H"
#include "WarpedCS.H"
#include "SchwarzChristoffelRampCS.H"

#include "generalFuncs.H"

// #define HALEM_PROC_SPEED
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#ifdef CH_Linux
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

// Possible pressure relationships for the initial condition
#define PRESSURE_ISENTROPIC 0
#define PRESSURE_CONSTANT   1

// amrRun is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrRun();

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

// setupVortices allows vortex parameters to be read in and used in this AMR
// computation example
void setupVortices(Vector<RealVect>& a_center,
                   Vector<Real>&     a_radius,
                   Vector<Real>&     a_strength,
                   int               a_verbosity,
                   std::string       a_vortexFile);

// One more function for MPI
void dumpmemoryatexit();

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
#ifdef CH_AIX
  H5dont_atexit();
#endif
  setChomboMPIErrorHandler();
#endif

  int rank, number_procs;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
  rank = 0;
  number_procs = 1;
#endif

  if (rank == 0)
  {
    pout() << " number_procs = " << number_procs << endl;
  }

  // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrRun...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

#ifdef TRAP_FPE
  enableFpExceptions ();
#endif

  // Run amrRun, i.e., do the computation
  amrRun();

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

void amrRun()
{
  // Read inputs that are prefixed with "godunov."
  ParmParse ppgodunov("godunov");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppgodunov.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // For all gas dynamics
  Real gamma = 1.4;
  ppgodunov.get("gamma",gamma);

  // Stop after this number of steps
  int nstop = 0;
  ppgodunov.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppgodunov.get("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppgodunov.query("domain_length",domainLength);

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i)
  {
    numCells[i] = 0;
  }
  ppgodunov.getarr("num_cells",numCells,0,SpaceDim);
  // Check that every component of numCells is positive and even
  for (int i = 0; i < SpaceDim; ++i)
    {
      CH_assert(numCells[i] > 0);
      CH_assert(numCells[i] % 2 == 0);
    }

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppgodunov.getarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim = 0; dim < SpaceDim; dim++)
  {
    isPeriodic[dim] = (isPeriodica[dim] == 1);
    if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
    {
      pout() << "Using Periodic BCs in direction: " << dim << endl;
    }
  }

  // Maximum AMR level limit
  int maxLevel = 0;
  ppgodunov.get("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppgodunov.getarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppgodunov.getarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppgodunov.get("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppgodunov.get ("refine_thresh",refineThresh);

  // Whether refinement threshold is scaled with dx
  int refinementIsScaledInt = 0;
  ppgodunov.query("refinement_is_scaled", refinementIsScaledInt);
  bool refinementIsScaled = (refinementIsScaledInt == 1);

  // Whether to tag on pressure instead of on density
  int tagPressureInt = 0;
  ppgodunov.query("tag_pressure", tagPressureInt);
  bool tagPressure = (tagPressureInt == 1);

  // Whether to tag on vorticity instead of on density
  int tagVorticityInt = 0;
  ppgodunov.query("tag_vorticity", tagVorticityInt);
  bool tagVorticity = (tagVorticityInt == 1);

  // Name of tagging strategy
  std::string tagStrategyName("DensityGradient");
  ppgodunov.query("tag_strategy", tagStrategyName);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppgodunov.get("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppgodunov.get("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppgodunov.get("fill_ratio",fillRatio);

  // Grid buffer size
  int gridBufferSize = 1;
  const int codeReqBufferSize = LevelGridMetrics::bufferSize4thO(
    refRatios,
    maxLevel,
    5);  // Num ghost 5 is hard coded in AMRLevelMappedCons
  const int haveUserBufferSize =
    ppgodunov.query("grid_buffer_size", gridBufferSize);
  if (haveUserBufferSize)
    {
      if (gridBufferSize < codeReqBufferSize)
        {
          pout() << "\nWARNING: Program requested grid buffer size: "
                 << codeReqBufferSize
                 << "\n         User requested grid buffer size   : "
                 << gridBufferSize
                 << "\nUsing a buffer size that is too small may corrupt "
            "the solution.\n";
          MayDay::Warning("Do not specify the grid buffer size to avoid this "
                          "warning");
          pout() << std::endl;
        }
    }
  else
    {
      gridBufferSize = codeReqBufferSize;
    }

  // Order of the normal predictor (CTU -> 0, PLM -> 1, PPM -> 2)
  std::string normalPred;
  int normalPredOrder = 0;
  ppgodunov.get("normal_predictor",normalPred);
  if (normalPred == "CTU" || normalPred == "ctu")
  {
    normalPredOrder = 0;
  }
  else if (normalPred == "PLM" || normalPred == "plm")
  {
    normalPredOrder = 1;
  }
  else if (normalPred == "PPM" || normalPred == "ppm")
  {
    normalPredOrder = 2;
  }
  else
  {
    MayDay::Error("Normal precdictor must by PLM or PPM");
  }

  // Use fourth order slopes:  default true
  int inFourthOrderSlopes = 1;
  bool useFourthOrderSlopes;
  ppgodunov.query("use_fourth_order_slopes",inFourthOrderSlopes);
  useFourthOrderSlopes = (inFourthOrderSlopes == 1);

  // Arbitrary fiat by petermc, 17 June 2008
  useFourthOrderSlopes = true;

  // Do slope limiting:  default true
  int inPrimLimiting = 1;
  bool usePrimLimiting;
  ppgodunov.query("use_prim_limiting",inPrimLimiting);
  usePrimLimiting = (inPrimLimiting == 1);

  // This should actually be 1 even if slope limiting is off
  int highOrderLimiterInt = 1;
  ppgodunov.query("high_order_limiter", highOrderLimiterInt);
  bool highOrderLimiter = (highOrderLimiterInt == 1);

  // NEW Kreiss-Oliger artificial viscosity
  int inArtVisc = 0;
  ppgodunov.query("use_art_visc", inArtVisc);
  bool useArtVisc = (inArtVisc == 1);

  Real ratioArtVisc = 0.;
  if (useArtVisc)
    {
      ppgodunov.query("ratio_art_visc", ratioArtVisc);
    }

  int inForwardEuler = 0;
  ppgodunov.query("forward_euler", inForwardEuler);
  bool forwardEuler = (inForwardEuler == 1);

  // Do slope limiting using characteristics
  int inCharLimiting = 0;
  bool useCharLimiting;
  ppgodunov.get("use_char_limiting",inCharLimiting);
  useCharLimiting = (inCharLimiting == 1);

  // Initial values are average:  default false
  int inInitialAverage = 0;
  bool initialAverage;
  ppgodunov.query("initial_average", inInitialAverage);
  initialAverage = (inInitialAverage == 1);

  // Do slope flattening:  default false
  int inFlattening = 0;
  bool useFlattening;
  ppgodunov.query("use_flattening", inFlattening);
  useFlattening = (inFlattening == 1);

  // Avoid PPM:  default false
  int inNoPPM = 0;
  bool noPPM;
  ppgodunov.query("no_ppm", inNoPPM);
  noPPM = (inNoPPM == 1);

  // Do deconvolution:  default true
  int inDoDeconvolution = 1;
  bool doDeconvolution;
  ppgodunov.query("do_deconvolution", inDoDeconvolution);
  doDeconvolution = (inDoDeconvolution == 1);

  // Do deconvolution:  default true
  int inDoFaceDeconvolution = 1;
  bool doFaceDeconvolution;
  ppgodunov.query("do_face_deconvolution", inDoFaceDeconvolution);
  doFaceDeconvolution = (inDoFaceDeconvolution == 1);

  // Apply artificial viscosity based on divergence.
  int inArtificialViscosity = 1;
  bool useArtificialViscosity;
  ppgodunov.get("use_artificial_viscosity",inArtificialViscosity);
  useArtificialViscosity = (inArtificialViscosity == 1);

  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.1;
  ppgodunov.get("artificial_viscosity",artificialViscosity);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppgodunov.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppgodunov.query("plot_interval",plotInterval);

  // CFL multiplier
  Real cfl = 0.8;
  ppgodunov.get("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppgodunov.get("initial_cfl",initialCFL);

  // Set up whether to use subcycling in time.
  int useSubcyclingInt = 1;
  ppgodunov.query("use_subcycling", useSubcyclingInt);
  bool useSubcycling = (useSubcyclingInt == 1);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppgodunov.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppgodunov.get("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppgodunov.get("dt_tolerance_factor",dtToleranceFactor);

  // Create and define IBC (initial and boundary condition) object
  PhysIBC* ibc = 0;

  // A minimum pressure needed to construct MOLPolytropicPhysics - used in slope
  // flattening
  Real smallPressure;

  // Don't use source term by default
  bool useSourceTerm = false;

  // Source term multiplier
  Real sourceTermScaling = 0.0;

  ProblemDomain probDomain(IntVect::Zero,
                           IntVect(D_DECL6(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1,
                                           numCells[3]-1,
                                           numCells[4]-1,
                                           numCells[5]-1)),
                           isPeriodic);

  // Coordinate system for mapped grids
  NewCoordSysFactory* coordSysFact = 0;

/*--------------------------------------------------------------------*
 * Determine the sample problem specified
 *--------------------------------------------------------------------*/

  std::string problemString;
  if (ppgodunov.contains("problem"))
  {
    ppgodunov.query("problem",problemString);

    // Print some parameters
    if (verbosity >= 2)
    {
      pout() << "problem = " << problemString << endl;
      pout() << "gamma = " << gamma << endl;
    }

    if (problemString == "ramp")
    {
      if (isPeriodic[0] || isPeriodic[1])
      {
        MayDay::Error("Neither x or y boundaries can be periodic");
      }

      // Ramp problem
      Real alpha = 30.0;
      ppgodunov.get("angle_deg",alpha);

      Real ms = 10.0;
      ppgodunov.get("shock_mach",ms);

      Real xcorner = 0.1;
      ppgodunov.get("xcorner",xcorner);

      if (verbosity >= 2)
      {
        pout() << "alpha = " << alpha << endl;
        pout() << "shock_mach = " << ms << endl;
        pout() << "xcorner = " << xcorner << endl;
      }

      // Define IBC for ramp problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       RampIBC* rampibc = new RampIBC;
//       rampibc->setFortranCommon(smallPressure,
//                                 gamma,
//                                 alpha,
//                                 ms,
//                                 xcorner,
//                                 artificialViscosity);
//       ibc = rampibc;
    }
    if (problemString == "mappedramp")
    {
      if (isPeriodic[0] || isPeriodic[1])
        {
          MayDay::Error("Neither x or y boundaries can be periodic");
        }

      // Ramp problem
      Real alpha = 30.0;
      ppgodunov.get("angle_deg",alpha);

      Real ms = 10.0;
      ppgodunov.get("shock_mach",ms);

      Real leadLength = 0.5;
      ppgodunov.get("lead_length", leadLength);

      Real rampLength = 4.0;
      ppgodunov.get("ramp_length", rampLength);

      Real X0Start = -leadLength/2;
      ppgodunov.get("shock_location", X0Start);


      if (verbosity >= 2)
        {
          pout() << "alpha = " << alpha << endl;
          pout() << "shock_mach = " << ms << endl;
          pout() << "lead_length = " << leadLength << endl;
          pout() << "ramp_length = " << rampLength << endl;
          pout() << "shock_location = " << X0Start << endl;
        }

      const Real alphaRad = alpha*Pi/180.;
      // Define the coordinate system for the ramp problem
      coordSysFact = new SchwarzChristoffelRampCSFactory(numCells[0],
                                                         alphaRad,
                                                         leadLength,
                                                         rampLength);

      // Reset the domain length in mapped space to 1.0
      domainLength = 1.0;
      if (verbosity >= 2)
        {
          pout() << "Computational domainLength fixed to 1.0 for mappedramp "
            "problem" << endl;
        }

      // Define IBC for ramp problem
      pout() << problemString << std::endl;
      RampMappedIBC* rampibc = new RampMappedIBC(smallPressure,
                                                 gamma,
                                                 alphaRad,
                                                 ms,
                                                 X0Start,
                                                 artificialViscosity);

      // The coordinate system needs to be available for any future IBC created
      // from this one.
      rampibc->setCoordSys(
        static_cast<NewFourthOrderCoordSys*>(
          coordSysFact->getCoordSys(
            probDomain, (domainLength/numCells[0])*RealVect::Unit)));
      ibc = rampibc;
    }
    else if (problemString == "channelShock")
    {
      if (isPeriodic[0])
      {
        MayDay::Error("The x boundary can't be periodic");
      }

      // Shock problem in a channel
      Real ms = 10.0;
      ppgodunov.get("shock_mach",ms);

      Real position = 0.1;
      ppgodunov.get("position",position);

      if (verbosity >= 2)
      {
        pout() << "shock_mach = " << ms << endl;
        pout() << "position = " << position << endl;
      }

      // Define IBC for channel shock problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       ChannelShockIBC* channelshockibc = new ChannelShockIBC;
//       channelshockibc->setFortranCommon(smallPressure,
//                                         gamma,
//                                         ms,
//                                         position,
//                                         artificialViscosity);
//       ibc = channelshockibc;
    }
    else if (problemString == "channelModiano")
    {
      if (isPeriodic[0])
      {
        MayDay::Error("The x boundary can't be periodic");
      }

      // Modiano problem in a channel
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      Real center = 0.5;
      ppgodunov.get("center",center);

      Real width = 0.5;
      ppgodunov.get("width",width);

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "center = " << center << endl;
        pout() << "width = " << width << endl;
      }

      // Define IBC for channel Modiano problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       ChannelModianoIBC* channelmodianoibc = new ChannelModianoIBC;
//       channelmodianoibc->setFortranCommon(smallPressure,
//                                           gamma,
//                                           ambientDensity,
//                                           deltaDensity,
//                                           center,
//                                           width,
//                                           artificialViscosity);
//       ibc = channelmodianoibc;
    }
    else if (problemString == "explosion")
    {
      // Explosion problem
      Real ms = 10.0;
      ppgodunov.get("shock_mach",ms);

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppgodunov.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "shock_mach = " << ms << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for explosion problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       ExplosionIBC* explosionibc = new ExplosionIBC;
//       explosionibc->setFortranCommon(smallPressure,
//                                      gamma,
//                                      ms,
//                                      center,
//                                      size,
//                                      velocity,
//                                      artificialViscosity);
//       ibc = explosionibc;
    }
    else if (problemString == "explosion1d")
    {
      // Explosion1d problem
      Real ms = 10.0;
      ppgodunov.get("shock_mach",ms);

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppgodunov.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      // use only component 0; set the other components to 0.
      velocity[0] = velocitypp[0];
      for (int i = 1; i < SpaceDim; i++)
      {
        velocity[i] = 0.;
      }

      if (verbosity >= 2)
      {
        pout() << "shock_mach = " << ms << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for explosion problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       Explosion1dIBC* explosion1dibc = new Explosion1dIBC;
//       explosion1dibc->setFortranCommon(smallPressure,
//                                        gamma,
//                                        ms,
//                                        center,
//                                        size,
//                                        velocity,
//                                        artificialViscosity);
//       ibc = explosion1dibc;
    }
    else if (problemString == "shockTube")
    {
      ParmParse ppshocktube("shockTube");
      Real center;
      ppshocktube.get("initial_center",center);

      int direction = 0;
      ppshocktube.query("direction", direction);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppshocktube.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
        {
          velocity[i] = velocitypp[i];
        }

      // define left and right states
      Real rhoLeft, rhoRight, eLeft, eRight;
      ppshocktube.get("rhoLeft", rhoLeft);
      ppshocktube.get("rhoRight", rhoRight);
      ppshocktube.get("eLeft", eLeft);
      ppshocktube.get("eRight", eRight);

      if (verbosity >= 2)
        {
          pout() << "rhoLeft = " << rhoLeft << endl;
          pout() << "rhoRight = " << rhoRight << endl;
          pout() << "eLeft = " << eLeft << endl;
          pout() << "eRight = " << eRight << endl;

          pout() << "initial_center = " << center << endl;
          pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                     velocity[1] << "  " <<,
                                                     velocity[2] << "  " <<,
                                                     velocity[3] << "  " <<,
                                                     velocity[4] << "  " <<,
                                                     velocity[5] << ) endl;
        }

      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       ShockTubeIBC* shocktubeibc = new ShockTubeIBC;
//       shocktubeibc->setFortranCommon(smallPressure,
//                                      gamma,
//                                      rhoLeft, rhoRight,
//                                      eLeft, eRight,
//                                      center,
//                                      direction,
//                                      velocity,
//                                      artificialViscosity);
//       ibc = shocktubeibc;
    }
    else if (problemString == "gaussian")
    {
      // Gaussian problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppgodunov.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       GaussianIBC* gaussianibc = new GaussianIBC;
//       gaussianibc->setFortranCommon(smallPressure,
//                                     gamma,
//                                     ambientDensity,
//                                     deltaDensity,
//                                     pressure,
//                                     center,
//                                     size,
//                                     velocity,
//                                     artificialViscosity);
//       ibc = gaussianibc;
    }
    else if (problemString == "gaussiansmooth")
    {
      // smoothed Gaussian problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      // for smoothing:  width of physical domain (unused?)
      vector<Real> widthpp(SpaceDim,1.);
      RealVect width;
      ppgodunov.queryarr("initial_width",widthpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
        {
          width[i] = widthpp[i];
        }

      // domainWidth = physical width of domain
      RealVect domainWidth = RealVect::Unit;
      int minWidthDir = domainWidth.minDir(true);
      Real widthMin = domainWidth[minWidthDir];

      Real size = 0.25;
      ppgodunov.get("initial_size",size);

      Real outerFlatBuffer = 0.;
      Real radmax = 0.;
      int hasRadMaxIB = ppgodunov.query("initial_buffer", outerFlatBuffer);
      int hasRadMax   = ppgodunov.query("initial_max_radius", radmax);
      if (hasRadMax)
        {
          if (hasRadMaxIB)
            {
              pout() << "Warning: initial_max_radius supersedes initial_buffer"
                     << std::endl;
            }
        }
      else
        {
          radmax = (1.0 - outerFlatBuffer) * 0.5 * widthMin;
        }

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        if (hasRadMaxIB)
          {
            pout() << "initial_buffer = " << outerFlatBuffer << endl;
          }
        pout() << "initial_size = " << size << endl;
        pout() << "initial_max_radius = " << radmax << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      GaussianMappedIBC* gaussianibc = new GaussianMappedIBC;
      gaussianibc->setFortranCommon(smallPressure,
                                    gamma,
                                    ambientDensity,
                                    deltaDensity,
                                    pressure,
                                    center,
                                    size,
                                    radmax,
                                    velocity,
                                    artificialViscosity);
      ibc = gaussianibc;
    }
    else if (problemString == "cosinepower")
    {
      // cosine power problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      // domainWidth = physical width of domain
      RealVect domainWidth = RealVect::Unit;
      int minWidthDir = domainWidth.minDir(true);
      Real widthMin = domainWidth[minWidthDir];

      int cosinePower = 6;
      ppgodunov.query("cosine_power", cosinePower);

      Real outerFlatBuffer = 0.;
      ppgodunov.query("initial_buffer", outerFlatBuffer);
      Real radmax = (1.0 - outerFlatBuffer) * 0.5 * widthMin;

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_buffer = " << outerFlatBuffer << endl;
        pout() << "initial_max_radius = " << radmax << endl;
        pout() << "initial_cosine_power = " << cosinePower << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for cosine-power problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       CosinePowerIBC* cospoweribc = new CosinePowerIBC;
//       cospoweribc->setFortranCommon(smallPressure,
//                                     gamma,
//                                     ambientDensity,
//                                     deltaDensity,
//                                     pressure,
//                                     center,
//                                     radmax,
//                                     cosinePower,
//                                     velocity,
//                                     artificialViscosity);
//       ibc = cospoweribc;
    }
    else if (problemString == "cosineeachpower")
    {
      // cosine power problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      // domainWidth = physical width of domain
      RealVect domainWidth = RealVect::Unit;
      int minWidthDir = domainWidth.minDir(true);
      Real widthMin = domainWidth[minWidthDir];

      int cosinePower = 6;
      ppgodunov.query("cosine_power", cosinePower);

      Real outerFlatBuffer = 0.;
      ppgodunov.query("initial_buffer", outerFlatBuffer);
      Real radmax = (1.0 - outerFlatBuffer) * 0.5 * widthMin;

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_buffer = " << outerFlatBuffer << endl;
        pout() << "initial_max_radius = " << radmax << endl;
        pout() << "initial_cosine_power = " << cosinePower << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for cosine-power problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       CosineEachPowerIBC* cospoweribc = new CosineEachPowerIBC;
//       cospoweribc->setFortranCommon(smallPressure,
//                                     gamma,
//                                     ambientDensity,
//                                     deltaDensity,
//                                     pressure,
//                                     center,
//                                     radmax,
//                                     cosinePower,
//                                     velocity,
//                                     artificialViscosity);
//       ibc = cospoweribc;
    }
    else if (problemString == "vortex")
    {
      // smoothed Gaussian problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString
               << "\", specified in input file" << endl << endl;
        return;
      }

      std::string vortexFile;
      ppgodunov.query("vortex_file", vortexFile);

      Vector<RealVect> center;
      Vector<Real> radius;
      Vector<Real> strength;
      setupVortices(center, radius, strength, verbosity, vortexFile);

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        for (int iVort = 0; iVort < center.size(); iVort++)
          {
            RealVect& thisCenter = center[iVort];
            pout() << "vortex_center = " << D_TERM6(thisCenter[0] << "  " <<,
                                                    thisCenter[1] << "  " <<,
                                                    thisCenter[2] << "  " <<,
                                                    thisCenter[3] << "  " <<,
                                                    thisCenter[4] << "  " <<,
                                                    thisCenter[5] << ) endl;
            pout() << "vortex_radius = " << radius[iVort] << endl;
            pout() << "vortex_strength = " << strength[iVort] << endl;
          }
      }

      // Define IBC for vortex problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       VortexIBC* vortexibc = new VortexIBC;
//       vortexibc->setFortranCommon(smallPressure,
//                                   gamma,
//                                   ambientDensity,
//                                   pressure,
//                                   center,
//                                   radius,
//                                   strength,
//                                   artificialViscosity);
//       ibc = vortexibc;
    }
    else if (problemString == "shear")
    {
      // shear-layer problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString
               << "\", specified in input file" << endl << endl;
        return;
      }

      Real widthinv = 42.;
      ppgodunov.get("shear_widthinv", widthinv);

      Real delt = 0.05;
      ppgodunov.get("shear_delt", delt);

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "shear_widthinv = " << widthinv << endl;
        pout() << "shear_delt = " << delt << endl;
      }

      // Define IBC for vortex problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       ShearIBC* shearibc = new ShearIBC;
//       shearibc->setFortranCommon(smallPressure,
//                                  gamma,
//                                  ambientDensity,
//                                  pressure,
//                                  widthinv,
//                                  delt,
//                                  artificialViscosity);
//       ibc = shearibc;
    }
    else if (problemString == "trig")
    {
      // trigonometric problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
      }

      // Define IBC for vortex problem
      TrigMappedIBC* trigibc = new TrigMappedIBC;
      trigibc->setFortranCommon(smallPressure,
                                gamma,
                                ambientDensity,
                                artificialViscosity);
      ibc = trigibc;
    }
    else if (problemString == "gaussian1d")
    {
      // 1D Gaussian problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppgodunov.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       Gaussian1dIBC* gaussianibc = new Gaussian1dIBC;
//       gaussianibc->setFortranCommon(smallPressure,
//                                     gamma,
//                                     ambientDensity,
//                                     deltaDensity,
//                                     pressure,
//                                     center,
//                                     size,
//                                     velocity,
//                                     artificialViscosity);
//       ibc = gaussianibc;
    }
    else if (problemString == "gaussian1dv")
    {
      // 1D Gaussian problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppgodunov.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       Gaussian1dvIBC* gaussianibc = new Gaussian1dvIBC;
//       gaussianibc->setFortranCommon(smallPressure,
//                                     gamma,
//                                     ambientDensity,
//                                     deltaDensity,
//                                     pressure,
//                                     center,
//                                     size,
//                                     velocity,
//                                     artificialViscosity);
//       ibc = gaussianibc;
    }
    else if (problemString == "gaussianp")
    {
      // Gaussian pressure problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real ambientPressure = 1.4;
      ppgodunov.get("ambient_pressure",ambientPressure);

      Real deltaPressure = 0.014;
      ppgodunov.get("delta_pressure",deltaPressure);

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppgodunov.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "ambient_pressure = " << ambientPressure << endl;
        pout() << "delta_pressure = " << deltaPressure << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       GaussianPBC* gaussianpbc = new GaussianPBC;
//       gaussianpbc->setFortranCommon(smallPressure,
//                                     gamma,
//                                     ambientDensity,
//                                     ambientPressure,
//                                     deltaPressure,
//                                     center,
//                                     size,
//                                     velocity,
//                                     artificialViscosity);
//       ibc = gaussianpbc;
    }
    else if (problemString == "sine")
    {
      // Sine problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for sine problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       SineIBC* sineibc = new SineIBC;
//       sineibc->setFortranCommon(smallPressure,
//                                 gamma,
//                                 ambientDensity,
//                                 deltaDensity,
//                                 pressure,
//                                 center,
//                                 velocity,
//                                 artificialViscosity);
//       ibc = sineibc;
    }
    else if (problemString == "wave")
    {
      // Plane wave problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
        {
          pressure = PRESSURE_ISENTROPIC;
        }
      else if (pressureString == "constant")
        {
          pressure = PRESSURE_CONSTANT;
        }

      if (pressure == -1)
        {
          pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
          return;
        }
      // petermc, 7 July 2008:  allow PRESSURE_CONSTANT.
      //      else
      //        {
      //          pressure = PRESSURE_ISENTROPIC;
      //        }

      vector<int> waveNumberpp(SpaceDim,0);
      waveNumberpp[0] = 1;
      IntVect waveNumber;
      ppgodunov.getarr("wave_number",waveNumberpp,0,SpaceDim);
      int norm2 = 0;
      for (int i = 0; i < SpaceDim; i++)
      {
        waveNumber[i] = waveNumberpp[i];
        norm2 += waveNumber[i]*waveNumber[i];
      }

      if (norm2 == 0)
      {
        pout() << "One component of the wave number must be non-zero" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "wave_number = " << D_TERM6(waveNumber[0] << "  " <<,
                                              waveNumber[1] << "  " <<,
                                              waveNumber[2] << "  " <<,
                                              waveNumber[3] << "  " <<,
                                              waveNumber[4] << "  " <<,
                                              waveNumber[5] << ) endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for plane wave problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       WaveIBC* waveibc = new WaveIBC;
//       waveibc->setFortranCommon(smallPressure,
//                                 gamma,
//                                 ambientDensity,
//                                 deltaDensity,
//                                 pressure,
//                                 waveNumber,
//                                 center,
//                                 velocity,
//                                 artificialViscosity);
//       ibc = waveibc;
    }
    else if (problemString == "source")
    {
      // Plane wave problem with a source
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
        {
          pressure = PRESSURE_ISENTROPIC;
        }
      else if (pressureString == "constant")
        {
          pressure = PRESSURE_CONSTANT;
        }

      useSourceTerm = true;
      ppgodunov.get("source_term_scaling",sourceTermScaling);

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }
      else
      {
        pressure = PRESSURE_ISENTROPIC;
      }

      vector<int> waveNumberpp(SpaceDim,0);
      waveNumberpp[0] = 1;
      IntVect waveNumber;
      ppgodunov.getarr("wave_number",waveNumberpp,0,SpaceDim);
      int norm2 = 0;
      for (int i = 0; i < SpaceDim; i++)
      {
        waveNumber[i] = waveNumberpp[i];
        norm2 += waveNumber[i]*waveNumber[i];
      }

      if (norm2 == 0)
      {
        pout() << "One component of the wave number must be non-zero" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "wave_number = " << D_TERM6(waveNumber[0] << "  " <<,
                                              waveNumber[1] << "  " <<,
                                              waveNumber[2] << "  " <<,
                                              waveNumber[3] << "  " <<,
                                              waveNumber[4] << "  " <<,
                                              waveNumber[5] << ) endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                 center[1] << "  " <<,
                                                 center[2] << "  " <<,
                                                 center[3] << "  " <<,
                                                 center[4] << "  " <<,
                                                 center[5] << ) endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                   velocity[1] << "  " <<,
                                                   velocity[2] << "  " <<,
                                                   velocity[3] << "  " <<,
                                                   velocity[4] << "  " <<,
                                                   velocity[5] << ) endl;
      }

      // Define IBC for plane wave with source problem
      pout() << problemString << std::endl;
      MayDay::Error("This problem does not have a mapped IBC file");
//       SourceIBC* sourceibc = new SourceIBC;
//       sourceibc->setFortranCommon(smallPressure,
//                                   gamma,
//                                   ambientDensity,
//                                   deltaDensity,
//                                   pressure,
//                                   waveNumber,
//                                   center,
//                                   velocity,
//                                   artificialViscosity);
//       ibc = sourceibc;
    }
    else
    {
      // The sample problem name given isn't valid
      pout() << "Invalid problem, \"" << problemString << "\", specified in input file" << endl << endl;
      return;
    }
  }
  else
  {
    // A sample problem must be specified
    pout() << "\"godunov.problem\" not specified in input file" << endl << endl;
    return;
  }

/*--------------------------------------------------------------------*
 * Create coordinate system factory and initialize
 *--------------------------------------------------------------------*/

  ParmParse coordSysPP("coordsys");
  string coordSysString("UNDEFINED");
  coordSysPP.query("type", coordSysString);

  if (coordSysString == "cylindrical")
    {
      RealVect stretch(RealVect::Unit);
      if (coordSysPP.contains("stretch"))
        {
          std::vector<Real> tempvect(SpaceDim, 1);
          coordSysPP.getarr("stretch", tempvect, 0, SpaceDim);
          for (int dir=0; dir<SpaceDim; dir++) stretch[dir] = tempvect[dir];
        }
      coordSysFact = new RThetaZCSFactory(stretch);
    }
  else if (coordSysString == "spherical")
    {
      Real stretch = 1.;
      coordSysPP.query("stretch", stretch);
      Real rMin = 1.;
      coordSysPP.query("min_radius", rMin);
      Real phiMin = Pi/4.;
      coordSysPP.query("min_inclination", phiMin);
      coordSysFact = new RThetaPhiCSFactory(stretch, rMin, phiMin);
    }
  else if (coordSysString == "affine")
    {
      RealVect transVect;
      std::vector<Real> b( SpaceDim, 0.0 );
      coordSysPP.getarr("translation_vector", b, 0, SpaceDim );
      D_TERM6(transVect[0]=b[0];,
              transVect[1]=b[1];,
              transVect[1]=b[2];,
              transVect[1]=b[3];,
              transVect[1]=b[4];,
              transVect[2]=b[5];)

        Vector<RealVect> transMtrx( SpaceDim );
      std::vector<Real> A( SpaceDim * SpaceDim );
      coordSysPP.getarr("transformation_matrix", A, 0, SpaceDim * SpaceDim );
      for (int j=0; j<SpaceDim; j++)
        for (int i=0; i<SpaceDim; i++)
          transMtrx[i][j] = A[i+j*SpaceDim];
      MayDay::Error("AffineCoordSys not implemented right now");
#if 0
      coordSysFact = new AffineCSFactory(transVect,
                                         transMtrx);
#endif
    }
  else if (coordSysString == "twisted")
    {

      Real radius;
      coordSysPP.get("radius", radius );
      //CH_assert(2*radius<=domainLength[domainLength.maxDir(false)]);
      CH_assert(radius>0);

      Real twist;
      coordSysPP.get("twist", twist );

      coordSysFact = new TwistedCSFactory(radius,
                                          twist);
    }
  else if (coordSysString == "warped")
    {

      RealVect scale;
      std::vector<Real> b( SpaceDim, 0.0 );
      coordSysPP.getarr("scale", b, 0, SpaceDim );
      D_TERM6(scale[0]=b[0];,
              scale[1]=b[1];,
              scale[1]=b[2];,
              scale[1]=b[3];,
              scale[1]=b[4];,
              scale[2]=b[5];)

        Real rtol = RTOL;
      coordSysPP.query( "relative_tolerance", rtol );
      CH_assert(rtol>0);

      Real atol = ATOL;
      coordSysPP.query( "absolute_tolerance", atol );
      CH_assert(atol>0);

      int imax = 100;
      coordSysPP.query( "maximum_iterations", imax );
      CH_assert(imax>0);

      coordSysFact = new WarpedCSFactory(scale,
                                         rtol,
                                         atol,
                                         imax);
    }
  else if (coordSysString == "cartesian")
    {
      RealVect origin(RealVect::Zero);
      if (coordSysPP.contains("origin"))
        {
          std::vector<Real> tempvect(SpaceDim, 0);
          coordSysPP.getarr("origin", tempvect, 0, SpaceDim);
          for (int dir=0; dir<SpaceDim; dir++) origin[dir] = tempvect[dir];
        }
      RealVect stretch(RealVect::Unit);
      if (coordSysPP.contains("stretch"))
        {
          std::vector<Real> tempvect(SpaceDim, 1);
          coordSysPP.getarr("stretch", tempvect, 0, SpaceDim);
          for (int dir=0; dir<SpaceDim; dir++) stretch[dir] = tempvect[dir];
        }

      coordSysFact = new CartesianCSFactory(origin, stretch);
    }
  // The AMRLevelFactory takes responsibility for destruction of coordSysFact

  // Set up output files
  std::string plotPrefix;
  if (ppgodunov.contains("plot_prefix"))
  {
    ppgodunov.query("plot_prefix",plotPrefix);
  }

  std::string chkPrefix;
  if (ppgodunov.contains("checkpoint_prefix"))
  {
    ppgodunov.query("checkpoint_prefix",chkPrefix);
  }

  if (verbosity >= 2)
  {
    pout() << "verbosity = " << verbosity << endl;

    pout() << "maximum_step = " << nstop << endl;
    pout() << "maximum_time = " << stopTime << endl;
    if (fixedDt > 0)
    {
      pout() << "fixed_dt = " << fixedDt << endl;
    }

    pout() << "number_of_cells = " << D_TERM6(numCells[0] << "  " <<,
                                              numCells[1] << "  " <<,
                                              numCells[2] << "  " <<,
                                              numCells[3] << "  " <<,
                                              numCells[4] << "  " <<,
                                              numCells[5] << ) endl;
    pout() << "is_period = " << D_TERM6(isPeriodic[0] << "  " <<,
                                        isPeriodic[1] << "  " <<,
                                        isPeriodic[2] << "  " <<,
                                        isPeriodic[3] << "  " <<,
                                        isPeriodic[4] << "  " <<,
                                        isPeriodic[5] << ) endl;

    pout() << "maximum_level = " << maxLevel << endl;
    pout() << "refinement_ratio = ";
    for (int i = 0; i < refRatios.size(); ++i)
    {
      pout() << refRatios[i] << " ";
    }
    pout() << endl;

    pout() << "regrid_interval = ";
    for (int i = 0; i < regridIntervals.size(); ++i)
    {
      pout() << regridIntervals[i] << " ";
    }
    pout() << endl;
    pout() << "tag_buffer_size = " << tagBufferSize << endl;
    pout() << "tag_strategy = " << tagStrategyName << endl;

    pout() << "refinement_threshold = " << refineThresh << endl;

    pout() << "blocking_factor = " << blockFactor << endl;
    pout() << "max_grid_size = " << maxGridSize << endl;
    pout() << "fill_ratio = " << fillRatio << endl;
    pout() << "grid buffer size = " << gridBufferSize << endl;

    pout() << "normal_predictor = ";
    if (normalPredOrder == 1)
    {
      pout() << "PLM" << endl;
    }
    else if (normalPredOrder == 2)
    {
      pout() << "PPM" << endl;
    }
    else
    {
      pout() << "Unknown (" << normalPredOrder << ")" << endl;
    }

    pout() << "slope_order = "
           << (useFourthOrderSlopes ? "4th" : "2nd") << endl;
    pout() << "use_primitive_slope_limiting = "
           << (usePrimLimiting ? "yes" : "no") << endl;
    pout() << "use_characteristic_slope_limiting = "
           << (useCharLimiting ? "yes" : "no") << endl;
    pout() << "initial_average = "
           << (initialAverage ? "yes" : "no") << endl;
    pout() << "use_slope_flattening = "
           << (useFlattening ? "yes" : "no") << endl;

    pout() << "use_artificial_viscosity = "
           << (useArtificialViscosity ? "yes" : "no") << endl;
    if (useArtificialViscosity)
    {
      pout() << "artificial_viscosity = " << artificialViscosity << endl;
    }

    pout() << "use_source_term = "
           << (useSourceTerm ? "yes" : "no") << endl;
    if (useSourceTerm)
    {
      pout() << "source_term_scaling = " << sourceTermScaling << endl;
    }

    pout() << "checkpoint_interval = " << checkpointInterval << endl;
    pout() << "plot_interval = " << plotInterval << endl;

    pout() << "CFL = " << cfl << endl;
    pout() << "initial_CFL = " << initialCFL << endl;

    pout() << "maximum_dt_growth = " << maxDtGrowth << endl;
    pout() << "dt_tolerance_factor = " << dtToleranceFactor << endl;
  }

  // Set up the physics for polytropic gas dynamics
  Real M0sq=-1.;
  ppgodunov.query("fourth_order_artificial_viscosity_parameter",M0sq);
  MOLPolytropicPhysics polytropicPhysics(smallPressure);
  polytropicPhysics.setPhysIBC(ibc);
  if (M0sq > 0.)
  polytropicPhysics.setFourthOrderArtificialViscosityParameter(M0sq);
  if (verbosity >= 2)
    {
      pout() << "fourth-order artificial viscosity parameters = " << M0sq
             << endl;
    }
  // Cast to physics base class pointer for technical reasons
  MOLPhysics* molPhysics = static_cast<MOLPhysics*> (&polytropicPhysics);

  int spaceOrder = (useFourthOrderSlopes) ? 4 : 2;

  // Strategies for stability, tagging.
  AMRLevelMappedStabilityStrategy* stabilityStrategy =
    new EulerEquationsMappedStabilityStrategy(1.3925);
  AMRLevelMappedTaggingStrategy* taggingStrategy;
  if (tagStrategyName == "gaussian_analytic")
    {
      taggingStrategy = new GaussianAnalyticMappedTaggingStrategy();
    }
  else
    {
      taggingStrategy = new DensityGradientMappedTaggingStrategy(
        refineThresh, refinementIsScaled, tagPressure, tagVorticity);
    }

  // Set up the AMRLevel... factory
  AMRLevelMappedConsFactory amrGodFact(coordSysFact, stabilityStrategy,
                                       taggingStrategy, plotPrefix);
  amrGodFact.spaceOrder(spaceOrder);
  amrGodFact.limitFaceValues(usePrimLimiting); // WAS limitFaceValues
  amrGodFact.highOrderLimiter(highOrderLimiter); // WAS limitFaceValues
  amrGodFact.initialAverage(initialAverage);
  amrGodFact.useFlattening(useFlattening);
  amrGodFact.noPPM(noPPM);
  amrGodFact.doDeconvolution(doDeconvolution);
  amrGodFact.doFaceDeconvolution(doFaceDeconvolution);
  amrGodFact.useArtificialViscosity(useArtificialViscosity);
  amrGodFact.artificialViscosity(artificialViscosity);
  amrGodFact.useArtVisc(useArtVisc);
  amrGodFact.ratioArtVisc(ratioArtVisc);
  amrGodFact.forwardEuler(forwardEuler);
  // amrGodFact.enforceMinVal(enforceMinVal, minVal);
  amrGodFact.CFL(cfl);
  amrGodFact.domainLength(domainLength);
  amrGodFact.refinementThreshold(refineThresh);
  amrGodFact.refinementIsScaled(refinementIsScaled);
  amrGodFact.tagPressure(tagPressure);
  amrGodFact.tagVorticity(tagVorticity);
  amrGodFact.tagBufferSize(tagBufferSize);
  amrGodFact.verbosity(verbosity);
  amrGodFact.initialDtMultiplier(initialCFL);
  // amrGodFact.IBC(ibc);
  amrGodFact.molPhysics(molPhysics);
  // molPhysics,
  // normalPredOrder,
  // useCharLimiting,
  // useFlattening,
  // useArtificialViscosity,
  // artificialViscosity,
  // useSourceTerm,
  // sourceTermScaling);

  { // scope of AMR amr;
  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel, refRatios, probDomain, &amrGodFact);

  if (fixedDt > 0)
  {
    amr.fixedDt(fixedDt);
  }

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  amr.blockFactor(blockFactor);
  amr.fillRatio(fillRatio);

  amr.gridBufferSize(gridBufferSize);

  // Set output parameters
  amr.checkpointInterval(checkpointInterval);
  amr.plotInterval(plotInterval);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(maxDtGrowth);
  amr.dtToleranceFactor(dtToleranceFactor);
  amr.useSubcyclingInTime(useSubcycling);

  // Set up output files
  if (plotPrefix.size() != 0)
    {
      amr.plotPrefix(plotPrefix);
    }

  if (chkPrefix.size() != 0)
    {
      amr.checkpointPrefix(chkPrefix);
    }

  amr.verbosity(verbosity);

  // Set up input files
  if (!ppgodunov.contains("restart_file"))
  {
    if (!ppgodunov.contains("fixed_hierarchy"))
    {
      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      amr.setupForNewAMRRun();
    }
    else
    {
      //      std::string gridFile;
      //      ppgodunov.query("fixed_hierarchy",gridFile);

      // initialize from a list of grids in "gridFile"
      int numLevels = maxLevel+1;
      Vector<Vector<Box> > amrGrids(numLevels);
//       setupFixedGrids(amrGrids,
//                       probDomain,
//                       maxLevel,
//                       maxGridSize,
//                       blockFactor,
//                       verbosity,
//                       gridFile);
      GenFuncs::readBoxes(amrGrids, ppgodunov, probDomain,
                          maxGridSize, blockFactor,
                          numLevels, refRatios, verbosity);
      amr.setupForFixedHierarchyRun(amrGrids,1);
    }
  }
  else // read from restart file
  {
    std::string restartFile;
    ppgodunov.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrRun restart only defined with hdf5");
#endif
    // This section added by petermc, 20 Jan 2010
    if (ppgodunov.contains("fixed_hierarchy"))
      {
        int numLevels = maxLevel+1;
        Vector<int> regridInterval(numLevels);
        for (int ilev = 0; ilev < numLevels; ilev++) regridInterval[ilev] = 0;
        amr.regridIntervals(regridInterval);
      }
  }

  // Run the computation
  amr.run(stopTime,nstop);

  // Output the last plot file and statistics
  amr.conclude();

  } // scope of amr
}

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile)
{
  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
  {
    a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));

    // Read in predefined grids
    ifstream is(a_gridFile.c_str(), ios::in);

    if (is.fail())
    {
      MayDay::Error("Cannot open grids file");
    }

    // Format of file:
    //   number of levels, then for each level (starting with level 1):
    //   number of grids on level, list of boxes

    int inNumLevels;
    is >> inNumLevels;

    CH_assert (inNumLevels <= a_maxLevel+1);

    if (a_verbosity >= 3)
    {
      pout() << "numLevels = " << inNumLevels << endl;
    }

    while (is.get() != '\n');

    a_amrGrids.resize(inNumLevels);

    // Check to see if coarsest level needs to be broken up
    domainSplit(a_domain,a_amrGrids[0],a_maxGridSize,a_blockFactor);

    if (a_verbosity >= 3)
    {
      pout() << "level 0: ";
      for (int n = 0; n < a_amrGrids[0].size(); n++)
      {
        pout() << a_amrGrids[0][0] << endl;
      }
    }

    // Now loop over levels, starting with level 1
    int ngrid;
    for (int lev = 1; lev < inNumLevels; lev++)
    {
      is >> ngrid;

      if (a_verbosity >= 3)
      {
        pout() << "level " << lev << " numGrids = " << ngrid << endl;
        pout() << "Grids: ";
      }

      while (is.get() != '\n');

      a_amrGrids[lev].resize(ngrid);

      for (int i = 0; i < ngrid; i++)
      {
        Box bx;
        is >> bx;

        while (is.get() != '\n');

        // Quick check on box size
        Box bxRef(bx);

        if (bxRef.longside() > a_maxGridSize)
        {
          pout() << "Grid " << bx << " too large" << endl;
          MayDay::Error();
        }

        if (a_verbosity >= 3)
        {
          pout() << bx << endl;
        }

        a_amrGrids[lev][i] = bx;
      } // End loop over boxes on this level
    } // End loop over levels
  }

  // Broadcast results to all the processors
  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}

// setupVortices allows vortex parameters to be read in and used in this AMR
// computation example
void setupVortices(Vector<RealVect>& a_center,
                   Vector<Real>&     a_radius,
                   Vector<Real>&     a_strength,
                   int               a_verbosity,
                   std::string       a_vortexFile)
{
  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
  {
    // Read in predefined grids
    ifstream is(a_vortexFile.c_str(), ios::in);

    if (is.fail())
      {
        MayDay::Error("Cannot open vortex file");
      }

    // Format of file:
    // number of vortices
    // (then for each vortex)
    // x y z of center
    // radius of vortex from x,y,z=0
    // strength

    int inNumVortices;
    // while (is.get() != '\n');
    is >> inNumVortices;

    a_center.resize(inNumVortices);
    a_radius.resize(inNumVortices);
    a_strength.resize(inNumVortices);

    // Now loop over vortices
    for (int iVort = 0; iVort < inNumVortices; iVort++)
      {
        Real x, y, z;
        //        while (is.get() != '\n');
        is >> x;
        is >> y;
        is >> z;
        a_center[iVort] = RealVect(D_DECL6(x, y, z, 0., 0., 0.));

        Real radius;
        //        while (is.get() != '\n');
        is >> radius;
        a_radius[iVort] = radius;

        Real strength;
        //        while (is.get() != '\n');
        is >> strength;
        a_strength[iVort] = strength;
      } // End loop over vortices
  }

  // Broadcast results to all the processors
  D_TERM6(
    broadcast(a_center[0], uniqueProc(SerialTask::compute)); ,
    broadcast(a_center[1], uniqueProc(SerialTask::compute)); ,
    broadcast(a_center[2], uniqueProc(SerialTask::compute)); ,
    broadcast(a_center[3], uniqueProc(SerialTask::compute)); ,
    broadcast(a_center[4], uniqueProc(SerialTask::compute)); ,
    broadcast(a_center[5], uniqueProc(SerialTask::compute)); );
  broadcast(a_radius, uniqueProc(SerialTask::compute));
  broadcast(a_strength, uniqueProc(SerialTask::compute));
}

#ifdef TRAP_FPE
#include <fenv.h>

// FE_INEXACT    inexact result
// FE_DIVBYZERO  division by zero
// FE_UNDERFLOW  result not representable due to underflow
// FE_OVERFLOW   result not representable due to overflow
// FE_INVALID    invalid operation

static void enableFpExceptions ()
{
  if (feclearexcept(FE_ALL_EXCEPT) != 0)
  {
    MayDay::Abort("feclearexcept failed");
  }

  int flags = FE_DIVBYZERO |
              FE_INVALID   |
//              FE_UNDERFLOW |
              FE_OVERFLOW  ;

  if (feenableexcept(flags) == -1)
  {
    MayDay::Abort("feenableexcept failed");
  }
}
#endif
