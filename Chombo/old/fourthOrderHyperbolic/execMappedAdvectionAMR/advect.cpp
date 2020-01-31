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
using namespace std;

#include "FABView.H"
// this lets us use dumpIVS, other dump functions
#include "DebugDump.H"

#include "PyParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"

#include "AMR.H"
#include "Scheduler.H"
#include "AMRLevel.H"
#include "AMRLevelMappedConsFactory.H"
#include "AMRLevelMappedCons.H"
#include "LevelGridMetrics.H"
#include "MayDay.H"

// Physics objects for linear advection.
#include "MOLAnalyticAdvectionPhysics.H"
#include "MOLAnalyticLonLatAdvectionPhysics.H"

// We use analytic functions for initial
// conditions and boundary conditions.
#include "AnalyticPeriodicAdvectionIBC.H"
#include "AnalyticAdvectionIBC.H"
#include "AnalyticLonLatAdvectionIBC.H"

// Strategies for stability, tagging.
#include "AnalyticAdvectionMappedStabilityStrategy.H"
#include "AnalyticLonLatAdvectionMappedStabilityStrategy.H"
#include "AdvectionGradientMappedTaggingStrategy.H"
#include "AdvectionAbsoluteMappedTaggingStrategy.H"
#include "AdvectionCosineBellMappedTaggingStrategy.H"

// Coordinate systems.
#include "CartesianCS.H"
#include "AnalyticCS.H"
#include "SingleBlockCSAdaptor.H"
#include "DoubleCartesianCS.H"
#include "TripleCartesianCS.H"
#include "TwistedCS.H"
#include "WarpedCS.H"
#include "CylinderEquiangularCS.H"
#include "CubedSphere2DCS.H"

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
                     string           a_gridFile);

// setupVortices allows vortex parameters to be read in and used in this AMR
// computation example
void setupVortices(Vector<RealVect>& a_center,
                   Vector<Real>&     a_radius,
                   Vector<Real>&     a_strength,
                   int               a_verbosity,
                   string       a_vortexFile);

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
    pout() << "Usage:  advect...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the input file.
  PyParse parms(inFile);

#ifdef TRAP_FPE
  enableFpExceptions();
#endif

  // Run amrRun, i.e., do the computation
  amrRun();

  Real end_memory = get_memory_usage_from_OS();
  pout() << endl
         << "Everything completed --- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << end_memory
         << " MB" << endl;

#if defined(CH_MPI)
  Real avg_memory, min_memory, max_memory;
  gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
amrRun()
{
  // Here's our parser.
  PyParse parser;

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  parser.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Stop after this number of steps
  int nstop = 0;
  parser.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  parser.get("max_time",stopTime);

  // Maximum AMR level limit
  int maxLevel = 0;
  parser.get("max_level", maxLevel);

  // Refinement ratios between levels
  vector<int> refRatios(maxLevel+1, 2);
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  parser.get("ref_ratio",refRatios);
  if (refRatios.size() <= maxLevel)
    refRatios.resize(maxLevel + 1, refRatios.back());

  // Number of coarse time steps from one regridding to the next
  vector<int> regridIntervals(maxLevel+1, 2);
  parser.query("regrid_interval", regridIntervals);
  if (regridIntervals.size() < maxLevel)
    regridIntervals.resize(maxLevel + 1, regridIntervals.back());

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  parser.query("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  parser.query("refine_thresh",refineThresh);

  // Whether refinement threshold is scaled with dx
  int refinementIsScaledInt = 0;
  parser.query("refinement_is_scaled", refinementIsScaledInt);
  bool refinementIsScaled = (refinementIsScaledInt == 1);

  // Whether gradient is relative (vs. absolute)
  int relativeGradientInt = 1;
  parser.query("relative_gradient", relativeGradientInt);
  bool relativeGradient = (relativeGradientInt == 1);

  // Minimum dimension of a grid
  int blockFactor = 1;
  parser.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  parser.query("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  parser.query("fill_ratio",fillRatio);

  // Grid buffer size
  int gridBufferSize = 1;
  const int codeReqBufferSize = LevelGridMetrics::bufferSize4thO(
    refRatios,
    maxLevel,
    5);  // Num ghost 5 is hard coded in AMRLevelMappedCons
  if (parser.contains("grid_buffer_size"))
  {
    parser.get("grid_buffer_size", gridBufferSize);
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
      pout() << endl;
    }
  }
  else
    {
      gridBufferSize = codeReqBufferSize;
    }

  // Order of the normal predictor (CTU -> 0, PLM -> 1, PPM -> 2)
  string normalPred = "PPM";
  int normalPredOrder = 0;
  parser.query("normal_predictor",normalPred);
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
    MayDay::Error("Normal predictor must by PLM or PPM");
  }

  // Use fourth order slopes:  default true
  int inFourthOrderSlopes = 1;
  bool useFourthOrderSlopes;
  parser.query("use_fourth_order_slopes",inFourthOrderSlopes);
  useFourthOrderSlopes = (inFourthOrderSlopes == 1);

  // Arbitrary fiat by petermc, 17 June 2008
  useFourthOrderSlopes = true;

  // Do slope limiting:  default true
  int inPrimLimiting = 1;
  bool usePrimLimiting;
  parser.query("use_prim_limiting",inPrimLimiting);
  usePrimLimiting = (inPrimLimiting == 1);

  // This should actually be 1 even if slope limiting is off
  int highOrderLimiterInt = 1;
  parser.query("high_order_limiter", highOrderLimiterInt);
  bool highOrderLimiter = (highOrderLimiterInt == 1);

  // NEW Kreiss-Oliger artificial viscosity
  int inArtVisc = 0;
  parser.query("use_art_visc", inArtVisc);
  bool useArtVisc = (inArtVisc == 1);

  Real ratioArtVisc = 0.;
  if (useArtVisc)
    {
      parser.query("ratio_art_visc", ratioArtVisc);
    }

  bool forwardEuler = false;
  parser.query("forward_euler", forwardEuler);

  // Do slope limiting using characteristics
  bool useCharLimiting = false;
  parser.query("use_char_limiting",useCharLimiting);

  // Initial values are averaged, since the solution is
  // given pointwise.
  bool initialAverage = true;

  // Do slope flattening:  default false
  bool useFlattening = false;
  parser.query("use_flattening", useFlattening);

  // Avoid PPM:  default false
  bool noPPM = false;
  parser.query("no_ppm", noPPM);

  // Do deconvolution:  default true
  bool doDeconvolution = true;
  parser.query("do_deconvolution", doDeconvolution);

  // Do deconvolution:  default true
  bool doFaceDeconvolution = true;
  parser.query("do_face_deconvolution", doFaceDeconvolution);

  // Apply artificial viscosity based on divergence.
  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.0;
  parser.query("artificial_viscosity",artificialViscosity);
  bool useArtificialViscosity = (artificialViscosity > 0.0);

  // Set up checkpointing
  int checkpointInterval = 0;
  parser.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  parser.query("plot_interval",plotInterval);

  // CFL multiplier
  Real cfl = 0.8;
  parser.get("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  parser.query("initial_cfl",initialCFL);

  // Set up whether to use subcycling in time.
  bool useSubcycling = true;
  parser.query("use_subcycling", useSubcycling);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  parser.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  parser.query("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  parser.query("dt_tolerance_factor",dtToleranceFactor);

  // Don't use source term by default
  bool useSourceTerm = false;

  // Source term multiplier
  Real sourceTermScaling = 0.0;

  // Create and define IBC (initial and boundary condition) object
  RefCountedPtr<ScalarFunction> phi0;
  parser.get("phi0", phi0);

  // Advection velocity field or stream function.
  RefCountedPtr<ScalarFunction> streamFunction;
  parser.query("stream_function", streamFunction);
  RefCountedPtr<VectorFunction> velocity;
  parser.query("velocity", velocity);
  if (velocity.isNull() and streamFunction.isNull())
    MayDay::Error("Either velocity or stream_function must be given.");

  // If we have both a stream function and a velocity function, prefer the
  // former.
  if (!streamFunction.isNull() and !velocity.isNull())
    velocity = RefCountedPtr<VectorFunction>();

  // Coordinate system for mapped grids
  ProblemDomain probDomain;
  MultiBlockCoordSysFactory* coordSysFact = NULL;
  vector<int> numCells(SpaceDim, 0);
  string coordSys;
  parser.query("coord_sys", coordSys);
  bool isPeriodic[SpaceDim];
  MOLPhysics* molPhysics = NULL;
  Real domainLength = 1.0;
  AMRLevelMappedStabilityStrategy* stabilityStrategy = NULL;
  Real stabilityFactor = 1.3925;

  string taggingStrategyString;
  parser.query("tagging_strategy", taggingStrategyString);
  if (taggingStrategyString == "") taggingStrategyString = "gradient";

  if (coordSys == "")
  {
    // Set the resolution of the coarsest level
    parser.get("num_cells",numCells);

    CH_assert(D_TERM(   (numCells[0] > 0),
                     && (numCells[1] > 0),
                     && (numCells[2] > 0)));
    CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                     && (numCells[1] % 2 == 0),
                     && (numCells[2] % 2 == 0)));

    // Set the physical size of the longest dimension of the domain
    parser.query("domain_length",domainLength);

    // For now, everything is periodic.
    for (int dim = 0; dim < SpaceDim; dim++)
    {
      isPeriodic[dim] = 1;
      if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
      {
        pout() << "Using Periodic BCs in direction: " << dim << endl;
      }
    }

    probDomain.define(IntVect::Zero,
                      IntVect(D_DECL(numCells[0]-1,
                      numCells[1]-1,
                      numCells[2]-1)),
                      isPeriodic);

    RealVect origin(RealVect::Zero);
    parser.query("origin", origin);
    RealVect stretch(RealVect::Unit);

    // Retrieve coordinate mappings and Jacobian, if present.
    RefCountedPtr<VectorFunction> X, Xi;
    RefCountedPtr<TensorFunction> J;
    parser.query("X", X);
    parser.query("Xi", Xi);
    parser.query("Jacobian", J);

    // If one of the three is present, all three must be present.
    if (!X.isNull() || !Xi.isNull() || !J.isNull())
    {
      if (X.isNull() || Xi.isNull() || J.isNull())
        MayDay::Error("X, Xi, and Jacobian must all be defined.");
    }

    // If no coordinate mapping is specified, we are using Cartesian
    // coordinates.
    NewCoordSysFactory* factory;
    if (X.isNull())
      factory = new CartesianCSFactory(origin, stretch);
    else
    {
      factory = new AnalyticCSFactory(X, Xi, J);
    }
    coordSysFact = new SingleBlockCSAdaptorFactory(factory);
    PhysIBC* ibc = new AnalyticPeriodicAdvectionIBC(phi0);

    // Set up the physics for advection.
    // FIXME: Doesn't work yet for stream function.
    MOLAnalyticAdvectionPhysics* advectionPhysics =
      new MOLAnalyticAdvectionPhysics(velocity);
    advectionPhysics->setPhysIBC(ibc);

    // Cast to physics base class pointer for technical reasons
    molPhysics = dynamic_cast<MOLPhysics*>(advectionPhysics);

    // This computes a stable time step.
    if (!velocity.isNull())
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, velocity);
    else
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, streamFunction);
  }
  else if ((coordSys == "double_cartesian") ||
           (coordSys == "triple_cartesian"))
  {
    parser.get("num_cells", numCells[0]);
    if (coordSys == "double_cartesian")
      numCells[0] *= 3;
    else
      numCells[0] *= 5;
    for (int d = 1; d < SpaceDim; ++d)
      numCells[d] = numCells[0];

    // Set the logical size of the domain.
    domainLength = (coordSys == "double_cartesian") ? 3.0 : 5.0;

    // This is a periodic domain/coordinate system--it doesn't require
    // explicit treatment of BCs.
    for (int dim = 0; dim < SpaceDim; dim++)
      isPeriodic[dim] = 0;

    probDomain.define(IntVect::Zero,
                      IntVect(D_DECL(numCells[0]-1,
                      numCells[1]-1,
                      numCells[2]-1)),
                      isPeriodic);

    // Set up coordinates and IBCs.
    if (coordSys == "double_cartesian")
      coordSysFact = new DoubleCartesianCSFactory();
    else
      coordSysFact = new TripleCartesianCSFactory();
    PhysIBC* ibc = new AnalyticPeriodicAdvectionIBC(phi0);

    // Set up the physics for advection.
    MOLAnalyticAdvectionPhysics* advectionPhysics =
      new MOLAnalyticAdvectionPhysics(velocity);
    advectionPhysics->setPhysIBC(ibc);

    // Cast to physics base class pointer for technical reasons
    molPhysics = dynamic_cast<MOLPhysics*>(advectionPhysics);

    // This computes a stable time step.
    if (!velocity.isNull())
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, velocity);
    else
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, streamFunction);
  }
  else if (coordSys == "twisted")
  {
    parser.get("num_cells", numCells);

    // Set the physical size of the longest dimension of the domain
    parser.query("domain_length",domainLength);

    // This is a periodic domain/coordinate system--it doesn't require
    // explicit treatment of BCs.
    for (int dim = 0; dim < SpaceDim; dim++)
      isPeriodic[dim] = 0;

    probDomain.define(IntVect::Zero,
                      IntVect(D_DECL(numCells[0]-1,
                      numCells[1]-1,
                      numCells[2]-1)),
                      isPeriodic);

    // Set up coordinates and IBCs.
    Real radius = 1.0, twist = 1.0;
    parser.query("radius", radius);
    parser.query("twist", twist);
    NewCoordSysFactory* factory = new TwistedCSFactory(radius, twist);
    coordSysFact = new SingleBlockCSAdaptorFactory(factory);
    PhysIBC* ibc = new AnalyticPeriodicAdvectionIBC(phi0);

    // Set up the physics for advection.
    MOLAnalyticAdvectionPhysics* advectionPhysics =
      new MOLAnalyticAdvectionPhysics(velocity);
    advectionPhysics->setPhysIBC(ibc);

    // Cast to physics base class pointer for technical reasons
    molPhysics = dynamic_cast<MOLPhysics*>(advectionPhysics);

    // This computes a stable time step.
    if (!velocity.isNull())
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, velocity);
    else
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, streamFunction);
  }
  else if (coordSys == "warped")
  {
    parser.get("num_cells", numCells);

    // Set the physical size of the longest dimension of the domain
    parser.query("domain_length",domainLength);

    // This is a periodic domain/coordinate system--it doesn't require
    // explicit treatment of BCs.
    for (int dim = 0; dim < SpaceDim; dim++)
      isPeriodic[dim] = 0;

    probDomain.define(IntVect::Zero,
                      IntVect(D_DECL(numCells[0]-1,
                      numCells[1]-1,
                      numCells[2]-1)),
                      isPeriodic);

    // Set up coordinates and IBCs.
    RealVect scale;
    parser.query("scale", scale);
    Real relTol = 1e-2, absTol = 1e-2;
    parser.query("relative_tolerance", relTol);
    parser.query("absolute_tolerance", absTol);
    int imax = 100;
    parser.query("maximum_iterations", imax);
    NewCoordSysFactory* factory = new WarpedCSFactory(scale, relTol, absTol, imax);
    coordSysFact = new SingleBlockCSAdaptorFactory(factory);
    PhysIBC* ibc = new AnalyticPeriodicAdvectionIBC(phi0);

    // Set up the physics for advection.
    MOLAnalyticAdvectionPhysics* advectionPhysics =
      new MOLAnalyticAdvectionPhysics(velocity);
    advectionPhysics->setPhysIBC(ibc);

    // Cast to physics base class pointer for technical reasons
    molPhysics = dynamic_cast<MOLPhysics*>(advectionPhysics);

    // This computes a stable time step.
    if (!velocity.isNull())
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, velocity);
    else
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, streamFunction);
  }
  else if (coordSys == "cylinder")
  {
    CH_assert(SpaceDim == 2);
    parser.get("num_cells", numCells[0]);
    numCells[1] = numCells[0];

    // The "domain length," which determines the grid spacing, is
    // made up of unit blocks.
    domainLength = 5.0 * 1.0;

    // This coordinate system isn't periodic.
    for (int dim = 0; dim < SpaceDim; dim++)
      isPeriodic[dim] = 0;

    // The problem domain spans [-2*numCells, 3*numCells-1] in each
    // direction.
    probDomain.define(IntVect(-2*numCells[0], -2*numCells[1]),
                      IntVect(3*numCells[0]-1, 3*numCells[1]-1),
                      isPeriodic);

    CylinderEquiangularCSFactory* cylFact = new CylinderEquiangularCSFactory();
    PhysIBC* ibc = new AnalyticAdvectionIBC(phi0);

    // Set the physical parameters of the cylinder.
    RealVect centerPoint;
    parser.query("center_point", centerPoint);
    cylFact->setCenterPoint(centerPoint);
    RealVect centralRectSize(0.5, 0.5);
    parser.query("central_rectangle_size", centralRectSize);
    cylFact->setCentralRectSize(centralRectSize);
    Real outerRadius = 1.0;
    parser.query("outer_radius", outerRadius);
    cylFact->setOuterRadius(outerRadius);

    // Set the generic coordinate system factory pointer.
    coordSysFact = cylFact;

    // Set up the physics for advection.
    MOLAnalyticAdvectionPhysics* advectionPhysics;
    if (!velocity.isNull())
      advectionPhysics = new MOLAnalyticAdvectionPhysics(velocity);
    else
      advectionPhysics = new MOLAnalyticAdvectionPhysics(streamFunction);
    advectionPhysics->setPhysIBC(ibc);

    // Cast to physics base class pointer for technical reasons
    molPhysics = dynamic_cast<MOLPhysics*>(advectionPhysics);

    // This computes a stable time step.
    if (!velocity.isNull())
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, velocity);
    else
      stabilityStrategy = new AnalyticAdvectionMappedStabilityStrategy(stabilityFactor, streamFunction);
  }
  else if (coordSys == "cubed_sphere")
  {
    CH_assert(SpaceDim == 2);
    parser.get("num_cells", numCells[0]);
    numCells[1] = numCells[0];
    numCells[0] *= 11;

    // The "domain length," which determines the grid spacing, is
    // pi/2 * each panel, including the space in between.
    domainLength = 11.0 * 0.5 * M_PI;

    // On a cubed sphere, the coordinate system isn't periodic
    // in the usual sense.
    for (int dim = 0; dim < SpaceDim; dim++)
      isPeriodic[dim] = 0;

    probDomain.define(IntVect::Zero,
                      IntVect(D_DECL(numCells[0]-1,
                      numCells[1]-1,
                      numCells[2]-1)),
                      isPeriodic);

    coordSysFact = new CubedSphere2DCSFactory();
    PhysIBC* ibc = new AnalyticLonLatAdvectionIBC(phi0);

    // Set up the physics for advection.
    MOLAnalyticLonLatAdvectionPhysics* advectionPhysics;
    if (!velocity.isNull())
      advectionPhysics = new MOLAnalyticLonLatAdvectionPhysics(velocity);
    else
      advectionPhysics = new MOLAnalyticLonLatAdvectionPhysics(streamFunction);
    advectionPhysics->setPhysIBC(ibc);

    // Cast to physics base class pointer for technical reasons
    molPhysics = dynamic_cast<MOLPhysics*>(advectionPhysics);

    // This computes a stable time step.
    if (!velocity.isNull())
      stabilityStrategy = new AnalyticLonLatAdvectionMappedStabilityStrategy(stabilityFactor, velocity);
    else
      stabilityStrategy = new AnalyticLonLatAdvectionMappedStabilityStrategy(stabilityFactor, streamFunction);
  }
  else
  {
    char errMesg[1024];
    snprintf(errMesg, 1024, "Unknown coordinate system: %s", coordSys.c_str());
    MayDay::Error(errMesg);
  }

  // Note: the AMRLevelFactory takes responsibility for destruction
  // of coordSysFact.

  // Set up output files
  string plotPrefix;
  parser.query("plot_prefix", plotPrefix);

  string chkPrefix;
  parser.query("checkpoint_prefix", chkPrefix);

  if (verbosity >= 2)
  {
    pout() << "verbosity = " << verbosity << endl;

    pout() << "maximum_step = " << nstop << endl;
    pout() << "maximum_time = " << stopTime << endl;
    if (fixedDt > 0)
    {
      pout() << "fixed_dt = " << fixedDt << endl;
    }

    pout() << "number_of_cells = " << D_TERM(numCells[0] << "  " <<,
                                             numCells[1] << "  " <<,
                                             numCells[2] << ) endl;
    pout() << "is_period = " << D_TERM(isPeriodic[0] << "  " <<,
                                       isPeriodic[1] << "  " <<,
                                       isPeriodic[2] << ) endl;

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

  int spaceOrder = (useFourthOrderSlopes) ? 4 : 2;

  AMRLevelMappedTaggingStrategy* taggingStrategy;
  if (taggingStrategyString == "gradient")
    {
      taggingStrategy =
        new AdvectionGradientMappedTaggingStrategy(refineThresh,
                                                   refinementIsScaled,
                                                   relativeGradient);
    }
  else if (taggingStrategyString == "absolute")
    {
      taggingStrategy =
        new AdvectionAbsoluteMappedTaggingStrategy(refineThresh);
    }
  else if (taggingStrategyString == "cosinebell")
    {
      Real ambientDensity = 0.0;
      parser.query("ambientDensity", ambientDensity);
      Real deltaDensity = 1000.0;
      parser.query("deltaDensity", deltaDensity);
      Real size = 0.3;
      parser.query("size", size);
      Real advectVelocity = M_PI / 6.0;
      parser.query("advectVelocity", advectVelocity);
      Real advectAngle = 0.0; //M_PI / 2.0;
      parser.query("advectAngle", advectAngle);
      RealVect center(M_PI/2.0, 0.);
      parser.query("center", center);
      taggingStrategy =
        new AdvectionCosineBellMappedTaggingStrategy(
                                                     ambientDensity,
                                                     deltaDensity,
                                                     size,
                                                     advectVelocity,
                                                     advectAngle,
                                                     center,
                                                     refineThresh);
    }
  else
    {
      pout() << "Unknown tagging strategy " << taggingStrategyString << endl;
      MayDay::Error("tagging_strategy must be 'gradient' or 'absolute' or 'cosinebell'");
    }

  // Set up the AMRLevelFactory.
  AMRLevelMappedConsFactory amrFact(coordSysFact, stabilityStrategy, taggingStrategy, plotPrefix);
  amrFact.spaceOrder(spaceOrder);
  amrFact.limitFaceValues(usePrimLimiting); // WAS limitFaceValues
  amrFact.highOrderLimiter(highOrderLimiter); // WAS limitFaceValues
  amrFact.initialAverage(initialAverage);
  amrFact.useFlattening(useFlattening);
  amrFact.noPPM(noPPM);
  amrFact.doDeconvolution(doDeconvolution);
  amrFact.doFaceDeconvolution(doFaceDeconvolution);
  amrFact.useArtificialViscosity(useArtificialViscosity);
  amrFact.artificialViscosity(artificialViscosity);
  amrFact.useArtVisc(useArtVisc);
  amrFact.ratioArtVisc(ratioArtVisc);
  amrFact.forwardEuler(forwardEuler);
  amrFact.CFL(cfl);
  amrFact.domainLength(domainLength);
  amrFact.refinementThreshold(refineThresh);
  amrFact.refinementIsScaled(refinementIsScaled);
  amrFact.tagBufferSize(tagBufferSize);
  amrFact.verbosity(verbosity);
  amrFact.initialDtMultiplier(initialCFL);
  amrFact.molPhysics(molPhysics);

  { // scope of AMR amr
    AMR amr;

    // Set up the AMR object
    amr.define(maxLevel, refRatios, probDomain, &amrFact);

    if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }

    // Set grid generation parameters
    amr.maxGridSize(maxGridSize);
    amr.blockFactor(blockFactor);
    amr.fillRatio(fillRatio);

    amr.gridBufferSize(gridBufferSize);

    // Set up periodic output.
    RefCountedPtr<Scheduler> scheduler(new Scheduler());
    if (plotInterval > 0)
    {
      RefCountedPtr<PlotterPeriodicFunction> plotter(new PlotterPeriodicFunction(plotPrefix));
      scheduler->schedule(plotter, plotInterval);
    }

    if (checkpointInterval > 0)
    {
      string chkPrefix;
      parser.get("chk_prefix", chkPrefix);
      RefCountedPtr<CheckpointPeriodicFunction> dumper(new CheckpointPeriodicFunction(chkPrefix));
      scheduler->schedule(dumper, checkpointInterval);
    }
    amr.schedule(scheduler);

    amr.regridIntervals(regridIntervals);
    amr.maxDtGrow(maxDtGrowth);
    amr.dtToleranceFactor(dtToleranceFactor);
    amr.useSubcyclingInTime(useSubcycling);

    // Set up output files
    if (!plotPrefix.empty())
    {
      amr.plotPrefix(plotPrefix);
    }

    if (!chkPrefix.empty())
    {
      amr.checkpointPrefix(chkPrefix);
    }

    amr.verbosity(verbosity);

    // Set up input files
    if (!parser.contains("restart_file"))
    {
      if (!parser.contains("boxes"))
      {
        // initialize from scratch for AMR run
        // initialize hierarchy of levels
        amr.setupForNewAMRRun();
      }
      else
      {
        // initialize from a set of Boxes.
        int numLevels = maxLevel+1;
        Vector<Vector<Box> > amrGrids(numLevels);
        parser.get("boxes", amrGrids);
        amr.setupForFixedHierarchyRun(amrGrids,1);

        // We don't regrid on fixed hierarchies.
        Vector<int> regridInterval(numLevels);
        for (int ilev = 0; ilev < numLevels; ilev++) regridInterval[ilev] = 0;
        amr.regridIntervals(regridInterval);
      }
    }
    else // read from restart file
    {
      string restartFile;
      parser.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
      HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      amr.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrRun restart only defined with hdf5");
#endif
      // If we have a fixed box layout, we don't regrid.
      if (parser.contains("boxes"))
      {
        int numLevels = maxLevel+1;
        Vector<int> regridInterval(numLevels);
        for (int ilev = 0; ilev < numLevels; ilev++) regridInterval[ilev] = 0;
        amr.regridIntervals(regridInterval);
      }
    }

    pout() << "AMR Setup completed ---- "
           << "mem: "
           << setw(8) << setprecision(3)
           << setiosflags(ios::fixed)
           << get_memory_usage_from_OS()
           << " MB" << endl;

    // Run the computation.
    amr.run(stopTime,nstop);

    pout() << "AMR Run completed ------ "
           << "mem: "
           << setw(8) << setprecision(3)
           << setiosflags(ios::fixed)
           << get_memory_usage_from_OS()
           << " MB" << endl;

    // Output the last plot file and statistics.
    amr.conclude();

    pout() << "AMR Conclude completed - "
           << "mem: "
           << setw(8) << setprecision(3)
           << setiosflags(ios::fixed)
           << get_memory_usage_from_OS()
           << " MB" << endl;
  } // scope of amr
}
//-----------------------------------------------------------------------

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
//-----------------------------------------------------------------------
void
setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                const ProblemDomain&  a_domain,
                int                   a_maxLevel,
                int                   a_maxGridSize,
                int                   a_blockFactor,
                int                   a_verbosity,
                string           a_gridFile)
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
//-----------------------------------------------------------------------

#ifdef TRAP_FPE
#include <fenv.h>

// FE_INEXACT    inexact result
// FE_DIVBYZERO  division by zero
// FE_UNDERFLOW  result not representable due to underflow
// FE_OVERFLOW   result not representable due to overflow
// FE_INVALID    invalid operation

//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
#endif
