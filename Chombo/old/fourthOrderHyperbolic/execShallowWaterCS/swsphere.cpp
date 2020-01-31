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
#include <iomanip>
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

#include "CubedSphere2DCS.H"
#include "MOLShallowWaterPhysics.H"
#include "MOLAdvectionPhysics.H"
#include "ShallowWaterAMRLevelMappedFactory.H"

#include "AdvectionCubedSphereIBC.H"
#include "ShallowWaterCubedSphereIBC.H"

#include "AMR.H"
#include "Scheduler.H"
#include "AMRLevel.H"
#include "AMRLevelMappedConsFactory.H"
#include "AMRLevelAdvectMappedFactory.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

// amrCS is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrCS();

// One more function for MPI
void dumpmemoryatexit();

///////////////////////////////////////////////////////////////////////////////

OldTimer Everything    ("gov Everything", 0);
OldTimer TimeReadInput ("gov Read Input",   Everything);
OldTimer TimeSetupAMR  ("gov Setup AMR",    Everything);
OldTimer TimeRun       ("gov Run",          Everything);
OldTimer TimeConclude  ("gov Conclude",     Everything);

///////////////////////////////////////////////////////////////////////////////

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
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
    pout() << "Usage:  swsphere.ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  PyParse parser(inFile);

  // Run amrCS, i.e., do the computation
  amrCS();

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
amrCS()
{
  // Access the parser.
  PyParse parser;
  CH_assert(parser.parsed());

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  if (parser.contains("verbosity"))
    parser.get("verbosity", verbosity);
  CH_assert(verbosity >= 0);

  // Stop after this number of steps
  int nstop = INT_MAX;
  if (parser.contains("max_step"))
    parser.get("max_step", nstop);

  // Stop when the simulation time get here
  Real stopTime = FLT_MAX;
  if (parser.contains("max_time"))
    parser.get("max_time", stopTime);

  // Set the resolution of the coarsest level
  int numCells;
  parser.get("num_cells", numCells);

  // Maximum AMR level limit
  int maxLevel = 0;
  if (parser.contains("max_level"))
    parser.get("max_level",maxLevel);

  // Refinement ratios between levels
  std::vector<int> refRatios(1, 2);
  if (parser.contains("ref_ratio"))
    parser.get("ref_ratio", refRatios);
  CH_assert(refRatios.size() >= 1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  if (parser.contains("regrid_interval"))
    parser.get("regrid_interval", regridIntervals);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  if (parser.contains("tag_buffer_size"))
    parser.get("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  if (parser.contains("refine_thresh"))
    parser.get("refine_thresh", refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  if (parser.contains("block_factor"))
    parser.get("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  if (parser.contains("max_grid_size"))
    parser.get("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  if (parser.contains("fill_ratio"))
    parser.get("fill_ratio",fillRatio);

  // The hyperbolic codes use a grid buffer of 1
  int gridBufferSize = 1;
  if (parser.contains("grid_buffer_size"))
    parser.get("grid_buffer_size", gridBufferSize);

  // Order of the normal predictor (CTU -> 0, PLM -> 1, PPM -> 2)
  // petermc, 10 Oct 2009, set default PPM
  std::string normalPred = "PPM";
  if (parser.contains("normal_predictor"))
    parser.get("normal_predictor",normalPred);
  int normalPredOrder;
  if (normalPred == "CTU" || normalPred == "ctu")
    normalPredOrder = 0;
  else if (normalPred == "PLM" || normalPred == "plm")
    normalPredOrder = 1;
  else if (normalPred == "PPM" || normalPred == "ppm")
    normalPredOrder = 2;
  else
    MayDay::Error("Normal predictor must be PLM or PPM");

  // Use fourth order slopes:  default true
  bool useFourthOrderSlopes = true;
  if (parser.contains("use_fourth_order_slopes"))
    parser.get("use_fourth_order_slopes",useFourthOrderSlopes);

  // Arbitrary fiat by petermc, 17 June 2008
  useFourthOrderSlopes = true;

  // Do slope limiting:  default true
  bool usePrimLimiting = true;
  if (parser.contains("use_prim_limiting"))
    parser.get("use_prim_limiting", usePrimLimiting);

  // NEW Kreiss-Oliger artificial viscosity
  Real ratioArtVisc = 0.;
  if (parser.contains("ratio_art_visc"))
    parser.get("ratio_art_visc", ratioArtVisc);
  bool useArtVisc = (ratioArtVisc > 0.0);

  bool forwardEuler = false;
  if (parser.contains("forward_euler"))
    parser.get("forward_euler", forwardEuler);

  // Do slope limiting using characteristics.  NOT USED.
  //  int inCharLimiting = 0;
  //  bool useCharLimiting;
  //  parser.get("use_char_limiting",inCharLimiting);
  //  useCharLimiting = (inCharLimiting == 1);

  // Do slope flattening:  default false
  bool useFlattening = false;
  if (parser.contains("use_flattening"))
    parser.get("use_flattening", useFlattening);

  // Avoid PPM:  default false
  bool noPPM = false;
  if (parser.contains("no_ppm"))
    parser.get("no_ppm", noPPM);

  // Do deconvolution:  default true
  bool doDeconvolution = true;
  if (parser.contains("do_deconvolution"))
    parser.get("do_deconvolution", doDeconvolution);

  // Do deconvolution:  default true
  bool doFaceDeconvolution = true;
  if (parser.contains("do_face_deconvolution"))
    parser.get("do_face_deconvolution", doFaceDeconvolution);

  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.0;
  if (parser.contains("artificial_viscosity"))
    parser.get("artificial_viscosity",artificialViscosity);
  bool useArtificialViscosity = (artificialViscosity > 0.0);

  // Set up checkpointing
  int checkpointInterval = 0;
  if (parser.contains("checkpoint_interval"))
    parser.get("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  if (parser.contains("plot_interval"))
    parser.get("plot_interval",plotInterval);

  // CFL multiplier
  Real cfl = 0.8;
  if (parser.contains("cfl"))
    parser.get("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  if (parser.contains("initial_cfl"))
    parser.get("initial_cfl", initialCFL);

  // Set up whether to use subcycling in time.
  bool useSubcycling = false;
  if (parser.contains("use_subcycling"))
    parser.get("use_subcycling", useSubcycling);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  if (parser.contains("fixed_dt"))
    parser.get("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  if (parser.contains("max_dt_growth"))
    parser.get("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  if (parser.contains("dt_tolerance_factor"))
    parser.get("dt_tolerance_factor",dtToleranceFactor);

  // Compute timestep from cell-centered velocities:  default false
//  bool dtFromCells = false;
//  if (parser.contains("dt_from_cells"))
//  parser.get("dt_from_cells", dtFromCells);

  // A minimum pressure needed to construct PolytropicPhysics - used in slope
  // flattening
  Real smallPressure = 1.0e-12;

  // Don't use source term by default
  bool useSourceTerm = false;

  // Source term multiplier
  Real sourceTermScaling = 0.0;

  // Problem domain
  ProblemDomain probDomain(IntVect::Zero,
                           IntVect(numCells, numCells));

  // Create coordinate system factory and initialize
  CubedSphere2DCSFactory coordSysFact;

  // Base class pointers for AMR algorithm, physics, initial conditions.
  AMRLevelMappedConsFactory* factory;
  MOLPhysics* physics;
  PhysMappedIBC* ibc;

  // What are the equations we're solving?
  string physicsStr;
  parser.get("physics", physicsStr);
  if (physicsStr == "shallow_water")
  {
    // We are solving the shallow water equations.
    factory = new ShallowWaterAMRLevelMappedFactory();
    MOLShallowWaterPhysics* swPhysics = new MOLShallowWaterPhysics(smallPressure);
    physics = dynamic_cast<MOLPhysics*>(swPhysics);

    // We need initial conditions for the height and the velocity.
    RefCountedPtr<ScalarFunction> h0;
    parser.get("height", h0);
    RefCountedPtr<VectorFunction> u0;
    parser.get("velocity", u0);

    // Set up the "IBC".
                ShallowWaterCubedSphereIBC* swIBC =
      new ShallowWaterCubedSphereIBC(h0, u0);
    ibc = dynamic_cast<PhysMappedIBC*>(swIBC);
  }
  else if (physicsStr == "advection")
  {
    // Set up the physics for advection.
    factory = new AMRLevelAdvectMappedFactory();
    MOLAdvectionPhysics* advectionPhysics = new MOLAdvectionPhysics();
    physics = dynamic_cast<MOLPhysics*>(advectionPhysics);

    // We need initial conditions for the solution and the velocity.
    RefCountedPtr<ScalarFunction> phi0;
    parser.get("solution", phi0);
    RefCountedPtr<VectorFunction> u0;
    parser.get("velocity", u0);

    // Set up the "IBC".
                AdvectionCubedSphereIBC* advectIBC =
      new AdvectionCubedSphereIBC(phi0, u0);
    ibc = dynamic_cast<PhysMappedIBC*>(advectIBC);
  }
  else
  {
    char err[1024];
    snprintf(err, 1024, "Unknown physics: %s", physicsStr.c_str());
    MayDay::Error(err);
  }

  // Print out some stuff.
  if (verbosity >= 2)
  {
    pout() << "verbosity = " << verbosity << endl;

    pout() << "maximum_step = " << nstop << endl;
    pout() << "maximum_time = " << stopTime << endl;
    if (fixedDt > 0)
    {
      pout() << "fixed_dt = " << fixedDt << endl;
    }

    pout() << "number_of_cells = " << numCells << endl;

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
    //    pout() << "use_characteristic_slope_limiting = "
    //           << (useCharLimiting ? "yes" : "no") << endl;
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

  factory->spaceOrder(spaceOrder);
  factory->limitFaceValues(usePrimLimiting); // WAS limitFaceValues
  factory->initialAverage(true); // Always average initial values.
  factory->useFlattening(useFlattening);
  factory->noPPM(noPPM);
  factory->doDeconvolution(doDeconvolution);
  factory->doFaceDeconvolution(doFaceDeconvolution);
  factory->useArtificialViscosity(useArtificialViscosity);
  factory->artificialViscosity(artificialViscosity);
  factory->useArtVisc(useArtVisc);
  factory->ratioArtVisc(ratioArtVisc);
  factory->forwardEuler(forwardEuler);
  // factory->enforceMinVal(enforceMinVal, minVal);
  factory->CFL(cfl);
  //factory->domainLength(domainLength);
  //factory->refinementThreshold(refineThresh);
  //factory->tagBufferSize(tagBufferSize);
  factory->verbosity(verbosity);
  factory->initialDtMultiplier(initialCFL);
  // factory->IBC(ibc);
  factory->molPhysics(physics);
//  factory->dtFromCells(dtFromCells);
  factory->coordinateSystemFactory(&coordSysFact);
  // swspherePhysics,
  // normalPredOrder,
  // useCharLimiting,
  // useFlattening,
  // useArtificialViscosity,
  // artificialViscosity,
  // useSourceTerm,
  // sourceTermScaling);

  // Set up output files -- do this with AMRLevelFactory
  // in order to be able to do mapped-grid output as well.
  string plotPrefix;
  if (parser.contains("plot_prefix"))
  {
    parser.get("plot_prefix",plotPrefix);
    factory->plotPrefix(plotPrefix);
  }

  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel, refRatios, probDomain, factory);

  if (fixedDt > 0)
    amr.fixedDt(fixedDt);

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

  amr.verbosity(verbosity);

  // Set up input files
  if (!parser.contains("restart_file"))
  {
    // Unimplemented
    CH_assert(false);

/*
    if (!parser.contains("fixed_hierarchy"))
    {
      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      amr.setupForNewAMRRun();
    }
    else
    {
      //      std::string gridFile;
      //      parser.get("fixed_hierarchy",gridFile);

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
      readBoxes(amrGrids, parser, probDomain,
                maxGridSize, blockFactor,
                numLevels, refRatios, verbosity);
      amr.setupForFixedHierarchyRun(amrGrids,1);
    }
*/
  }
  else
  {
    std::string restartFile;
    parser.get("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrCS restart only defined with hdf5");
#endif
  }

  // End timing AMR solver setup
  TimeSetupAMR.stop();

#ifndef CH_NTIMER
  pout() << "AMR Setup completed ---- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage_from_OS()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeSetupAMR.wc_time()
         << " sec (wall-clock)" << endl;

  if (verbosity >= 1)
  {
    pout() << endl;
  }
#endif

  // Run and time the computation
  TimeRun.start();
  amr.run(stopTime,nstop);
  TimeRun.stop();

#ifndef CH_NTIMER
  if (verbosity >= 1)
  {
    pout() << endl;
  }

  pout() << "AMR Run completed ------ "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage_from_OS()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeRun.wc_time()
         << " sec (wall-clock)" << endl;
#endif

  // Output the last plot file and statistics - time the process
  TimeConclude.start();
  amr.conclude();
  TimeConclude.stop();

  delete factory;
  delete physics;
  delete ibc;

#ifndef CH_NTIMER
  pout() << "AMR Conclude completed - "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage_from_OS()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeConclude.wc_time()
         << " sec (wall-clock)" << endl;
#endif
}
//-----------------------------------------------------------------------

/*
// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
*/

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
