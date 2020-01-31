#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "AMR.H"

#include "ParmParse.H"
#include "parstream.H"
#include "AMRLevelClaw.H"
#include "AMRLevelClawFactory.H"
#include "ClawPatch.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

//// #define TRAP_FPE
//#ifdef TRAP_FPE
//extern "C" {
//#include <fpu_control.h>
//}
///* IM: Invalid operation mask
// * DM: Denormalized operand mask
// * ZM: Zero-divide mask
// * OM: Overflow mask
// * UM: Underflow mask
// * PM: Precision (inexact result) mask
//  ---(pm is kinda stupid)
//*/
//static void __attribute__ ((constructor)) trapfpe(void)
//{
//  fpu_control_t cw =
//    _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
//  _FPU_SETCW(cw);
//}
//#endif

void amrClaw();

// amrClaw is a function (as opposed to inline in main()) to get
// around MPI scoping problems
// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
  // setChomboMPIErrorHandler();
#endif

  // Check for an input file
  char* inFile = NULL;

#ifdef TRAP_FPE
  //  trapfpe();
#endif

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrClaw...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  // Run amrClaw, i.e., do the computation
  amrClaw();

  pout() << " amrClaw Finished!" << endl;
#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

void amrClaw()
{
  // Read inputs that are prefixed with "claw."
  ParmParse ppclaw("claw");

  // Determine the sample problem specified
  //int problem = -1;
  std::string problemString;

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppclaw.query("verbosity",verbosity);
 CH_assert(verbosity >= 0);

  // Parameters specific to different sample problems

  // Stop after this number of steps
  int nstop = 0;
  ppclaw.query("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppclaw.query("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  // Real domainLength = 1.0;
  // ppclaw.query("domain_length",domainLength);

  // Maximum AMR level limit
  int maxLevel = 0;
  ppclaw.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppclaw.queryarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppclaw.queryarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppclaw.query("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppclaw.query ("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppclaw.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppclaw.query("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppclaw.query("fill_ratio",fillRatio);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppclaw.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppclaw.query("plot_interval",plotInterval);

  // CFL multiplier
  Real cfl = 0.9;
  ppclaw.query("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  if (ppclaw.query("initial_cfl",initialCFL))
  {
      pout() << "WARNING : initial_cfl is not used by ChomboClaw" << endl;
  }

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppclaw.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 100;
  ppclaw.query("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppclaw.query("dt_tolerance_factor",dtToleranceFactor);

  // Now get clawpack parameters
  // RefCountedPtr<ClawPatch> clawPatchFact = new ClawPatch;
  ClawPatch clawPatch;
  clawPatch.get_inputParams(); // Read more variables from input file

  // figure out what to do with initial_dt
  Real initial_dt = clawPatch.get_initial_dt();
  if (fixedDt > 0)
  {
      if (fabs(fixedDt - initial_dt) > 1e-8)
      {
          pout() << "Setting claw.initial_dt = claw.fixed_dt" << endl;
          clawPatch.set_initial_dt(fixedDt);
      }
  }
  else
  {
      if (initial_dt < 0)
      {
          MayDay::Error("You must set either claw.initial_dt or claw.fixed_dt");
      }
  }



  // Set the resolution of the coarsest level
  Vector<int> numCells(SpaceDim);

  numCells[0] = clawPatch.get_mx();
  numCells[1] = clawPatch.get_my();
#if CH_SPACEDIM == 3
  numCells[2] = clawPatch.get_mz();
#endif

 CH_assert(D_TERM(   (numCells[0] > 0),
                   && (numCells[1] > 0),
                   && (numCells[2] > 0)));
 CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                   && (numCells[1] % 2 == 0),
                   && (numCells[2] % 2 == 0)));


  // Get domain length
  Real a[SpaceDim], b[SpaceDim];
  a[0] = clawPatch.get_xlower();
  a[1] = clawPatch.get_ylower();
  b[0] = clawPatch.get_xupper();
  b[1] = clawPatch.get_yupper();
#if CH_SPACEDIM == 3
  a[2] = clawPatch.get_zlower();
  b[2] = clawPatch.get_zupper();
#endif

  int maxcells = 0;
  int maxdim;
  for (int i = 0; i < SpaceDim; i++)
  {
      if (numCells[i] > maxcells)
      {
          maxcells = numCells[i];
          maxdim = i;
      }
  }
  Real domainLength = b[maxdim] - a[maxdim];
  Real h = domainLength/maxcells;

  for (int i = 0; i < SpaceDim; i++)
  {
      Real hi = (b[i] - a[i])/numCells[i];
      if (fabs(hi - h) > 1e-12)
      {
          MayDay::Error("You must have (dx=dy=dz) in computational space");
      }
  }


  const int* mthbc = clawPatch.get_mthbc();
  bool isPeriodic[SpaceDim];
  for (int i = 0; i < SpaceDim; i++)
  {
      isPeriodic[i] = mthbc[2*i] == 2 && mthbc[2*i+1] == 2;
  }


  // Print the parameters

  pout() << "maximum step = " << nstop << endl;
  pout() << "maximum time = " << stopTime << endl;

  pout() << "number of cells = " << D_TERM(numCells[0] << "  " <<,
                                               numCells[1] << "  " <<,
                                               numCells[2] << ) endl;

  pout() << "maximum level = " << maxLevel << endl;

  pout() << "refinement ratio = ";
  for (int i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
  pout() << endl;

  pout() << "regrid interval = ";
  for (int i = 0; i < regridIntervals.size(); ++i) pout() << regridIntervals[i] << " ";
  pout() << endl;

  pout() << "refinement threshold = " << refineThresh << endl;

  pout() << "blocking factor = " << blockFactor << endl;
  pout() << "max grid size = " << maxGridSize << endl;
  pout() << "fill ratio = " << fillRatio << endl;

  pout() << "checkpoint interval = " << checkpointInterval << endl;
  pout() << "plot interval = " << plotInterval << endl;
  pout() << "desired CFL = " << cfl << endl;
  // pout() << "initial CFL = " << initialCFL << endl;
  if (fixedDt > 0)
    {
      pout() << "fixed dt = " << fixedDt << endl;
    }
  // pout() << "maximum dt growth = " << maxDtGrowth << endl;
  // pout() << "dt tolerance factor = " << dtToleranceFactor << endl;


  // print out clawpack related stuff
  clawPatch.print_inputParams();

  // Set up problem data as well for clawpack
  setprob_();

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  // Set up the AMRLevel... factory
  //
  AMRLevelClawFactory amrClawFact;
  amrClawFact.CFL(cfl);
  amrClawFact.domainLength(domainLength);
  amrClawFact.refinementThreshold(refineThresh);
  amrClawFact.tagBufferSize(tagBufferSize);
  amrClawFact.verbosity(verbosity);
  amrClawFact.initialDtMultiplier(initialCFL);
  amrClawFact.clawPatch(clawPatch);
  amrClawFact.set_stateNames(clawPatch.get_stateNames());

  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrClawFact);

  if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  amr.blockFactor(blockFactor);
  amr.fillRatio(fillRatio);

  // The hyperbolic codes use a grid buffer of 1
  amr.gridBufferSize(1);

  // Set output parameters
  amr.checkpointInterval(checkpointInterval);
  amr.plotInterval(plotInterval);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(maxDtGrowth);
  amr.dtToleranceFactor(dtToleranceFactor);

  // Set up output files
  if (ppclaw.contains("plot_prefix"))
    {
      std::string prefix;
      ppclaw.query("plot_prefix",prefix);
      amr.plotPrefix(prefix);
    }

  if (ppclaw.contains("chk_prefix"))
    {
      std::string prefix;
      ppclaw.query("chk_prefix",prefix);
      amr.checkpointPrefix(prefix);
    }

  amr.verbosity(verbosity);

  // Set up input files
  if (!ppclaw.contains("restart_file"))
    {
      if (!ppclaw.contains("fixed_hierarchy"))
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
          amr.setupForNewAMRRun();
        }
      else
        {
          std::string gridFile;
          ppclaw.query("fixed_hierarchy",gridFile);

          // initialize from a list of grids in "gridFile"
          Vector<Vector<Box> > amrGrids(maxLevel+1);
          setupFixedGrids(amrGrids,
                          probDomain,
                          maxLevel,
                          maxGridSize,
                          blockFactor,
                          verbosity,
                          gridFile);
          amr.setupForFixedHierarchyRun(amrGrids,1);
        }
    }
  else
    {
      std::string restartFile;
      ppclaw.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
      HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      amr.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrClaw restart only defined with hdf5");
#endif
    }

  // Run the computation
  amr.run(stopTime,nstop);

  // Output the last plot file and statistics
  amr.conclude();
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
//  // Run this task on one processor
//  if (procID() == uniqueProc(SerialTask::compute))
//  {
//    a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));
//
//    // Read in predefined grids
//    ifstream is(a_gridFile.c_str(), ios::in);
//
//    if (is.fail())
//    {
//      MayDay::Error("Cannot open grids file");
//    }
//
//    // Format of file:
//    //   number of levels, then for each level (starting with level 1):
//    //   number of grids on level, list of boxes
//
//    int inNumLevels;
//    is >> inNumLevels;
//
//   CH_assert (inNumLevels <= a_maxLevel+1);
//
//    if (a_verbosity >= 3)
//    {
//      pout() << "numLevels = " << inNumLevels << endl;
//    }
//
//    while (is.get() != '\n');
//
//    a_amrGrids.resize(inNumLevels);
//
//    // Check to see if coarsest level needs to be broken up
//    domainSplit(a_domain,a_amrGrids[0],a_maxGridSize,a_blockFactor);
//
//    if (a_verbosity >= 3)
//    {
//      pout() << "level 0: ";
//      for (int n = 0; n < a_amrGrids[0].size(); n++)
//      {
//        pout() << a_amrGrids[0][0] << endl;
//      }
//    }
//
//    // Now loop over levels, starting with level 1
//    int ngrid;
//    for (int lev = 1; lev < inNumLevels; lev++)
//    {
//      is >> ngrid;
//
//      if (a_verbosity >= 3)
//      {
//        pout() << "level " << lev << " numGrids = " << ngrid << endl;
//        pout() << "Grids: ";
//      }
//
//      while (is.get() != '\n');
//
//      a_amrGrids[lev].resize(ngrid);
//
//      for (int i = 0; i < ngrid; i++)
//      {
//        Box bx;
//        is >> bx;
//
//        while (is.get() != '\n');
//
//        // Quick check on box size
//        Box bxRef(bx);
//
//        if (bxRef.longside() > a_maxGridSize)
//        {
//          pout() << "Grid " << bx << " too large" << endl;
//          MayDay::Error();
//        }
//
//        if (a_verbosity >= 3)
//        {
//          pout() << bx << endl;
//        }
//
//        a_amrGrids[lev][i] = bx;
//      } // End loop over boxes on this level
//    } // End loop over levels
//  }
//
//  // Broadcast results to all the processors
//  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}
