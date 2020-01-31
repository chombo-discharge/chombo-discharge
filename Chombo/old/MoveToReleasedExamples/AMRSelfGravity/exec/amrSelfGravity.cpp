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
using std::ifstream;
using std::ios;

#include "ParmParse.H"
#include "parstream.H"

#include "AMRIO.H"
#include "CH_HDF5.H"
#include "MayDay.H"
#include "CH_Timer.H"
#include "memusage.H"

#include "AMR.H"
#include "AMRLevelSelfGravity.H"
#include "AMRLevelSelfGravityFactory.H"

#include "SelfGravityPhysics.H"
#include "RefCellTagger.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "TestsIBC.H"
#include "DustCollapseIBC.H"
#include "RampIBC.H"



// amrSelfGravity is a function to get around MPI scoping problems
void amrSelfGravity();

void setupProblem(SelfGravityPhysics* a_selfGravityPhysics);

//
void setupAMRFactory(AMRLevelSelfGravityFactory& a_amrGodFact
                     ,SelfGravityPhysics*         a_selfGravityPhysics);

//
void setupAMRObject(int& a_maxStep, Real& a_maxTime, AMR& amr,
                    const AMRLevelSelfGravityFactory& amrGodFact);

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

#ifdef CH_Linux
//#define TRAP_FPE  //(should be off by default)
#endif
#ifdef TRAP_FPE
// funproto
static void __attribute__ ((constructor)) trapfpe(void);
#endif

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
  // setChomboMPIErrorHandler();
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

  OldTimer::TimerInit(rank);

#ifdef TRAP_FPE
  trapfpe();
#endif

  // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrSelfGravity...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  // Run amrSelfGravity, i.e., do the computation

  amrSelfGravity();

  Real vmpeak, vmhwm;
  getPeakMemoryFromOS(vmpeak, vmhwm);
  pout() << " getPeakMemoryFromOS(vmpeak, vmhwm)=" << vmpeak << " " << vmhwm << endl;
  reduce_print_avg_min_max("peakRSS", vmpeak);
  reduce_print_avg_min_max("peakVM",  vmhwm);

  Real memtrackCurrentMemory;
  Real memtrackPeakMemory;
  memtrackStamp(memtrackCurrentMemory, memtrackPeakMemory);
  reduce_print_avg_min_max("MT", memtrackPeakMemory);

  print_memory_line("after all");

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

void amrSelfGravity()
{
  // Read inputs that are prefixed with "selfGravity."

  ParmParse ppcharm("charm");

  // max simulation step
  int maxStep = 0;

  // max simulation time
  Real maxTime = 0.0;

  SelfGravityPhysics* selfGravityPhysics = new SelfGravityPhysics();

  // sets up FORTRAN /common/, by reading in all necessary parameters
  // and defines IBC objects.

  setupProblem(selfGravityPhysics);

  // Set up the AMR Factory; read in all the required parameters
  AMRLevelSelfGravityFactory amrGodFact;
  setupAMRFactory(amrGodFact,selfGravityPhysics);

  AMR amr;

  // Set up the AMR object and get the probl. domain back
  setupAMRObject(maxStep,maxTime,amr,amrGodFact);



  // Run the computation

  amr.run(maxTime,maxStep);


  // Output the last plot file and statistics

  amr.conclude();

}


void setupProblem(SelfGravityPhysics*  selfGravityPhysics)
{
  // Read inputs that are prefixed with "selfGravity."
  ParmParse ppcharm("charm");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppcharm.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  /* PROBLEM PARAMETERS */

  /* Hydro */

  // For all gas dynamics
  Real gamma = 5./3.;
  ppcharm.query("gamma",gamma);

  // tolerance for error in Riemann solver
  Real riemannSolTol = 1.e-6;
  ppcharm.query("Rs_tolerance",riemannSolTol);

  //
  int maxRsIter = 10;
  ppcharm.query("max_rs_iter",maxRsIter);

  //  threshold for switch from energy to entropy cons. eq.
  Real maxMach = 5.e1;
  ppcharm.query("max_mach",maxMach);

  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.0;
  ppcharm.query("artificial_viscosity",artificialViscosity);

  // Create initial and boundary condition (IBC) object and initialize
  std::string problemString;
  if (ppcharm.contains("problem"))
  {
    ppcharm.query("problem",problemString);
    if (problemString == "gastests")
    {
      // specify problem further
      std::string testString;
      ppcharm.query("gas_test",testString);

      // For "explosion"
      Real ms = ten;
      ppcharm.query("shock_mach",ms);

      // For "explosion"
      Real size = fourth;
      ppcharm.query("size",size);

      // For "explosion" and "wave"
      vector<Real> centerpp(SpaceDim,zero);
      RealVect center;
      ppcharm.queryarr("center",centerpp,0,SpaceDim);
      for (int i=0; i<SpaceDim; i++) center[i]=centerpp[i];

      vector<Real> velocitypp(SpaceDim,zero);
      RealVect velocity;
      ppcharm.queryarr("velocity",velocitypp,0,SpaceDim);
      for (int i=0; i<SpaceDim; i++) velocity[i]=velocitypp[i];

      TestsIBC* tests = new TestsIBC();
      tests->setTestProblem(testString);
      tests->setFortranCommon(gamma,
                              ms,
                              center,
                              size,
                              velocity,
                              artificialViscosity,
                              riemannSolTol,
                              maxRsIter,
                              maxMach);
      //
      selfGravityPhysics->setPhysIBC(tests);

      delete tests;
    }
    else if (problemString == "ramp")
    {
      // Ramp problem
      Real alpha = 30.0;
      ppcharm.get("angle_deg",alpha);

      Real ms = ten;
      ppcharm.query("shock_mach",ms);

      Real xcorner = 0.1;
      ppcharm.get("xcorner",xcorner);

      if (verbosity >= 2)
      {
        pout() << "alpha = " << alpha << endl;
        pout() << "shock_mach = " << ms << endl;
        pout() << "xcorner = " << xcorner << endl;
      }
      Real smallPressure = 1.e-6;

      RampIBC* rampibc = new RampIBC();
      rampibc->setFortranCommon(smallPressure,
                                gamma,
                                alpha,
                                ms,
                                xcorner,
                                artificialViscosity,
                                riemannSolTol,
                                maxRsIter,
                                maxMach);

      //
      selfGravityPhysics->setPhysIBC(rampibc);

      delete rampibc;
    }
    else if (problemString == "dustcollapse")
    {
      DustCollapseIBC* dustCollapse = new DustCollapseIBC();

      Real cloudRadius= half;
      ppcharm.query("cloud_radius",cloudRadius);

      Real cloudDensity= one;
      ppcharm.query("cloud_density",cloudDensity);

      dustCollapse->setFortranCommon(gamma,
                                  cloudRadius,
                                  cloudDensity,
                                  artificialViscosity,
                                  riemannSolTol,
                                  maxRsIter,
                                  maxMach);
      //
      selfGravityPhysics->setPhysIBC(dustCollapse);

      delete dustCollapse;
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
    pout() << "\"charm.problem\" not specified in input file" << endl << endl;
    return;
  }


  // Print the parameters
  if ( verbosity >= 2 )
  {
    // Stop after this number of steps
    int nstop;
    ppcharm.get("max_step",nstop);
    // Stop when the simulation time get here
    Real stopTime;
    ppcharm.get("max_time",stopTime);

    pout() << "      " <<  endl;
    pout() << " problem = " << problemString << endl;
    pout() << "      " <<  endl;
    pout() << " max # steps      " << nstop    << endl;
    pout() << " max time         " << stopTime << endl;
    pout() << "      " <<  endl;
    pout() << " Riemann Solver: " <<  endl;
    pout() << "      " <<  endl;
    pout() << " gamma               = " << gamma << endl;
    pout() << " Tolerance           = " << riemannSolTol << endl;
    pout() << " Max Rs iter         = " << maxRsIter   << endl;
    pout() << " Mach n for Entr Sol = " << maxMach  << endl;
    pout() << "      " <<  endl;
    pout() << "      " <<  endl;
    pout() << " IC:  " <<  endl;
    pout() << "      " <<  endl;

    // display inputs
    if (problemString == "gastests")
    {
      Real ms,size;
      ppcharm.query("shock_mach",ms);
      ppcharm.query("size",size);
      vector<Real> center(SpaceDim), velocity(SpaceDim);
      ppcharm.queryarr("center",center,0,SpaceDim);
      ppcharm.queryarr("velocity",velocity,0,SpaceDim);

      pout() << " shock mach = " << ms << endl;
      pout() << " center = " << D_TERM(center[0] << "  " <<,
                                       center[1] << "  " <<,
                                       center[2] << ) endl;
      pout() << " size = " << size << endl;
      pout() << " velocity = " << D_TERM(velocity[0] << "  " <<,
                                         velocity[1] << "  " <<,
                                         velocity[2] << ) endl;
      pout() << "      " <<  endl;
    }
    else if (problemString == "dustcollapse")
    {
      Real cloudRadius, cloudDensity;
      ppcharm.query("cloud_radius",cloudRadius);
      ppcharm.query("cloud_density",cloudDensity);

      pout() << " cloudRadius = " << cloudRadius << endl;
      pout() << " cloudDensity= " << cloudDensity << endl;
      pout() << "      " <<  endl;
    }
    else if (problemString == "ramp")
    {
    }
    else
    {
      MayDay::Error("II-setupProblem \"problem\" for IBC is invalid");
    }
  }
}



void setupAMRFactory(AMRLevelSelfGravityFactory&  amrGodFactory,
                     SelfGravityPhysics*          selfGravityPhysics)
{
  // Read inputs that are prefixed with "charm."
  ParmParse ppcharm("charm");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppcharm.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppcharm.query("domain_length",domainLength);


  /* AMR */

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppcharm.query("tag_buffer_size",tagBufferSize);

  int maxInitRefLevel = 0;
  ppcharm.query ("max_init_ref_level",maxInitRefLevel);

  int gradientRefine = 0;
  ppcharm.query ("use_gradient_refine",gradientRefine);
  const bool useGradientRefine = (gradientRefine>0);
  Real gradRefineThreshold  = hundred;
  ppcharm.query ("grad_refine_thresh",gradRefineThreshold);

  vector<int> gradVarVector(2);
  ppcharm.queryarr("grad_var_interv",gradVarVector,0,2);
  Interval gradVarInterv(gradVarVector[0],gradVarVector[1]);

  int shockRefine = 0;
  ppcharm.query ("use_shock_refine",shockRefine);
  const bool useShockRefine = (shockRefine>0);
  Real presJumpThreshold = hundred;
  ppcharm.query ("pres_jump_thresh",presJumpThreshold);

  int overDenseRefine = 0;
  ppcharm.query ("use_over_dense_refine",overDenseRefine);
  const bool useOverDenseRefine = (overDenseRefine>0);
  Real cellMassThreshold = hundred;
  ppcharm.query ("cell_mass_thresh",cellMassThreshold);

  int jeansRefine = 0;
  ppcharm.query ("use_jeans_refine",jeansRefine);
  const bool useJeansRefine = (jeansRefine>0);
  Real jeansResolThreshold = four;
  ppcharm.query ("jeans_resol_thresh",jeansResolThreshold);

  Vector<RefineMode> refMode;
  vector<int> inRefMode;
  int irm = ppcharm.queryarr("refine_region_mode",inRefMode,0,maxInitRefLevel);
  if (irm==1)
  {
    refMode.resize(inRefMode.size());
    for (int i=0; i<inRefMode.size(); i++)
    {
      switch(inRefMode[i])
      {
        case 0: refMode[i] = FIX; break;
        case 1: refMode[i] = AND; break;
        case 2: refMode[i] = OR; break;
        default: MayDay::Error("setupAMRFactory:: Wrong input refine Mode");
      }
    }
  }

  Vector<Box> refRegion;
  bool useRegionRefine = false;
  for (int l=0; l<maxInitRefLevel; l++)
  {
    char boxLoPt[20]; sprintf(boxLoPt,"%s%d","ref_box_lo_",l);
    char boxHiPt[20]; sprintf(boxHiPt,"%s%d","ref_box_hi_",l);
    vector<int> loPt(SpaceDim,0);
    vector<int> hiPt(SpaceDim,0);
    int i = ppcharm.queryarr(boxLoPt,loPt,0,SpaceDim);
    i *= ppcharm.queryarr(boxHiPt,hiPt,0,SpaceDim);

    if (i==1)
    {
      if (l==0) useRegionRefine=true;
      IntVect lo(D_DECL(loPt[0], loPt[1], loPt[2]) );
      IntVect hi(D_DECL(hiPt[0], hiPt[1], hiPt[2]) );

      // remove that layer of cells which will be added as buffer zone.
      refRegion.push_back(grow(Box(lo,hi),-tagBufferSize));

      pout() << " lev " << l << " ref box " << Box(lo,hi) << endl;
      pout() << " lev " << l << " ref t box " << grow(Box(lo,hi),-tagBufferSize) << endl;
    }
  }

  /* Hydro */

  // Order of the normal predictor (PLM -> 1, PPM -> 2)
  std::string normalPred;
  int normalPredOrder;
  ppcharm.get("normal_predictor",normalPred);
  if (normalPred == "PLM" || normalPred == "plm")
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

  // Do 4th order (1) or 2nd order (0) slope computations
  int inFourthOrderSlopes = 1;
  ppcharm.query("use_fourth_order_slopes",inFourthOrderSlopes);
  const bool useFourthOrderSlopes = (inFourthOrderSlopes == 1);

  // Do slope limiting
  int inPrimLimiting = 1;
  ppcharm.get("use_prim_limiting",inPrimLimiting);
  const bool usePrimLimiting = (inPrimLimiting == 1);

  // Do slope limiting using characteristics
  int inCharLimiting = 0;
  ppcharm.get("use_char_limiting",inCharLimiting);
  const bool useCharLimiting = (inCharLimiting == 1);

  // Do slope flattening - only valid with 4th order slopes
  int inFlattening = 1;
  ppcharm.query("use_flattening",inFlattening);
  const bool useFlattening = (inFlattening == 1);

  CH_assert(useFourthOrderSlopes || !useFlattening);

  // Apply artificial viscosity
  int inUseArtificialViscosity = 1;
  ppcharm.query("use_artificial_viscosity",inUseArtificialViscosity);
  const bool useArtificialViscosity = (inUseArtificialViscosity == 1);

  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.1;
  ppcharm.query("artificial_viscosity",artificialViscosity);

  /* CFL's */

  // CFL multiplier
  Real cfl = 0.8;
  ppcharm.query("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppcharm.query("initial_cfl",initialCFL);

  /* Gravity */

  const StencilType arrStencil[3] =
  {
    TwoPts,FourPts,TenPts
  };
  int iStencil=0;
  ppcharm.query("force_stencil",iStencil);
  StencilType forceStencil = arrStencil[iStencil];

  int iuseDeltaPhiCorr=0;
  ppcharm.query("use_delta_phi_corr",iuseDeltaPhiCorr);
  const bool useDeltaPhiCorr = iuseDeltaPhiCorr>=1 ? true : false;


  RefCellTagger* refCellTagger = new RefCellTagger();
  refCellTagger->setRefineGradient(useGradientRefine,gradRefineThreshold,
                                   gradVarInterv);
  refCellTagger->setRefineShocks(useShockRefine,presJumpThreshold);
  refCellTagger->setRefineOverdense(useOverDenseRefine,cellMassThreshold);
  refCellTagger->setRefineJeans(useJeansRefine,jeansResolThreshold);
  refCellTagger->setRefineRegion(useRegionRefine,refRegion,refMode);

  /* Set up the AMRLevel... factory */

  amrGodFactory.define(cfl, domainLength, verbosity,
                    tagBufferSize,
                    maxInitRefLevel, initialCFL,
                    static_cast<GodunovPhysics*>(selfGravityPhysics),
                    normalPredOrder, useFourthOrderSlopes,
                    usePrimLimiting, useCharLimiting,
                    useFlattening ,useArtificialViscosity ,artificialViscosity,
                    refCellTagger, useDeltaPhiCorr, forceStencil);

  /* Print out the parameters */

  if ( verbosity >= 2 )
  {
    pout() << "      " <<  endl;
    pout() << " Parameters: " << endl;
    pout() << "      " <<  endl;
    pout() << " tagBufferSize = " << tagBufferSize << endl;
    pout() << " CFL Hydro     = " << cfl << endl;
    if (useArtificialViscosity)
      pout() << " Artificial Viscosity = " << artificialViscosity << endl;
    pout() << " initial CFL   = " << initialCFL << endl;
    pout() << " maxInitRefLevel  = " << maxInitRefLevel << endl;
  }

}


void setupAMRObject(int& a_maxStep, Real& a_maxTime, AMR& amr,
                    const AMRLevelSelfGravityFactory& amrGodFactory)
{
  // Read inputs that are prefixed with "charm."
  ParmParse ppcharm("charm");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppcharm.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  /* Problem domain */

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  ppcharm.getarr("num_cells",numCells,0,SpaceDim);

  CH_assert(D_TERM(   (numCells[0] > 0),
                && (numCells[1] > 0),
                && (numCells[2] > 0)));
  CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                && (numCells[1] % 2 == 0),
                && (numCells[2] % 2 == 0)));

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppcharm.queryarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
  {
    isPeriodic[dim] = (isPeriodica[dim] == 1);
    if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
      pout() << " Using Periodic BCs in direction: " << dim << endl;
  }


  /* AMR */

  // Maximum AMR level limit
  int maxLevel = 0;
  ppcharm.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppcharm.getarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppcharm.getarr("regrid_interval",regridIntervals,0,numReadLevels);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppcharm.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppcharm.query("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppcharm.query("fill_ratio",fillRatio);


  /* Timing and Output */

  // Stop after this number of steps
  ppcharm.get("max_step",a_maxStep);

  // Stop when the simulation time get here
  ppcharm.get("max_time",a_maxTime);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppcharm.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppcharm.query("plot_interval",plotInterval);

  /* dt control */

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppcharm.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppcharm.query("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppcharm.query("dt_tolerance_factor",dtToleranceFactor);

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrGodFactory);

  amr.verbosity(verbosity);

  //  amr.setQueueAgent(&queueAgent);

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
  if (ppcharm.contains("plot_prefix"))
  {
    std::string prefix;
    ppcharm.query("plot_prefix",prefix);
    amr.plotPrefix(prefix);
  }

  if (ppcharm.contains("chk_prefix"))
  {
    std::string prefix;
    ppcharm.query("chk_prefix",prefix);
    amr.checkpointPrefix(prefix);
  }

  // Read problem specified
  std::string problem;
  ppcharm.query("problem",problem);

  /*  print-out the parameters */

  if ( verbosity >= 2 )
  {
    pout() << "     " << endl;
    pout() << " Micellania:   " << endl;
    pout() << "     " << endl;
    pout() << " number of cells = " << D_TERM(numCells[0] << "  " <<,
                                              numCells[1] << "  " <<,
                                              numCells[2] << ) endl;

    pout() << " maximum level   = " << maxLevel << endl;

    pout() << " refinement ratio= ";
    for (int i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
    pout() << endl;

    pout() << " regrid interval = ";
    for (int i = 0; i < regridIntervals.size(); ++i) pout() <<
                        regridIntervals[i] << " ";
    pout() << endl;

    pout() << " blocking factor = " << blockFactor << endl;
    pout() << " max grid size   = " << maxGridSize << endl;
    pout() << " fill ratio      = " << fillRatio << endl;

    pout() << " checkpoint interval = " << checkpointInterval << endl;
    pout() << " plot interval       = " << plotInterval << endl;

    if (fixedDt > 0)
      pout() << " fixed dt           = " << fixedDt << endl;
    pout() << " maximum dt growth  = " << maxDtGrowth << endl;
    pout() << " dt tolerance factor= " << dtToleranceFactor << endl;
  }

  /* Set up input files */

  if (!ppcharm.contains("restart_file"))
  {
    if (!ppcharm.contains("fixed_hierarchy"))
    {
      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      amr.setupForNewAMRRun();
    }
    else
    {
      std::string gridFile;
      ppcharm.query("fixed_hierarchy",gridFile);

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
    ppcharm.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("setupAMRObject: restart only defined with hdf5");
#endif
  }
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




#ifdef TRAP_FPE
// Previous versions of glibc require the following code:
extern "C"
{
#include <fpu_control.h>
}
/* IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask */
static void __attribute__ ((constructor)) trapfpe(void)
{
  //pout() << " Turning on floating-point traps! " << endl;
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM);
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_OM | _FPU_MASK_UM);
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_UM);
  fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_DM | _FPU_MASK_UM);
  //fpu_control_t cw = _FPU_DEFAULT;
   _FPU_SETCW(cw);
   /* On x86, this expands to: */
   /* unsigned int cw = 0x037f & ~(0x01 | 0x04 | 0x08); */
   /* __asm__ ("fldcw %0" : : "m" (*&cw));              */
}
#endif
