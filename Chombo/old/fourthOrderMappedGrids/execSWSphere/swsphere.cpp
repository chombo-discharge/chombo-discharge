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

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif

#include "AMR.H"
#include "AMRLevelShallowWaterMappedFactory.H"
#include "MOLShallowWaterPhysics.H"
#include "LevelSWSourceTerm.H"

#include "MultiBlockCoordSys.H"
#include "CylinderEquiangularCS.H"
#include "CubedSphere2DCS.H"
#include "DoubleCartesianCS.H"
#include "TripleCartesianCS.H"

#include "PhysMappedIBC.H"
#include "BalancedFlowIBC.H"
#include "RossbyHaurwitzIBC.H"
#include "BarotropicInstabilityIBC.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "NamespaceHeader.H"

#ifdef USE_ARRAYVIEW

extern "C"
{
#include <fpu_control.h>
}
/* IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask
  ---(pm is kinda stupid)
*/
static void __attribute__ ((constructor)) trapfpe(void)
{
  fpu_control_t cw =
    _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
  _FPU_SETCW(cw);
}
#endif

// Coordinate system code
enum CoordSysCode
{
  CYLINDERSPOKES,
  CYLINDEREQUIANGULAR,
  CYLINDERTRANSITION,
  CUBEDSPHERE2D,
  CUBEDSPHERE2DFLAT,
  DOUBLECARTESIAN,
  TRIPLECARTESIAN
};

// amrShallowWater is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrShallowWater();

// One more function for MPI
void dumpmemoryatexit();

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
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);


  // Run amrShallowWater, i.e., do the computation
  amrShallowWater();

#ifdef CH_MPI

  MPI_Finalize();
#endif
}

///////////////////////////////////////////////////////////////////////////////

void amrShallowWater()
{
  // Read inputs that are prefixed with "swsphere."
  ParmParse ppswsphere("swsphere");

  // Determine the sample problem specified
  //XXX -- not used
  //XXXint problem = -1;
  std::string problemString;

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppswsphere.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Parameters specific to different sample problems
  int realDim = SpaceDim;

  // Stop after this number of steps
  int nstop = 0;
  ppswsphere.get("max_step", nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppswsphere.get("max_time", stopTime);

  // order of spatial  accuracy to use (2 or 4)
  int spaceOrder = 4;
  ppswsphere.query("spaceOrder", spaceOrder);

  // Initial values are average:  default TRUE
  int inInitialAverage = 1;
  ppswsphere.query("initial_average", inInitialAverage);
  bool initialAverage = (inInitialAverage == 1);

  // limit face values in advection:  default TRUE
  int inLimitFaceValues = 1;
  ppswsphere.query("limitFaceValues", inLimitFaceValues);
  bool limitFaceValues = (inLimitFaceValues == 1);

  // limit extrapolants in advection:  default TRUE
  //  int inLimitExtrapolants = 1;
  //  ppswsphere.query("limitExtrapolants", inLimitExtrapolants);
  //  bool limitExtrapolants = (inLimitExtrapolants == 1);

  // enforce a min value on advected quantity:  default FALSE
  int inEnforceMinVal = 0;
  ppswsphere.query("enforceMinVal", inEnforceMinVal);
  bool enforceMinVal = (inEnforceMinVal == 1);

  Real minVal = 0.0;
  if (enforceMinVal)
    {
      ppswsphere.get("minVal", minVal);
    }


  // forward Euler integration:  default FALSE
  int inForwardEuler = 0;
  ppswsphere.query("forward_euler", inForwardEuler);
  bool forwardEuler = (inForwardEuler == 1);

  // use hyperviscosity:  default FALSE
  int inUseHyperviscosity = 0;
  ppswsphere.query("useHyperviscosity", inUseHyperviscosity);
  bool useHyperviscosity = (inUseHyperviscosity == 1);

  // limit face values in advection
  Real hyperviscosity = -1.0;
  ppswsphere.query("hyperviscosity", hyperviscosity);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppswsphere.query("domain_length",domainLength);

  // Set the location of the lower left corner
  std::vector<Real> x0a(SpaceDim,0.0);
  ppswsphere.queryarr("x0",x0a,0,SpaceDim);
  RealVect x0;
  for ( int d=0 ; d<SpaceDim ; ++d ) x0[d] = x0a[d] ;

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  ppswsphere.queryarr("num_cells",numCells,0,SpaceDim);

 CH_assert(D_TERM(   (numCells[0] > 0),
                && (numCells[1] > 0),
                && (numCells[2] > 0)));
 CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                && (numCells[1] % 2 == 0),
                && (numCells[2] % 2 == 0)));

 // IGNORE dimensions other than the first one.
 int lengthCells = numCells[0];

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppswsphere.queryarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      isPeriodic[dim] = (isPeriodica[dim] == 1);
      if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
        pout() << "Using Periodic BCs in direction: " << dim << endl;
    }

  // Maximum AMR level limit
  int maxLevel = 0;
  ppswsphere.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppswsphere.queryarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppswsphere.queryarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppswsphere.query("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppswsphere.query ("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppswsphere.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppswsphere.query("max_grid_size",maxGridSize);
  int maxBaseGridSize = 0;
  ppswsphere.query("max_base_grid_size",maxBaseGridSize);

  Real fillRatio = 0.75;
  ppswsphere.query("fill_ratio",fillRatio);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppswsphere.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppswsphere.query("plot_interval",plotInterval);

  // Write mapped-grid geometry info:  default TRUE
  int inWriteMap = 1;
  ppswsphere.query("write_map", inWriteMap);
  bool writeMap = (inWriteMap == 1);

  // Write error:  default FALSE
  int inWriteError = 0;
  ppswsphere.query("write_error", inWriteError);
  bool writeError = (inWriteError == 1);

  // CFL multiplier
  Real cfl = 0.8;
  ppswsphere.query("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppswsphere.query("initial_cfl",initialCFL);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppswsphere.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppswsphere.query("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppswsphere.query("dt_tolerance_factor",dtToleranceFactor);

  // Print the parameters
  pout() << "maximum step = " << nstop << endl;
  pout() << "maximum time = " << stopTime << endl;

  pout() << "spatial order = " << spaceOrder << endl;
  pout() << "initial_average = "
         << (initialAverage ? "yes" : "no") << endl;

  pout() << "limit face values = " << limitFaceValues << endl;

  pout() << "enforceMinVal = " << enforceMinVal;
  if (enforceMinVal)
    {
      pout() << "; minVal = " << minVal;
    }
  pout() << endl;

  pout() << "use hyperviscosity = " << useHyperviscosity << endl;

  pout() << "hyperviscosity = " << hyperviscosity << endl;

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
  pout() << "max base grid size = " << maxBaseGridSize << endl ;
  pout() << "fill ratio = " << fillRatio << endl;

  pout() << "checkpoint interval = " << checkpointInterval << endl;
  pout() << "plot interval = " << plotInterval << endl;
  pout() << "CFL = " << cfl << endl;
  pout() << "initial CFL = " << initialCFL << endl;
  if (fixedDt > 0)
    {
      pout() << "fixed dt = " << fixedDt << endl;
    }
  pout() << "maximum dt growth = " << maxDtGrowth << endl;
  pout() << "dt tolerance factor = " << dtToleranceFactor << endl;

  // Create initial and boundary condition (IBC) object and initialize
  CoordSysCode csCode;
  RealVect centerPoint = RealVect::Zero;
  RealVect centralRectSize = RealVect::Unit;
  Real outerRadius = 1.5;

  vector<Real> centerPointVect(SpaceDim, 0.);
  if (ppswsphere.queryarr("center_point", centerPointVect, 0, SpaceDim))
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        centerPoint[idir] = centerPointVect[idir];
    }

  vector<Real> centralRectSizeVect(SpaceDim, 0.);
  if (ppswsphere.queryarr("central_rect_size", centralRectSizeVect, 0, SpaceDim))
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        centralRectSize[idir] = centralRectSizeVect[idir];
    }

  ppswsphere.query("outer_radius", outerRadius);

  ppswsphere.get("problem", problemString);

  PhysShallowWaterMappedIBC* ibc = NULL;

  // Steady-state Geostrophically Balanced Flow
  if (problemString == "balancedflow")
    {
      ParmParse balancedflowPP("balancedflow");

      Real bgVelocity = M_PI / 6.0;
      balancedflowPP.query("vel", bgVelocity);

      Real bgHeight = 2.94e4 / 9.80616 / 6.37122e6;
      balancedflowPP.query("height", bgHeight);

      BalancedFlowIBC* bfibc = new BalancedFlowIBC(bgVelocity, bgHeight);

      ibc = bfibc;
    }
  else if (problemString == "rossby")
    {
      ParmParse rossbyPP("rossby");

      Real a = 6.37122e6; // sphere radius, in meters
      Real day = 86164.1; // seconds in a sidereal day

      Real w = 7.848e-6 * day;
      rossbyPP.query("w", w);

      Real K = w;
      rossbyPP.query("K", K);

      Real h0 = 8e3 / a;
      rossbyPP.query("h0", h0);

      int R = 4;
      rossbyPP.query("R", R);

      RossbyHaurwitzIBC* rhibc = new RossbyHaurwitzIBC(w, K, h0, R);

      ibc = rhibc;
    }
  else if (problemString == "barotropic")
    {
      ParmParse barotropicPP("barotropic");

      Real a = 6.37122e6; // sphere radius, in meters
      Real day = 86400.; // 86164.1; // seconds in a sidereal day

      Real phi0 = M_PI / 7.;
      barotropicPP.query("phi0", phi0);

      Real phi1 = M_PI / 2. - phi0;
      barotropicPP.query("phi1", phi1);

      Real phi2 = M_PI / 4.;
      barotropicPP.query("phi2", phi2);

      Real umax = 80. * (day/a); // 80 m/s
      barotropicPP.query("umax", umax);

      Real hmean = 10000. / a; // 10 km
      barotropicPP.query("hmean", hmean);

      Real hhat = 120. / a; // 120 m
      barotropicPP.query("hhat", hhat);

      Real BIalpha = 1./3.;
      barotropicPP.query("alpha", BIalpha);

      Real BIbeta = 1./15.;
      barotropicPP.query("beta", BIbeta);

      BarotropicInstabilityIBC* biibc =
        new BarotropicInstabilityIBC(phi0, phi1, phi2, umax,
                                     hmean, hhat, BIalpha, BIbeta);
      ibc = biibc;
    }
  else
    {
      MayDay::Error("Invalid problemString");
    }

  // Constant of gravity
  Real gravity = 9.80616 * 86400.0 * 86400.0 / 6.37122e6;
  ppswsphere.query("gravity", gravity);

  // Rotation rate of the Earth
  Real omega = 7.292e-5 * 86400.0;
  ppswsphere.query("omega", omega);

  // Grid angle
  Real alpha = 0.0;
  ppswsphere.query("alpha", alpha);

  // Default is to use source term.
  bool useSourceTerm = true;
  {
    int useSourceTermInt = (useSourceTerm) ? 1 : 0;
    ppswsphere.query("use_source", useSourceTermInt);
    useSourceTerm = (useSourceTermInt == 1);
  }

  LevelSourceTerm* sourceTerm = NULL;

  // Set constants
  ibc->setFortranCommon(gravity, omega, alpha);

  // // Source term multiplier
  //  Real sourceTermScaling = 0.0;

  // Problem domain
  // ProblemDomain probDomain(IntVect::Zero,
  //                         IntVect(numCells, numCells));

  // Create physics object
  MOLShallowWaterPhysics* swPhysics =
    new MOLShallowWaterPhysics();

  swPhysics->setPhysIBC(ibc);

  if (useSourceTerm)
    {
      LevelSWSourceTerm* swSourceTerm = new LevelSWSourceTerm();
      // Cast to source term base class pointer for technical reasons
      sourceTerm = dynamic_cast<LevelSourceTerm*>(swSourceTerm);
    }

  // Create coordinate system factory and initialize
  ParmParse coordSysPP("coordsys");
  string coordSysString("cartesian");
  coordSysPP.get("type", coordSysString);
  if (coordSysString == "cylinderspokes")
    csCode = CYLINDERSPOKES;
  else if (coordSysString == "cylinderequiangular")
    csCode = CYLINDEREQUIANGULAR;
  else if (coordSysString == "cylindertransition")
    csCode = CYLINDERTRANSITION;
  else if (coordSysString == "cubedsphere2d")
    csCode = CUBEDSPHERE2D;
  else if (coordSysString == "cubedsphere2dflat")
    csCode = CUBEDSPHERE2DFLAT;
  else if (coordSysString == "doublecartesian")
    csCode = DOUBLECARTESIAN;
  else if (coordSysString == "triplecartesian")
    csCode = TRIPLECARTESIAN;

  // Set up the coordinate system... factory
  //
  Real dx = domainLength / (Real(numCells[0]));
  RealVect dxVect = dx * RealVect::Unit;
  Vector<int> ratios(refRatios);
  MultiBlockCoordSysFactory* coordSysFactPtr;
  IntVect levelDomainLo, levelDomainHi;
  Real fullDomainLength = domainLength;
  switch (csCode)
    {
    case CYLINDERSPOKES:
      {
        // smdPtr = new CylindricalSpokesDomain(domain, RealVect::Zero,
        //                                      bxWidth,
        //                                      outerRadius);
        //           ghostFactor = radius + 2;
        MayDay::Error("CylinderSpokes not implemented yet.");
        break;
      }
    case CYLINDEREQUIANGULAR:
      {
        CylinderEquiangularCSFactory* cylinderCSFactPtr =
          new CylinderEquiangularCSFactory;
        cylinderCSFactPtr->setCenterPoint(centerPoint);
        cylinderCSFactPtr->setCentralRectSize(centralRectSize);
        cylinderCSFactPtr->setOuterRadius(outerRadius);
        coordSysFactPtr = cylinderCSFactPtr;
        levelDomainLo = IntVect(D_DECL(-2*lengthCells,
                                       -2*lengthCells,
                                       0));
        levelDomainHi = IntVect(D_DECL(3*lengthCells-1,
                                       3*lengthCells-1,
                                       lengthCells-1));
        fullDomainLength = 5. * domainLength;
        break;
      }
    case CYLINDERTRANSITION:
      {
        // smdPtr = new CylindricalTransitionDomain(domain, RealVect::Zero,
        //                                          bxWidth, outerRadius);
        // ghostFactor = radius + 1;
        MayDay::Error("CylinderTransition not implemented yet.");
        break;
      }
    case CUBEDSPHERE2D:
      {
        // smdPtr = new SphericalDomain(domain, RealVect::Zero,
        //                              bxWidth, outerRadius);
        // ghostFactor = radius + 2;
        CubedSphere2DCSFactory* cubedSphere2DCSFactPtr =
          new CubedSphere2DCSFactory;
        // cubedSphere2DCSFactPtr->setCenterPoint(centerPoint);
        // cubedSphere2DCSFactPtr->setCentralRectSize(centralRectSize);
        // cubedSphere2DCSFactPtr->setOuterRadius(outerRadius);
        coordSysFactPtr = cubedSphere2DCSFactPtr;
        levelDomainLo = IntVect::Zero;
        levelDomainHi = IntVect(D_DECL(11*lengthCells-1,
                                       lengthCells-1,
                                       0)); // third dimension unused
        fullDomainLength = 11. * M_PI / 2. * domainLength;
        realDim = 3; // but SpaceDim == 2
        break;
      }
    case CUBEDSPHERE2DFLAT:
      {
        CubedSphere2DCSFactory* cubedSphere2DCSFactPtr =
          new CubedSphere2DCSFactory;
        coordSysFactPtr = cubedSphere2DCSFactPtr;
        ((CubedSphere2DCSFactory*) coordSysFactPtr)->setFlatMap(true);
        levelDomainLo = IntVect::Zero;
        levelDomainHi = IntVect(D_DECL(11*lengthCells-1,
                                       lengthCells-1,
                                       0)); // third dimension unused
        fullDomainLength = 11. * M_PI / 2. * domainLength;
        break;
      }
    case DOUBLECARTESIAN:
      {
        DoubleCartesianCSFactory* doubleCartesianCSFactPtr =
          new DoubleCartesianCSFactory;
        coordSysFactPtr = doubleCartesianCSFactPtr;
        levelDomainLo = IntVect::Zero;
        levelDomainHi = (3*lengthCells-1) * IntVect::Unit;
        fullDomainLength = 3. * domainLength;
        break;
      }
    case TRIPLECARTESIAN:
      {
        TripleCartesianCSFactory* tripleCartesianCSFactPtr =
          new TripleCartesianCSFactory;
        coordSysFactPtr = tripleCartesianCSFactPtr;
        levelDomainLo = IntVect::Zero;
        levelDomainHi = (5*lengthCells-1) * IntVect::Unit;
        fullDomainLength = 5. * domainLength;
        break;
      }
    }


  Box levelDomainBox(levelDomainLo, levelDomainHi);
  ProblemDomain levelDomain(levelDomainBox);
  MultiBlockCoordSys* coordSysPtr =
    coordSysFactPtr->getCoordSys(levelDomain, dxVect);
  const Vector<Box>& blockBoxes = coordSysPtr->mappingBlocks();
  int nblocks = blockBoxes.size();
  Vector<Box> allBoxes;
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      Vector<Box> thisBlockBoxes;
      domainSplit(blockBoxes[iblock], thisBlockBoxes,
                  maxGridSize, blockFactor);
      allBoxes.append(thisBlockBoxes);
    }
  //  Vector<int> allProcs(allBoxes.size());
  //  LoadBalance(allProcs, allBoxes);
  //  DisjointBoxLayout grids(allBoxes, allProcs);

  // { // scoping trick

  // Set up the AMRLevel... factory
  //
  AMRLevelShallowWaterMappedFactory amrShallowWaterFact;
  amrShallowWaterFact.spaceOrder(spaceOrder);
  amrShallowWaterFact.initialAverage(initialAverage);
  amrShallowWaterFact.limitFaceValues(limitFaceValues);
  amrShallowWaterFact.enforceMinVal(enforceMinVal, minVal);
  // amrShallowWaterFact.useHyperviscosity(useHyperviscosity);
  //  if (hyperviscosity>=0)
  //     amrShallowWaterFact.hyperviscosity(hyperviscosity);
  amrShallowWaterFact.CFL(cfl);
  amrShallowWaterFact.domainLength(fullDomainLength); // WAS (domainLength);
  amrShallowWaterFact.forwardEuler(forwardEuler);
  amrShallowWaterFact.refinementThreshold(refineThresh);
  amrShallowWaterFact.tagBufferSize(tagBufferSize);
  amrShallowWaterFact.verbosity(verbosity);
  amrShallowWaterFact.initialDtMultiplier(initialCFL);
  amrShallowWaterFact.molPhysics(swPhysics);
  // amrShallowWaterFact.IBC(ibc);
  amrShallowWaterFact.coordinateSystemFactory(coordSysFactPtr);
  amrShallowWaterFact.useSourceTerm(useSourceTerm);
  if (useSourceTerm)
    {
      amrShallowWaterFact.sourceTerm(sourceTerm);
    }

  // Set up output files -- do this with AMRLevelFactory
  // in order to be able to do mapped-grid output as well.
  std::string prefix;
  if (ppswsphere.contains("plot_prefix"))
    {
      ppswsphere.query("plot_prefix",prefix);
      amrShallowWaterFact.plotPrefix(prefix);
    }
  amrShallowWaterFact.writeMap(writeMap);
  amrShallowWaterFact.writeError(writeError);

  { // scope of AMR amr;
  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel, refRatios, levelDomain, &amrShallowWaterFact);

  if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }

  if (ppswsphere.contains("plot_prefix"))
    {
      ppswsphere.query("plot_prefix",prefix);
      amr.plotPrefix(prefix);
      amrShallowWaterFact.plotPrefix(prefix);
    }

  if (ppswsphere.contains("chk_prefix"))
    {
      std::string chkprefix;
      ppswsphere.query("chk_prefix",chkprefix);
      amr.checkpointPrefix(chkprefix);
    }

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  if ( maxBaseGridSize != 0)
  {
    //defaults to maxGridSize
    amr.maxBaseGridSize(maxBaseGridSize);
  }
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


  amr.verbosity(verbosity);

  // Set up input files
  if (!ppswsphere.contains("restart_file"))
    {
      if (!ppswsphere.contains("fixed_hierarchy"))
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
          amr.setupForNewAMRRun();
        }
      else
        {
          std::string gridFile;
          ppswsphere.query("fixed_hierarchy",gridFile);

          int numLevels = maxLevel+1;
          // initialize from a list of grids in "gridFile"
          Vector<Vector<Box> > amrGrids(numLevels);
          amrGrids[0] = allBoxes;
          //          setupFixedGrids(amrGrids,
          //                          probDomain,
          //                          maxLevel,
          //                          maxGridSize,
          //                          blockFactor,
          //                          verbosity,
          //                          gridFile);
          //          readBoxes(amrGrids, ppswsphere, probDomain,
          //                    maxGridSize, blockFactor,
          //                    numLevels, refRatios, verbosity);
          amr.setupForFixedHierarchyRun(amrGrids,1);
        }
    }
  else
    {
      std::string restartFile;
      ppswsphere.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
      HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      amr.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrAdvect restart only defined with hdf5");
#endif
    }

  amr.run(stopTime,nstop);

  amr.conclude();

  } // end scope of amr

  // clean up memory
  // as things stand now, coordSysFact is deleted by amrShallowWaterFactory's destructor
  delete ibc;
  delete swPhysics;
  delete coordSysPtr;

  // delete coordSysFactPtr;

}

/*
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
*/

#include "NamespaceFooter.H"

