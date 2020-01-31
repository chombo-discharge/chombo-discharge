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
#include <string>

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif

#include "AMR.H"
#include "AMRLevelAdvectMappedFactory.H"
#include "MOLAdvectionPhysics.H"
#include "CircleAdvectMultiMappedIBC.H"
#include "GaussianAdvectMultiMappedIBC.H"
#include "PhysAdvectMappedIBC.H"

#include "MultiBlockCoordSys.H"
#include "CylinderEquiangularCS.H"
#include "CubedSphere2DCS.H"
#include "DoubleCartesianCS.H"
#include "TripleCartesianCS.H"

#include "AdvectCosineBellIBC.H"
#include "AdvectConstantIBC.H"
#include "DeformationalFlowIBC.H"
#include "GaussianHillsFlowIBC.H"
#include "CosineBellsFlowIBC.H"
#include "SlottedCylindersFlowIBC.H"
#include "CorrelatedCosineBellsFlowIBC.H"
#include "CosineBellsDivergentFlowIBC.H"

// #include "generalFuncs.H"

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

void amrAdvect();

// amrAdvect is a function (as opposed to inline in main()) to get
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

// One more function for MPI
void dumpmemoryatexit();

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
  setChomboMPIErrorHandler();
#endif

  // Check for an input file
  char* inFile = NULL;
#ifdef USE_ARRAYVIEW
  trapfpe();
#endif

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrAdvect...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  // Run amrAdvect, i.e., do the computation
  amrAdvect();

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
#endif

  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void amrAdvect()
{

  // Read inputs that are prefixed with "advect."
  ParmParse ppadvect("advect");

  // Determine the sample problem specified
  //XXX -- not used
  //XXXint problem = -1;
  std::string problemString;


  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppadvect.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Parameters specific to different sample problems
  int realDim = SpaceDim;

  // Stop after this number of steps
  int nstop = 0;
  ppadvect.query("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppadvect.query("max_time",stopTime);

  // order of spatial  accuracy to use (2 or 4)
  int spaceOrder = 4;
  ppadvect.query("spaceOrder", spaceOrder);

  // Initial values are average:  default TRUE
  int inInitialAverage = 1;
  ppadvect.query("initial_average", inInitialAverage);
  bool initialAverage = (inInitialAverage == 1);

  // limit face values in advection:  default TRUE
  int inLimitFaceValues = 1;
  ppadvect.query("limitFaceValues", inLimitFaceValues);
  bool limitFaceValues = (inLimitFaceValues == 1);

  // limit extrapolants in advection:  default TRUE
  //  int inLimitExtrapolants = 1;
  //  ppadvect.query("limitExtrapolants", inLimitExtrapolants);
  //  bool limitExtrapolants = (inLimitExtrapolants == 1);

  // enforce a min value on advected quantity:  default FALSE
  int inEnforceMinVal = 0;
  ppadvect.query("enforceMinVal", inEnforceMinVal);
  bool enforceMinVal = (inEnforceMinVal == 1);

  Real minVal = 0.0;
  if (enforceMinVal)
    {
      ppadvect.get("minVal", minVal);
    }


  // forward Euler integration:  default FALSE
  int inForwardEuler = 0;
  ppadvect.query("forward_euler", inForwardEuler);
  bool forwardEuler = (inForwardEuler == 1);

  // use hyperviscosity:  default FALSE
  int inUseHyperviscosity = 0;
  ppadvect.query("useHyperviscosity", inUseHyperviscosity);
  bool useHyperviscosity = (inUseHyperviscosity == 1);

  // limit face values in advection
  Real hyperviscosity = -1.0;
  ppadvect.query("hyperviscosity", hyperviscosity);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppadvect.query("domain_length",domainLength);

  // Set the location of the lower left corner
  std::vector<Real> x0a(SpaceDim,0.0);
  ppadvect.queryarr("x0",x0a,0,SpaceDim);
  RealVect x0;
  for ( int d=0 ; d<SpaceDim ; ++d ) x0[d] = x0a[d] ;

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  ppadvect.queryarr("num_cells",numCells,0,SpaceDim);
  // Check that every component of numCells is positive
  for (int i = 0; i < SpaceDim; ++i)
    {
      CH_assert(numCells[i] > 0);
    }
  // IGNORE dimensions other than the first one.
  int lengthCells = numCells[0];

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppadvect.queryarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      isPeriodic[dim] = (isPeriodica[dim] == 1);
      if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
        pout() << "Using Periodic BCs in direction: " << dim << endl;
    }

  // Maximum AMR level limit
  int maxLevel = 0;
  ppadvect.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppadvect.queryarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppadvect.queryarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppadvect.query("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppadvect.query ("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppadvect.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppadvect.query("max_grid_size",maxGridSize);
  int maxBaseGridSize = 0;
  ppadvect.query("max_base_grid_size",maxBaseGridSize);

  Real fillRatio = 0.75;
  ppadvect.query("fill_ratio",fillRatio);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppadvect.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppadvect.query("plot_interval",plotInterval);

  // Write mapped-grid geometry info:  default TRUE
  int inWriteMap = 1;
  ppadvect.query("write_map", inWriteMap);
  bool writeMap = (inWriteMap == 1);

  // Write error:  default FALSE
  int inWriteError = 0;
  ppadvect.query("write_error", inWriteError);
  bool writeError = (inWriteError == 1);

  // CFL multiplier
  Real cfl = 0.8;
  ppadvect.query("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppadvect.query("initial_cfl",initialCFL);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppadvect.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppadvect.query("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppadvect.query("dt_tolerance_factor",dtToleranceFactor);

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

  pout() << "number of cells =";
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      pout() << " " << numCells[idir];
    }
  pout() << endl;

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
  if (ppadvect.queryarr("center_point", centerPointVect, 0, SpaceDim))
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        centerPoint[idir] = centerPointVect[idir];
    }

  vector<Real> centralRectSizeVect(SpaceDim, 0.);
  if (ppadvect.queryarr("central_rect_size", centralRectSizeVect, 0, SpaceDim))
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        centralRectSize[idir] = centralRectSizeVect[idir];
    }

  ppadvect.query("outer_radius", outerRadius);


  ppadvect.get("problem", problemString);

  PhysAdvectMappedIBC* ibc = NULL;

  if (problemString == "cosinebelladvection")
    {
      ParmParse cosinebellPP("cosinebell");
      Real ambientDensity = 0.0;
      cosinebellPP.query("ambient_density", ambientDensity);
      Real deltaDensity = 1000.0;
      cosinebellPP.query("delta_density", deltaDensity);
      Real size = 0.3;
      cosinebellPP.query("size", size);
      Real advectVelocity = M_PI / 6.0;
      cosinebellPP.query("vel", advectVelocity);
      Real advectAngle = 0.0; //M_PI / 2.0;
      cosinebellPP.query("angle", advectAngle);
      Real evalTime = 0.0;
      cosinebellPP.query("time", evalTime);

      AdvectCosineBellIBC* advectibc = new AdvectCosineBellIBC(ambientDensity,
                                                               deltaDensity,
                                                               centerPoint,
                                                               size,
                                                               advectVelocity,
                                                               advectAngle,
                                                               evalTime);

      ibc = advectibc;
    }
  else if (problemString == "constant")
    {
      ParmParse constantPP("constant");
      Real fieldVal = 1.0;
      constantPP.query("field", fieldVal);
      Real advectVelocity = M_PI / 6.0;
      constantPP.query("vel", advectVelocity);
      Real advectAngle = 0.0; //M_PI / 2.0;
      constantPP.query("angle", advectAngle);
      AdvectConstantIBC* advectibc = new AdvectConstantIBC(fieldVal,
                                                           advectVelocity,
                                                           advectAngle);
      ibc = advectibc;
    }
  else if (problemString == "deformationalflow")
    {
      ParmParse deformationalflowPP("deformationalflow");
      Real ambientDensity = 0.1;
      deformationalflowPP.query("ambient_density", ambientDensity);
      Real deltaDensity = 0.9;
      deformationalflowPP.query("delta_density", deltaDensity);
      Real size = 0.5;
      deformationalflowPP.query("size", size);
      // this is used ONLY in initialization:
      // const double InitialTheta = m_evalTime * m_advectVelocity;
      Real advectVelocity = M_PI / 6.0;
      deformationalflowPP.query("vel", advectVelocity);
      Real advectAngle = 0.0; //M_PI / 2.0;
      deformationalflowPP.query("angle", advectAngle);
      Real evalTime = 0.0;
      deformationalflowPP.query("time", evalTime);

      DeformationalFlowIBC* dfibc = new DeformationalFlowIBC(ambientDensity,
                                                             deltaDensity,
                                                             centerPoint,
                                                             size,
                                                             advectVelocity,
                                                             advectAngle,
                                                             evalTime);

      ibc = dfibc;
    }
  else if (problemString == "cosinebellsdivergentflow")
    {
     ParmParse cosinebellsdivergentflowPP("cosinebellsdivergentflow");
      Real hmax = 1.;
      cosinebellsdivergentflowPP.query("hmax", hmax);
      Real radius = 0.5;
      cosinebellsdivergentflowPP.query("radius", radius);
      int bells = 2;
      cosinebellsdivergentflowPP.query("bells", bells);

      Vector<Real> longitude(bells);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      cosinebellsdivergentflowPP.queryarr("longitude", longitude, 0, bells);

      Vector<Real> latitude(bells, 0.);
      cosinebellsdivergentflowPP.queryarr("latitude", latitude, 0, bells);

      Real background = 0.1;
      cosinebellsdivergentflowPP.query("background", background);
      Real delta = 0.9;
      cosinebellsdivergentflowPP.query("delta", delta);

      Real period = 5.;
      cosinebellsdivergentflowPP.query("period", period);
      Real k = 1.;
      cosinebellsdivergentflowPP.query("k", k);
      Real evalTime = 0.0;
      cosinebellsdivergentflowPP.query("time", evalTime);

      CosineBellsDivergentFlowIBC* dfibc = new CosineBellsDivergentFlowIBC(
                                                             hmax,
                                                         radius,
                                                         longitude,
                                                         latitude,
                                                         background,
                                                         delta,
                                                         period,
                                                         k,
                                                         evalTime);
      ibc = dfibc;
    }
  else if (problemString == "cosinebellsflow")
    {
      ParmParse cosinebellsflowPP("cosinebellsflow");
      Real hmax = 1.;
      cosinebellsflowPP.query("hmax", hmax);
      Real radius = 0.5;
      cosinebellsflowPP.query("radius", radius);
      int bells = 2;
      cosinebellsflowPP.query("bells", bells);

      Vector<Real> longitude(bells);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      cosinebellsflowPP.queryarr("longitude", longitude, 0, bells);

      Vector<Real> latitude(bells, 0.);
      cosinebellsflowPP.queryarr("latitude", latitude, 0, bells);

      Real background = 0.1;
      cosinebellsflowPP.query("background", background);
      Real delta = 0.9;
      cosinebellsflowPP.query("delta", delta);

      Real period = 5.;
      cosinebellsflowPP.query("period", period);
      Real kappa = 2.;
      cosinebellsflowPP.query("kappa", kappa);
      Real evalTime = 0.0;
      cosinebellsflowPP.query("time", evalTime);

      CosineBellsFlowIBC* dfibc = new CosineBellsFlowIBC(hmax,
                                                         radius,
                                                         longitude,
                                                         latitude,
                                                         background,
                                                         delta,
                                                         period,
                                                         kappa,
                                                         evalTime);
      ibc = dfibc;
    }
  else if (problemString == "correlatedcosinebellsflow")
    {
      ParmParse cosinebellsflowPP("correlatedcosinebellsflow");

      int ncoeffs = 3;
      cosinebellsflowPP.query("ncoeffs", ncoeffs);
      // coefficients, starting with constant term
      Vector<Real> coeffs(ncoeffs);
      coeffs[0] = 0.9;
      coeffs[1] = 0.;
      coeffs[2] = -0.8;
      cosinebellsflowPP.queryarr("coeffs", coeffs, 0, ncoeffs);

      Real hmax = 1.;
      cosinebellsflowPP.query("hmax", hmax);
      Real radius = 0.5;
      cosinebellsflowPP.query("radius", radius);
      int bells = 2;
      cosinebellsflowPP.query("bells", bells);

      Vector<Real> longitude(bells);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      cosinebellsflowPP.queryarr("longitude", longitude, 0, bells);

      Vector<Real> latitude(bells, 0.);
      cosinebellsflowPP.queryarr("latitude", latitude, 0, bells);

      Real background = 0.1;
      cosinebellsflowPP.query("background", background);
      Real delta = 0.9;
      cosinebellsflowPP.query("delta", delta);

      Real period = 5.;
      cosinebellsflowPP.query("period", period);
      Real kappa = 2.;
      cosinebellsflowPP.query("kappa", kappa);
      Real evalTime = 0.0;
      cosinebellsflowPP.query("time", evalTime);

      CorrelatedCosineBellsFlowIBC* dfibc =
        new CorrelatedCosineBellsFlowIBC(coeffs,
                                         hmax,
                                         radius,
                                         longitude,
                                         latitude,
                                         background,
                                         delta,
                                         period,
                                         kappa,
                                         evalTime);
      ibc = dfibc;
    }
  else if (problemString == "gaussianhillsflow")
    {
      ParmParse gaussianhillsflowPP("gaussianhillsflow");
      Real hmax = 0.95;
      gaussianhillsflowPP.query("hmax", hmax);
      Real width = 5.;
      gaussianhillsflowPP.query("width", width);
      int hills = 2;
      gaussianhillsflowPP.query("hills", hills);

      Vector<Real> longitude(hills);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      gaussianhillsflowPP.queryarr("longitude", longitude, 0, hills);

      Vector<Real> latitude(hills, 0.);
      gaussianhillsflowPP.queryarr("latitude", latitude, 0, hills);

      Real period = 5.;
      gaussianhillsflowPP.query("period", period);
      Real kappa = 2.;
      gaussianhillsflowPP.query("kappa", kappa);
      Real evalTime = 0.0;
      gaussianhillsflowPP.query("time", evalTime);

      GaussianHillsFlowIBC* dfibc = new GaussianHillsFlowIBC(hmax,
                                                             width,
                                                             longitude,
                                                             latitude,
                                                             period,
                                                             kappa,
                                                             evalTime);
      ibc = dfibc;
    }
  else if (problemString == "slottedcylindersflow")
    {
      ParmParse slottedcylindersflowPP("slottedcylindersflow");

      Real hmax = 1.;
      slottedcylindersflowPP.query("hmax", hmax);
      Real radius = 0.5;
      slottedcylindersflowPP.query("radius", radius);
      int cylinders = 2;
      slottedcylindersflowPP.query("cylinders", cylinders);

      Vector<Real> longitude(cylinders);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      slottedcylindersflowPP.queryarr("longitude", longitude, 0, cylinders);

      Vector<Real> latitude(cylinders, 0.);
      slottedcylindersflowPP.queryarr("latitude", latitude, 0, cylinders);

      Real background = 0.1;
      slottedcylindersflowPP.query("background", background);

      Real period = 5.;
      slottedcylindersflowPP.query("period", period);
      Real kappa = 2.;
      slottedcylindersflowPP.query("kappa", kappa);
      Real evalTime = 0.0;
      slottedcylindersflowPP.query("time", evalTime);

      SlottedCylindersFlowIBC* dfibc =
        new SlottedCylindersFlowIBC(
                                         hmax,
                                         radius,
                                         longitude,
                                         latitude,
                                         background,
                                         period,
                                         kappa,
                                         evalTime);
      ibc = dfibc;
    }
  else if (problemString == "circle")
    {
      CircleAdvectMultiMappedIBC* advectibc = new CircleAdvectMultiMappedIBC;
      ParmParse circlePP("circle");
      std::string velType;
      circlePP.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          circlePP.getarr("vel", vela,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              vel[idir] = vela[idir];
            }
          advectibc->setUniformVel(vel);
        }
      else if (velType == "solidBody")
        {
          Real omega;
          circlePP.get("omega", omega);
          RealVect rotationCenter;
          Vector<Real> rotCtra(SpaceDim);
          circlePP.getarr("rotationCenter", rotCtra,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              rotationCenter[idir] = rotCtra[idir];
            }
          advectibc->setSolidBodyRotation(rotationCenter, omega);
        }
      else
        {
          MayDay::Error("CircleAdvectMultiMappedIBC -- bad velType");
        }

      Real r0;
      circlePP.query("radius",r0);

      // Set the location of the center of the circle
      std::vector<Real> centera(SpaceDim,0.0);
      circlePP.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;

      advectibc->setParams(r0,center,x0);

      ibc  = advectibc;
    }
  else if (problemString == "gaussian")
    {
      GaussianAdvectMultiMappedIBC* advectibc = new GaussianAdvectMultiMappedIBC;
      ParmParse gaussianPP("gaussian");
      std::string velType;
      gaussianPP.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          gaussianPP.getarr("vel", vela,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              vel[idir] = vela[idir];
            }
          advectibc->setUniformVel(vel);
        }
      else if (velType == "solidBody")
        {
          Real omega;
          gaussianPP.get("omega", omega);
          RealVect rotationCenter;
          Vector<Real> rotCtra(SpaceDim);
          gaussianPP.getarr("rotationCenter", rotCtra,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              rotationCenter[idir] = rotCtra[idir];
            }
          advectibc->setSolidBodyRotation(rotationCenter, omega);
        }
      else if (velType == "transosc")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          gaussianPP.getarr("vel", vela,0,SpaceDim);
          if (gaussianPP.contains("transvel"))
            {
              gaussianPP.getarr("transvel", vela,0,SpaceDim);
            }
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              vel[idir] = vela[idir];
            }
          Real oscAmp;
          gaussianPP.get("oscAmp", oscAmp);
          advectibc->setTranslatingOscillation(vel,oscAmp);
        }
      else
        {
          MayDay::Error("GaussianAdvectMappedIBC -- bad velType");
        }

      Real r0;
      gaussianPP.query("radius",r0);

      // Set the location of the center of the gaussian
      std::vector<Real> centera(SpaceDim,0.0);
      gaussianPP.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;

      Real evalTime = 0.0;
      gaussianPP.query("time", evalTime);

      advectibc->setParams(r0,center,x0,evalTime);

      ibc  = advectibc;
    }
  else
    {
      MayDay::Error("Invalid problemString");
    }
  /*
  else if (problemString == "circle")
    {
      CircleAdvectMappedIBC* advectibc = new CircleAdvectMappedIBC;
      ParmParse circlePP("circle");
      std::string velType;
      circlePP.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          circlePP.getarr("vel", vela,0,SpaceDim);
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)

          advectibc->setUniformVel(vel);
        }
      else if (velType == "solidBody")
        {
          Real omega;
          circlePP.get("omega", omega);
          RealVect rotationCenter;
          Vector<Real> rotCtra(SpaceDim);
          circlePP.getarr("rotationCenter", rotCtra,0,SpaceDim);
          D_TERM(rotationCenter[0]=rotCtra[0];,
                 rotationCenter[1]=rotCtra[1];,
                 rotationCenter[2]=rotCtra[2];)

         advectibc->setSolidBodyRotation(rotationCenter, omega);
        }
      else
        {
          MayDay::Error("CircleAdvectMappedIBC -- bad velType");
        }

      Real r0;
      circlePP.query("radius",r0);

      // Set the location of the center of the circle
      std::vector<Real> centera(SpaceDim,0.0);
      circlePP.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;



      advectibc->setParams(r0,center,x0);

      ibc  = advectibc;
    }
  */

  MOLAdvectionPhysics* advectionPhysics = new MOLAdvectionPhysics();
  advectionPhysics->setPhysIBC(ibc);
  // godunovPhysics = static_cast<MOLPhysics*> (advectionPhysics);

  //  ProblemDomain probDomain (IntVect::Zero,
  //                            IntVect(D_DECL(numCells[0]-1,
  //                                           numCells[1]-1,
  //                                           numCells[2]-1)),
  //                            isPeriodic);

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
        levelDomainLo = IntVect(D_DECL6(-2*lengthCells,
                                        -2*lengthCells,
                                        0, 0, 0, 0));
        levelDomainHi = IntVect(D_DECL6(3*lengthCells-1,
                                        3*lengthCells-1,
                                        lengthCells-1,
                                        lengthCells-1,
                                        lengthCells-1,
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
        levelDomainHi = IntVect(D_DECL6(11*lengthCells-1,
                                        lengthCells-1,
                                        0, 0, 0, 0)); // dimensions > 2 unused
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
        levelDomainHi = IntVect(D_DECL6(11*lengthCells-1,
                                        lengthCells-1,
                                        0, 0, 0, 0)); // dimensions > 2 unused
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
  AMRLevelAdvectMappedFactory amrAdvectFact;
  amrAdvectFact.spaceOrder(spaceOrder);
  amrAdvectFact.initialAverage(initialAverage);
  amrAdvectFact.limitFaceValues(limitFaceValues);
  amrAdvectFact.enforceMinVal(enforceMinVal, minVal);
  // amrAdvectFact.useHyperviscosity(useHyperviscosity);
  //  if (hyperviscosity>=0)
  //     amrAdvectFact.hyperviscosity(hyperviscosity);
  amrAdvectFact.CFL(cfl);
  amrAdvectFact.domainLength(fullDomainLength); // WAS (domainLength);
  amrAdvectFact.forwardEuler(forwardEuler);
  amrAdvectFact.refinementThreshold(refineThresh);
  amrAdvectFact.tagBufferSize(tagBufferSize);
  amrAdvectFact.verbosity(verbosity);
  amrAdvectFact.initialDtMultiplier(initialCFL);
  amrAdvectFact.molPhysics(advectionPhysics);
  // amrAdvectFact.IBC(ibc);
  amrAdvectFact.coordinateSystemFactory(coordSysFactPtr);

  // Set up output files -- do this with AMRLevelFactory
  // in order to be able to do mapped-grid output as well.
  std::string prefix;
  if (ppadvect.contains("plot_prefix"))
    {
      ppadvect.query("plot_prefix",prefix);
      amrAdvectFact.plotPrefix(prefix);
    }
  amrAdvectFact.writeMap(writeMap);
  amrAdvectFact.writeError(writeError);

  { // scope of AMR amr;
  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel, refRatios, levelDomain, &amrAdvectFact);

  if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }

  if (ppadvect.contains("plot_prefix"))
    {
      ppadvect.query("plot_prefix",prefix);
      amr.plotPrefix(prefix);
      amrAdvectFact.plotPrefix(prefix);
    }

  if (ppadvect.contains("chk_prefix"))
    {
      std::string chkprefix;
      ppadvect.query("chk_prefix",chkprefix);
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
  if (!ppadvect.contains("restart_file"))
    {
      if (!ppadvect.contains("fixed_hierarchy"))
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
          amr.setupForNewAMRRun();
        }
      else
        {
          std::string gridFile;
          ppadvect.query("fixed_hierarchy",gridFile);

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
          //          readBoxes(amrGrids, ppadvect, probDomain,
          //                    maxGridSize, blockFactor,
          //                    numLevels, refRatios, verbosity);
          amr.setupForFixedHierarchyRun(amrGrids,1);
        }
    }
  else
    {
      std::string restartFile;
      ppadvect.query("restart_file",restartFile);

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
  // as things stand now, coordSysFact is deleted by amrAdvectFactory's destructor
  delete ibc;
  delete advectionPhysics;
  delete coordSysPtr;

  // delete coordSysFactPtr;
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

#include "NamespaceFooter.H"
