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
#include "GaussianAdvectMultiMappedIBC.H"
#include "PhysAdvectMappedIBC.H"
#include "AMRLevelMappedConsFactory.H"
#include "LevelGridMetrics.H"
#include "MayDay.H"
#include "MOLAdvectionPhysics.H"
#include "AdvectionAbsoluteMappedTaggingStrategy.H"
#include "NoPyAdvectionMappedStabilityStrategy.H"

#include "DoubleCartesianCS.H"
#include "TripleCartesianCS.H"
#include "CylinderEquiangularCS.H"
#include "CylinderSpokesCS.H"
#include "CubedSphere2DCS.H"
#include "XPointCoordSys.H"
#include "CubedSphereSolidCS.H"

#include "AdvectCosineBellIBC.H"
#include "AdvectConstantIBC.H"
#include "DeformationalFlowIBC.H"
#include "GaussianHillsFlowIBC.H"
#include "CosineBellsFlowIBC.H"
#include "SlottedCylindersFlowIBC.H"
#include "CorrelatedCosineBellsFlowIBC.H"
#include "MultiBlockMeshRefine.H"

#include "DebugDump.H"

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
  TRIPLECARTESIAN,
  XPOINT,
  CUBEDSPHERESOLID
};

void amrAdvect();

// amrAdvect is a function (as opposed to inline in main()) to get
// around MPI scoping problems

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
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void amrAdvect()
{
  // Determine the sample problem specified
  //XXX -- not used
  //XXXint problem = -1;
  std::string problemString;

  // Read inputs that are prefixed with "run."
  ParmParse ppRun("run");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppRun.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Parameters specific to different sample problems
  // int realDim = SpaceDim;

  // Stop after this number of steps
  int nstop = 0;
  ppRun.query("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppRun.query("max_time",stopTime);

  // Read inputs that are prefixed with "amrlevelfactory."
  ParmParse ppAmrlevelfactory("amrlevelfactory");

  // order of spatial  accuracy to use (2 or 4)
  int spaceOrder = 4;
  ppAmrlevelfactory.query("space_order", spaceOrder);

  // Initial values are average:  default TRUE
  bool initialAverage = true;
  ppAmrlevelfactory.query("initial_average", initialAverage);

  // limit face values in advection:  default TRUE
  bool limitFaceValues = true;
  ppAmrlevelfactory.query("limit_face_values", limitFaceValues);

  // enforce a min value on advected quantity:  default FALSE
  bool enforceMinVal = false;
  ppAmrlevelfactory.query("enforce_min_val", enforceMinVal);

  Real minVal = 0.0;
  if (enforceMinVal)
    {
      ppAmrlevelfactory.get("min_val", minVal);
    }

  // forward Euler integration:  default FALSE
  bool forwardEuler = false;
  ppAmrlevelfactory.query("forward_euler", forwardEuler);

  // use artificial dissipation:  default FALSE
  bool useArtificialDissipation = false;
  ppAmrlevelfactory.query("use_artificial_dissipation", useArtificialDissipation);

  // coefficient of artificial dissipation
  Real artificialDissipation = 1.;
  ppAmrlevelfactory.query("artificial_dissipation", artificialDissipation);

  //  // use hyperviscosity:  default FALSE
  //  bool useHyperviscosity = false;
  //  ppAmrlevelfactory.query("use_hyperviscosity", useHyperviscosity);

  //  Real hyperviscosity = -1.0;
  //  ppAmrlevelfactory.query("hyperviscosity", hyperviscosity);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppRun.query("domain_length",domainLength);

  // Set the location of the lower left corner
  std::vector<Real> x0a(SpaceDim,0.0);
  ppRun.queryarr("x0",x0a,0,SpaceDim);
  RealVect x0;
  for ( int d=0 ; d<SpaceDim ; ++d ) x0[d] = x0a[d] ;

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  ppRun.queryarr("num_cells",numCells,0,SpaceDim);
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
  ppRun.queryarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      isPeriodic[dim] = (isPeriodica[dim] == 1);
      if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
        pout() << "Using Periodic BCs in direction: " << dim << endl;
    }

  // Read inputs that are prefixed with "amr."
  ParmParse ppAmr("amr");

  // Maximum AMR level limit
  int maxLevel = 0;
  ppAmr.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppAmr.queryarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppAmr.queryarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppAmrlevelfactory.query("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppAmrlevelfactory.query("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppAmr.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppAmr.query("max_grid_size",maxGridSize);
  int maxBaseGridSize = 0;
  ppAmr.query("max_base_grid_size",maxBaseGridSize);

  Real fillRatio = 0.75;
  ppAmr.query("fill_ratio",fillRatio);

  int gridBufferSize = 1;
  ppAmr.query("grid_buffer_size",gridBufferSize);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppAmr.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppAmr.query("plot_interval",plotInterval);

  // Write mapped-grid geometry info:  default TRUE
  bool writeMap = true;
  ppAmrlevelfactory.query("write_map", writeMap);

  // Write J:  default TRUE
  bool writeJ = true;
  ppAmrlevelfactory.query("write_J", writeJ);

  // Write exact solution:  default FALSE
  bool writeExact = false;
  ppAmrlevelfactory.query("write_exact", writeExact);

  // Write error:  default FALSE
  bool writeError = false;
  ppAmrlevelfactory.query("write_error", writeError);

  // CFL multiplier
  Real cfl = 0.8;
  ppAmrlevelfactory.query("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppAmrlevelfactory.query("initial_cfl",initialCFL);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppAmr.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppAmr.query("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppAmr.query("dt_tolerance_factor",dtToleranceFactor);

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

  pout() << "use artificial dissipation = " << useArtificialDissipation << endl;
  pout() << "artificial dissipation = " << artificialDissipation << endl;

  //  pout() << "use hyperviscosity = " << useHyperviscosity << endl;
  //  pout() << "hyperviscosity = " << hyperviscosity << endl;

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

  // Read inputs that are prefixed with "run.ibc."
  ParmParse ppRunIBC("run.ibc");

  ppRunIBC.get("problem", problemString);

  PhysAdvectMappedIBC* ibc = NULL;

  if (problemString == "cosinebelladvection")
    {
      Real ambientDensity = 0.0;
      ppRunIBC.query("ambient_density", ambientDensity);
      Real deltaDensity = 1000.0;
      ppRunIBC.query("delta_density", deltaDensity);
      // center point of cosine bell:
      // (longitude, latitude) in cubed sphere;
      // real coordinates in other coordinate systems.
      RealVect centerPoint;
      vector<Real> centerPointVect(SpaceDim, 0.);
      if (ppRunIBC.queryarr("center_point", centerPointVect, 0, SpaceDim))
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            centerPoint[idir] = centerPointVect[idir];
        }
      Real size = 0.3;
      ppRunIBC.query("size", size);
      Real advectVelocity = M_PI / 6.0;
      ppRunIBC.query("vel", advectVelocity);
      Real advectAngle = 0.0; //M_PI / 2.0;
      ppRunIBC.query("angle", advectAngle);
      Real evalTime = 0.0;
      ppRunIBC.query("time", evalTime);
      AdvectCosineBellIBC* advectibc = new AdvectCosineBellIBC(ambientDensity,
                                                               deltaDensity,
                                                               centerPoint,
                                                               size,
                                                               advectVelocity,
                                                               advectAngle,
                                                               evalTime);
      std::string velType;
      ppRunIBC.query("velType", velType);
      if (velType == "uniform")
        {
          advectibc->setUniformVel();
          vector<Real> velUniformVect(SpaceDim, 0.);
          if (ppRunIBC.queryarr("vel_uniform", velUniformVect, 0, SpaceDim))
            {
              RealVect velUniform = RealVect(velUniformVect);
              advectibc->setVelUniform(velUniform);
            }
        }
      else if (velType == "solidBody")
        { // this is the default anyway
          advectibc->setSolidBodyRotation();
        }

      ibc = advectibc;
    }
  else if (problemString == "constant")
    {
      Real fieldVal = 1.0;
      ppRunIBC.query("field", fieldVal);
      Real advectVelocity = M_PI / 6.0;
      ppRunIBC.query("vel", advectVelocity);
      Real advectAngle = 0.0; //M_PI / 2.0;
      ppRunIBC.query("angle", advectAngle);
      AdvectConstantIBC* advectibc = new AdvectConstantIBC(fieldVal,
                                                           advectVelocity,
                                                           advectAngle);
      ibc = advectibc;
    }
  else if (problemString == "deformationalflow")
    {
      Real ambientDensity = 0.1;
      ppRunIBC.query("ambient_density", ambientDensity);
      Real deltaDensity = 0.9;
      ppRunIBC.query("delta_density", deltaDensity);
      Real size = 0.5;
      ppRunIBC.query("size", size);
      RealVect centerPoint;
      vector<Real> centerPointVect(SpaceDim, 0.);
      if (ppRunIBC.queryarr("center_point", centerPointVect, 0, SpaceDim))
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            centerPoint[idir] = centerPointVect[idir];
        }
      // this is used ONLY in initialization:
      // const double InitialTheta = m_evalTime * m_advectVelocity;
      Real advectVelocity = M_PI / 6.0;
      ppRunIBC.query("vel", advectVelocity);
      Real advectAngle = 0.0; //M_PI / 2.0;
      ppRunIBC.query("angle", advectAngle);
      Real evalTime = 0.0;
      ppRunIBC.query("time", evalTime);

      DeformationalFlowIBC* dfibc = new DeformationalFlowIBC(ambientDensity,
                                                             deltaDensity,
                                                             centerPoint,
                                                             size,
                                                             advectVelocity,
                                                             advectAngle,
                                                             evalTime);

      ibc = dfibc;
    }
  else if (problemString == "gaussian")
    {
      GaussianAdvectMultiMappedIBC* advectibc = new GaussianAdvectMultiMappedIBC;
      std::string velType;
      ppRunIBC.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          ppRunIBC.getarr("vel", vela,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              vel[idir] = vela[idir];
            }
          advectibc->setUniformVel(vel);
        }
      else if (velType == "solidBody")
        {
          Real omega;
          ppRunIBC.get("omega", omega);
          RealVect rotationCenter;
          Vector<Real> rotCtra(SpaceDim);
          ppRunIBC.getarr("rotationCenter", rotCtra,0,SpaceDim);
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
          ppRunIBC.getarr("vel", vela,0,SpaceDim);
          if (ppRunIBC.contains("transvel"))
            {
              ppRunIBC.getarr("transvel", vela,0,SpaceDim);
            }
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              vel[idir] = vela[idir];
            }
          Real oscAmp;
          ppRunIBC.get("oscAmp", oscAmp);
          advectibc->setTranslatingOscillation(vel,oscAmp);
        }
      else
        {
          MayDay::Error("GaussianAdvectMappedIBC -- bad velType");
        }

      Real r0;
      ppRunIBC.query("radius",r0);

      // Set the location of the center of the gaussian
      std::vector<Real> centera(SpaceDim,0.0);
      ppRunIBC.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;

      Real evalTime = 0.0;
      ppRunIBC.query("time", evalTime);

      advectibc->setParams(r0,center,x0,evalTime);

      ibc  = advectibc;
    }
  else if (problemString == "cosinebellsflow")
    {
      Real hmax = 1.;
      ppRunIBC.query("hmax", hmax);
      Real radius = 0.5;
      ppRunIBC.query("radius", radius);
      int bells = 2;
      ppRunIBC.query("bells", bells);

      Vector<Real> longitude(bells);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      ppRunIBC.queryarr("longitude", longitude, 0, bells);

      Vector<Real> latitude(bells, 0.);
      ppRunIBC.queryarr("latitude", latitude, 0, bells);

      Real background = 0.1;
      ppRunIBC.query("background", background);
      Real delta = 0.9;
      ppRunIBC.query("delta", delta);

      Real period = 5.;
      ppRunIBC.query("period", period);
      Real kappa = 2.;
      ppRunIBC.query("kappa", kappa);
      Real evalTime = 0.0;
      ppRunIBC.query("time", evalTime);

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
      int ncoeffs = 3;
      ppRunIBC.query("ncoeffs", ncoeffs);
      // coefficients, starting with constant term
      Vector<Real> coeffs(ncoeffs);
      coeffs[0] = 0.9;
      coeffs[1] = 0.;
      coeffs[2] = -0.8;
      ppRunIBC.queryarr("coeffs", coeffs, 0, ncoeffs);

      Real hmax = 1.;
      ppRunIBC.query("hmax", hmax);
      Real radius = 0.5;
      ppRunIBC.query("radius", radius);
      int bells = 2;
      ppRunIBC.query("bells", bells);

      Vector<Real> longitude(bells);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      ppRunIBC.queryarr("longitude", longitude, 0, bells);

      Vector<Real> latitude(bells, 0.);
      ppRunIBC.queryarr("latitude", latitude, 0, bells);

      Real background = 0.1;
      ppRunIBC.query("background", background);
      Real delta = 0.9;
      ppRunIBC.query("delta", delta);

      Real period = 5.;
      ppRunIBC.query("period", period);
      Real kappa = 2.;
      ppRunIBC.query("kappa", kappa);
      Real evalTime = 0.0;
      ppRunIBC.query("time", evalTime);

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
      Real hmax = 0.95;
      ppRunIBC.query("hmax", hmax);
      Real width = 5.;
      ppRunIBC.query("width", width);
      int hills = 2;
      ppRunIBC.query("hills", hills);

      Vector<Real> longitude(hills);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      ppRunIBC.queryarr("longitude", longitude, 0, hills);

      Vector<Real> latitude(hills, 0.);
      ppRunIBC.queryarr("latitude", latitude, 0, hills);

      Real period = 5.;
      ppRunIBC.query("period", period);
      Real kappa = 2.;
      ppRunIBC.query("kappa", kappa);
      Real evalTime = 0.0;
      ppRunIBC.query("time", evalTime);

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
      Real hmax = 1.;
      ppRunIBC.query("hmax", hmax);
      Real radius = 0.5;
      ppRunIBC.query("radius", radius);
      int cylinders = 2;
      ppRunIBC.query("cylinders", cylinders);

      Vector<Real> longitude(cylinders);
      longitude[0] = (5*M_PI) / 6.;
      longitude[1] = (7*M_PI) / 6.;
      ppRunIBC.queryarr("longitude", longitude, 0, cylinders);

      Vector<Real> latitude(cylinders, 0.);
      ppRunIBC.queryarr("latitude", latitude, 0, cylinders);

      Real background = 0.1;
      ppRunIBC.query("background", background);

      Real period = 5.;
      ppRunIBC.query("period", period);
      Real kappa = 2.;
      ppRunIBC.query("kappa", kappa);
      Real evalTime = 0.0;
      ppRunIBC.query("time", evalTime);

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
  else
    {
      MayDay::Error("Invalid problemString");
    }

  MOLAdvectionPhysics* advectionPhysics = new MOLAdvectionPhysics();
  advectionPhysics->setPhysIBC(ibc);

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
  else if (coordSysString == "xpoint")
    csCode = XPOINT;
  else if (coordSysString == "cubedspheresolid")
    csCode = CUBEDSPHERESOLID;

  // Set up the coordinate system... factory
  //
  Real dx = domainLength / (Real(numCells[0]));
  RealVect dxVect = dx * RealVect::Unit;
  Vector<int> ratios(refRatios);
  RefCountedPtr<MultiBlockCoordSysFactory> coordSysFactPtr;
  IntVect levelDomainLo, levelDomainHi;
  Real fullDomainLength = domainLength;
  switch (csCode)
    {
    case CYLINDERSPOKES:
      {
        RefCountedPtr<CylinderSpokesCSFactory> cylinderCSFactPtr =
          RefCountedPtr<CylinderSpokesCSFactory>
          (new CylinderSpokesCSFactory());

        // Read inputs that are prefixed with "cylinder."
        RealVect centralRectSize = 2. * RealVect::Unit;
        vector<Real> centralRectSizeVect(SpaceDim, 0.);
        if (coordSysPP.queryarr("central_rect_size", centralRectSizeVect, 0, SpaceDim))
          {
            for (int idir = 0; idir < SpaceDim; idir++)
              centralRectSize[idir] = centralRectSizeVect[idir];
          }

        Real outerRadius = 3.;
        coordSysPP.query("outer_radius", outerRadius);

        // I HAD this set to centerPoint, but now centerPoint is for IBC.
        cylinderCSFactPtr->setCenterPoint(RealVect::Zero);
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
        fullDomainLength = 5 * domainLength;
        break;
      }
    case CYLINDEREQUIANGULAR:
      {
        RefCountedPtr<CylinderEquiangularCSFactory> cylinderCSFactPtr =
          RefCountedPtr<CylinderEquiangularCSFactory>
          (new CylinderEquiangularCSFactory());

        RealVect centralRectSize = RealVect::Unit;
        vector<Real> centralRectSizeVect(SpaceDim, 0.);
        if (coordSysPP.queryarr("central_rect_size", centralRectSizeVect, 0, SpaceDim))
          {
            for (int idir = 0; idir < SpaceDim; idir++)
              centralRectSize[idir] = centralRectSizeVect[idir];
          }

        Real outerRadius = 1.5;
        coordSysPP.query("outer_radius", outerRadius);

        // I HAD this set to centerPoint, but now centerPoint is for IBC.
        cylinderCSFactPtr->setCenterPoint(RealVect::Zero);
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
        fullDomainLength = 5 * domainLength;
        break;
      }
    case CYLINDERTRANSITION:
      {
        MayDay::Error("CylinderTransition not implemented yet.");
        break;
      }
    case CUBEDSPHERE2D:
      {
        RefCountedPtr<CubedSphere2DCSFactory> cubedSphere2DCSFactPtr =
          RefCountedPtr<CubedSphere2DCSFactory>
          (new CubedSphere2DCSFactory());
        coordSysFactPtr = cubedSphere2DCSFactPtr;
        levelDomainLo = IntVect::Zero;
        levelDomainHi = IntVect(D_DECL6(11*lengthCells-1,
                                        lengthCells-1,
                                        0, 0, 0, 0)); // dimensions > 2 unused
        fullDomainLength = 11. * M_PI / 2. * domainLength;
        // realDim = 3; // but SpaceDim == 2
        break;
      }
    case CUBEDSPHERE2DFLAT:
      {
        RefCountedPtr<CubedSphere2DCSFactory> cubedSphere2DCSFactPtr =
          RefCountedPtr<CubedSphere2DCSFactory>
          (new CubedSphere2DCSFactory());
        cubedSphere2DCSFactPtr->setFlatMap(true);
        coordSysFactPtr = cubedSphere2DCSFactPtr;
        levelDomainLo = IntVect::Zero;
        levelDomainHi = IntVect(D_DECL6(11*lengthCells-1,
                                        lengthCells-1,
                                        0, 0, 0, 0)); // dimensions > 2 unused
        fullDomainLength = 11. * M_PI / 2. * domainLength;
        break;
      }
    case DOUBLECARTESIAN:
      {
        RefCountedPtr<DoubleCartesianCSFactory> doubleCartesianCSFactPtr =
          RefCountedPtr<DoubleCartesianCSFactory>
          (new DoubleCartesianCSFactory());
        coordSysFactPtr = doubleCartesianCSFactPtr;
        levelDomainLo = IntVect::Zero;
        levelDomainHi = (3*lengthCells-1) * IntVect::Unit;
        fullDomainLength = 3 * domainLength;
        break;
      }
    case TRIPLECARTESIAN:
      {
        RefCountedPtr<TripleCartesianCSFactory> tripleCartesianCSFactPtr =
          RefCountedPtr<TripleCartesianCSFactory>
          (new TripleCartesianCSFactory());
        coordSysFactPtr = tripleCartesianCSFactPtr;
        levelDomainLo = IntVect::Zero;
        levelDomainHi = (5*lengthCells-1) * IntVect::Unit;
        fullDomainLength = 5 * domainLength;
        break;
      }
    case XPOINT:
      {
        RefCountedPtr<XPointCoordSysFactory> xpointCoordSysFactPtr =
          RefCountedPtr<XPointCoordSysFactory>
          (new XPointCoordSysFactory());
        coordSysFactPtr = xpointCoordSysFactPtr;
        levelDomainLo = IntVect::Zero;
        levelDomainHi = IntVect(D_DECL6(15*lengthCells-1,
                                        lengthCells-1,
                                        0, 0, 0, 0)); // dimensions > 2 unused
        fullDomainLength = 15 * domainLength;
        break;
      }
    case CUBEDSPHERESOLID:
      {
        RefCountedPtr<CubedSphereSolidCSFactory> cubedSphereSolidCSFactPtr =
          RefCountedPtr<CubedSphereSolidCSFactory>
          (new CubedSphereSolidCSFactory());

        Real r0 = 1.;
        coordSysPP.query("r0", r0);

        Real r1 = 3.;
        coordSysPP.query("r1", r1);

        cubedSphereSolidCSFactPtr->setR0(r0);
        cubedSphereSolidCSFactPtr->setR1(r1);

        coordSysFactPtr = cubedSphereSolidCSFactPtr;
        levelDomainLo = -2*lengthCells * IntVect::Unit;
        levelDomainHi = (3*lengthCells-1) * IntVect::Unit;
        fullDomainLength = 5 * domainLength;
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

  // This computes a stable time step.
  Real stabilityFactor = 1.3925;
  NoPyAdvectionMappedStabilityStrategy* stabilityStrategy =
    new NoPyAdvectionMappedStabilityStrategy(stabilityFactor);

  AdvectionAbsoluteMappedTaggingStrategy* taggingStrategy =
    new AdvectionAbsoluteMappedTaggingStrategy(refineThresh);

  // Set up the AMRLevel... factory.
  // Its destructor will delete stabilityStrategy and taggingStrategy.

  AMRLevelMappedConsFactory amrAdvectFact(&(*coordSysFactPtr),
                                          stabilityStrategy,
                                          taggingStrategy);

  amrAdvectFact.spaceOrder(spaceOrder);
  amrAdvectFact.limitFaceValues(limitFaceValues);
  amrAdvectFact.initialAverage(initialAverage);
  amrAdvectFact.enforceMinVal(enforceMinVal, minVal);
  amrAdvectFact.useArtificialDissipation(useArtificialDissipation);
  amrAdvectFact.artificialDissipation(artificialDissipation);
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
  amrAdvectFact.writeJ(writeJ);
  amrAdvectFact.writeExact(writeExact);
  amrAdvectFact.writeError(writeError);
  amrAdvectFact.writeMap(writeMap);
  // amrAdvectFact.IBC(ibc);
  // amrAdvectFact.coordinateSystemFactory(coordSysFactPtr);

  // Set up output files -- do this with AMRLevelFactory
  // in order to be able to do mapped-grid output as well.
  std::string prefix;
  if (ppAmr.contains("plot_prefix"))
    {
      ppAmr.query("plot_prefix",prefix);
      amrAdvectFact.plotPrefix(prefix);
    }

  RefCountedPtr<MultiBlockMeshRefine> mbmr =
    RefCountedPtr<MultiBlockMeshRefine>
    (new MultiBlockMeshRefine(levelDomain,
                              dxVect,
                              coordSysFactPtr,
                              ratios,
                              fillRatio,
                              blockFactor,
                              gridBufferSize,
                              maxGridSize));

  { // scope of AMR amr;
  AMR amr;

  amr.setMeshRefine(mbmr);
  // Set up the AMR object
  amr.define(maxLevel, refRatios, levelDomain, &amrAdvectFact);


  if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }

  if (ppAmr.contains("plot_prefix"))
    {
      ppAmr.query("plot_prefix",prefix);
      amr.plotPrefix(prefix);
      amrAdvectFact.plotPrefix(prefix);
    }

  if (ppAmr.contains("chk_prefix"))
    {
      std::string chkprefix;
      ppAmr.query("chk_prefix",chkprefix);
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
  amr.gridBufferSize(gridBufferSize);

  // Set output parameters
  amr.checkpointInterval(checkpointInterval);
  amr.plotInterval(plotInterval);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(maxDtGrowth);
  amr.dtToleranceFactor(dtToleranceFactor);


  amr.verbosity(verbosity);

  // Set up input files
  if (!ppAmr.contains("restart_file"))
    {
      if (!ppAmr.contains("fixed_hierarchy"))
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
          amr.setupForNewAMRRun();
        }
      else
        {
          std::string gridFile;
          ppAmr.query("fixed_hierarchy",gridFile);

          int numLevels = maxLevel+1;
          // initialize from a list of grids in "gridFile"
          Vector<Vector<Box> > amrGrids(numLevels);
          amrGrids[0] = allBoxes;
          amr.setupForFixedHierarchyRun(amrGrids,1);
        }
    }
  else
    {
      std::string restartFile;
      ppAmr.query("restart_file",restartFile);

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
  delete ibc;
  delete advectionPhysics;
  delete coordSysPtr;
  // as things stand now, coordSysPtr is deleted by amrAdvectFact's destructor
  // delete coordSysFactPtr;
}

#include "NamespaceFooter.H"
