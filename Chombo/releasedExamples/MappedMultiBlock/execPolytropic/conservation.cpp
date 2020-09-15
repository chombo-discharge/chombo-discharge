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
#include "GaussianMappedIBC.H"
#include "CubedSpherePolyIBC.H"
#include "CubedSphereShellCS.H"
#include "MOLMappedPolytropicPhysics.H"


#include "RampMappedIBC.H"
#include "GaussianMappedIBC.H"
#include "TrigMappedIBC.H"

#include "GaussianMappedMultiBlockIBC.H"

#include "DoubleCartesianCS.H"
#include "TripleCartesianCS.H"
#include "CylinderEquiangularCS.H"
#include "CubedSphere2DCS.H"
#include "SchwarzChristoffelRampCS.H"
#include "SingleBlockCSAdaptor.H"

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

#include "UsingNamespace.H"

#ifdef CH_Linux
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
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
    RAMP
  };

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
      pout() << "Usage:  conservation...ex <inputfile>" << endl;
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
  /*
    Read inputs that are prefixed with "run."
  */
  ParmParse ppRun("run");

  // This determines the amount of diagnostic output generated
  int verbosity = 0;
  ppRun.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Stop after this number of steps
  int nstop = 0;
  ppRun.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppRun.get("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppRun.query("domain_length",domainLength);

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i)
    {
      numCells[i] = 0;
    }
  ppRun.getarr("num_cells",numCells,0,SpaceDim);
  // Check that every component of numCells is positive and even
  for (int i = 0; i < SpaceDim; ++i)
    {
      CH_assert(numCells[i] > 0);
      CH_assert(numCells[i] % 2 == 0);
    }

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppRun.getarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim = 0; dim < SpaceDim; dim++)
    {
      isPeriodic[dim] = (isPeriodica[dim] == 1);
      if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
        {
          pout() << "Using Periodic BCs in direction: " << dim << endl;
        }
    }

  /*
    Read inputs that are prefixed with "amr."
    These are used in the AMR class.
  */
  ParmParse ppAmr("amr");

  // Maximum AMR level limit
  int maxLevel = 0;
  ppAmr.get("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppAmr.getarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppAmr.getarr("regrid_interval",regridIntervals,0,numReadLevels);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppAmr.get("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppAmr.get("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppAmr.get("fill_ratio",fillRatio);

  // Grid buffer size
  int gridBufferSize = 1;
  const int codeReqBufferSize = LevelGridMetrics::bufferSize4thO(
                                                                 refRatios,
                                                                 maxLevel,
                                                                 5);  // Num ghost 5 is hard coded in AMRLevelMappedCons
  const int haveUserBufferSize =
    ppAmr.query("grid_buffer_size", gridBufferSize);
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

  // Set up checkpointing
  int checkpointInterval = 0;
  ppAmr.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppAmr.query("plot_interval",plotInterval);

  // Set up whether to use subcycling in time.
  bool useSubcycling = true;
  ppAmr.query("use_subcycling", useSubcycling);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppAmr.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppAmr.get("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppAmr.get("dt_tolerance_factor",dtToleranceFactor);

  // Set up output files
  std::string plotPrefix;
  if (ppAmr.contains("plot_prefix"))
    {
      ppAmr.query("plot_prefix",plotPrefix);
    }

  std::string chkPrefix;
  if (ppAmr.contains("checkpoint_prefix"))
    {
      ppAmr.query("checkpoint_prefix",chkPrefix);
    }

  /*
    Read inputs that are prefixed with "amrlevelfactory."
    These are used in the AMRLevelFactory class.
  */
  ParmParse ppAmrlevelfactory("amrlevelfactory");

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppAmrlevelfactory.get("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppAmrlevelfactory.get ("refine_thresh",refineThresh);

  // Whether refinement threshold is scaled with dx
  bool refinementIsScaled = false;
  ppAmrlevelfactory.query("refinement_is_scaled", refinementIsScaled);

  // Whether to tag on pressure instead of on density
  bool tagPressure = false;
  ppAmrlevelfactory.query("tag_pressure", tagPressure);

  // Whether to tag on vorticity instead of on density
  bool tagVorticity = false;
  ppAmrlevelfactory.query("tag_vorticity", tagVorticity);

  // Name of tagging strategy
  std::string tagStrategyName("DensityGradient");
  ppAmrlevelfactory.query("tag_strategy", tagStrategyName);

  // Use fourth order slopes:  default true
  bool useFourthOrderSlopes = true;
  ppAmrlevelfactory.query("use_fourth_order_slopes", useFourthOrderSlopes);

  // Do slope limiting:  default true
  bool usePrimLimiting = true;
  ppAmrlevelfactory.query("use_prim_limiting", usePrimLimiting);

  // This should actually be 1 even if slope limiting is off
  bool highOrderLimiter = true;
  ppAmrlevelfactory.query("high_order_limiter", highOrderLimiter);

  // NEW Kreiss-Oliger artificial viscosity
  bool useArtVisc = false;
  ppAmrlevelfactory.query("use_art_visc", useArtVisc);

  Real ratioArtVisc = 0.;
  if (useArtVisc)
    {
      ppAmrlevelfactory.query("ratio_art_visc", ratioArtVisc);
    }

  bool forwardEuler = false;
  ppAmrlevelfactory.query("forward_euler", forwardEuler);

  // Initial values are average:  default false
  bool initialAverage = false;
  ppAmrlevelfactory.query("initial_average", initialAverage);

  // Do slope flattening:  default false
  bool useFlattening = false;
  ppAmrlevelfactory.query("use_flattening", useFlattening);

  // Avoid PPM:  default false
  bool noPPM = false;
  ppAmrlevelfactory.query("no_ppm", noPPM);

  // Do deconvolution:  default true
  bool doDeconvolution = true;
  ppAmrlevelfactory.query("do_deconvolution", doDeconvolution);

  // Do face deconvolution:  default true
  bool doFaceDeconvolution = true;
  ppAmrlevelfactory.query("do_face_deconvolution", doFaceDeconvolution);

  // Apply artificial viscosity based on divergence:  default true
  bool useArtificialViscosity = true;
  ppAmrlevelfactory.get("use_artificial_viscosity", useArtificialViscosity);

  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.1;
  if (useArtificialViscosity)
    {
      ppAmrlevelfactory.get("artificial_viscosity", artificialViscosity);
    }

  // CFL multiplier
  Real cfl = 0.8;
  ppAmrlevelfactory.get("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppAmrlevelfactory.get("initial_cfl",initialCFL);

  // Create and define IBC (initial and boundary condition) object
  PhysIBC* ibc = 0;

  // A minimum pressure needed to construct MOLPolytropicPhysics - used in slope
  // flattening.
  // It is returned by PhysMappedIBC functions, and then given to
  // MOLMappedPolytropicPhysics.
  Real smallPressure;

  ProblemDomain probDomain(IntVect::Zero,
                           IntVect(D_DECL6(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1,
                                           numCells[3]-1,
                                           numCells[4]-1,
                                           numCells[5]-1)),
                           isPeriodic);

  ParmParse ppPolytropic("polytropic");

  // For all gas dynamics
  Real gamma = 1.4;
  ppPolytropic.query("gamma",gamma);

  /*--------------------------------------------------------------------*
   * Determine the sample problem specified
   *--------------------------------------------------------------------*/

  // Read inputs that are prefixed with "run.ibc."
  ParmParse ppRunIBC("run.ibc");

  // argh, these are used to construct both IBC and CS for ramp problem
  Real alphaRad, leadLength, rampLength;

  std::string problemString;
  if (ppRunIBC.contains("problem"))
    {
      ppRunIBC.query("problem", problemString);

      // Print some parameters
      if (verbosity >= 2)
        {
          pout() << "problem = " << problemString << endl;
          pout() << "gamma = " << gamma << endl;
        }

      if (problemString == "gaussiansmooth")
        {
          // smoothed Gaussian problem
          Real ambientDensity = 1.4;
          ppRunIBC.get("ambient_density",ambientDensity);

          Real deltaDensity = 0.014;
          ppRunIBC.get("delta_density",deltaDensity);

          int pressure = -1;
          std::string pressureString;
          ppRunIBC.get("initial_pressure",pressureString);
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
          ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
          for (int i = 0; i < SpaceDim; i++)
            {
              center[i] = centerpp[i];
            }

          // for smoothing:  width of physical domain (unused?)
          vector<Real> widthpp(SpaceDim,1.);
          RealVect width;
          ppRunIBC.queryarr("initial_width",widthpp,0,SpaceDim);
          for (int i = 0; i < SpaceDim; i++)
            {
              width[i] = widthpp[i];
            }

          // domainWidth = physical width of domain
          RealVect domainWidth = RealVect::Unit;
          int minWidthDir = domainWidth.minDir(true);
          Real widthMin = domainWidth[minWidthDir];

          Real size = 0.25;
          ppRunIBC.get("initial_size",size);

          Real outerFlatBuffer = 0.;
          Real radmax = 0.;
          int hasRadMaxIB = ppRunIBC.query("initial_buffer", outerFlatBuffer);
          int hasRadMax   = ppRunIBC.query("initial_max_radius", radmax);
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
          ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
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
          GaussianMappedMultiBlockIBC* gaussianibc = new GaussianMappedMultiBlockIBC;
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
      else if (problemString == "trig")
        {
          // trigonometric problem
          Real ambientDensity = 1.4;
          ppRunIBC.get("ambient_density",ambientDensity);

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
      else if (problemString == "cosinebell")
        {
          pout() << "doing cosine bell problem for cubed sphere" << endl;

          // Define IBC for cubed sphere poly 
          Real rho0, pres0, radiusDegrees, delta,  magnitude, angle, shift;
          int  power, numBump;
          ppRunIBC.get("ambient_density", rho0);
          ppRunIBC.get("ambient_pressure", pres0);
          ppRunIBC.get("density_perturbation", delta);
          ppRunIBC.get("bump_radius", radiusDegrees);
          ppRunIBC.get("bump_number", numBump);
          ppRunIBC.get("bump_power" , power);
          ppRunIBC.get("adv_magnitude" , magnitude);
          Vector<Real> longitudeDegrees(numBump), latitudeDegrees(numBump);
          ppRunIBC.getarr("bump_longitude", longitudeDegrees, 0, numBump);
          ppRunIBC.getarr("bump_latitude",   latitudeDegrees, 0, numBump);

          pout() << "ambient density      = " << rho0    << endl;
          pout() << "ambient pressure     = " << pres0   << endl;
          pout() << "density perturbation = " << delta   << endl;
          pout() << "bump radius          = " << radiusDegrees << " degrees" << endl;
          pout() << "bump power           = " << power   << endl;
          pout() << "velocity magnitude   = " << magnitude << endl;
          pout() << "bump number          = " << numBump << endl;
          Real radiusRadians = (M_PI/180.) * radiusDegrees;
          Vector<Real> longitudeRadians(numBump), latitudeRadians(numBump);
          for (int ibump = 0; ibump < numBump ; ibump++)
            {
              pout() << "bump [" << ibump << "] : "
                     << "longitude " << longitudeDegrees[ibump] << " degrees, "
                     << "latitude " << latitudeDegrees[ibump] << " degrees"
                     << endl;
              latitudeRadians[ibump] = (M_PI/180.) * latitudeDegrees[ibump];
            }
          
          CubedSpherePolyIBC* timecube = new CubedSpherePolyIBC(gamma,
                                                                rho0,
                                                                pres0,
                                                                radiusRadians,
                                                                longitudeRadians,
                                                                latitudeRadians,
                                                                delta,
                                                                power,
                                                                magnitude, 
                                                                angle,
                                                                shift);
          ibc = timecube;
        }
      else if (problemString == "mappedramp")
        {
          if (isPeriodic[0] || isPeriodic[1])
            {
              MayDay::Error("Neither x or y boundaries can be periodic");
            }

          // Ramp problem
          Real alpha = 30.0;
          ppRunIBC.get("angle_deg",alpha);

          Real ms = 10.0;
          ppRunIBC.get("shock_mach",ms);

          leadLength = 0.5;
          ppRunIBC.get("lead_length", leadLength);

          rampLength = 4.0;
          ppRunIBC.get("ramp_length", rampLength);

          Real X0Start = -leadLength/2;
          ppRunIBC.get("shock_location", X0Start);

          if (verbosity >= 2)
            {
              pout() << "alpha = " << alpha << endl;
              pout() << "shock_mach = " << ms << endl;
              pout() << "lead_length = " << leadLength << endl;
              pout() << "ramp_length = " << rampLength << endl;
              pout() << "shock_location = " << X0Start << endl;
            }
          
          alphaRad = alpha*Pi/180.;

          // Reset the domain length in mapped space to 1.0
          domainLength = 1.0;
          if (verbosity >= 2)
            {
              pout() << "Computational domainLength fixed to 1.0 for mappedramp problem" << endl;
            }

          // Define IBC for ramp problem
          pout() << problemString << std::endl;
          RampMappedIBC* rampibc = new RampMappedIBC(smallPressure,
                                                     gamma,
                                                     alphaRad,
                                                     ms,
                                                     X0Start,
                                                     artificialViscosity);
          ibc = rampibc;
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
      pout() << "\"run.ibc.problem\" not specified in input file" << endl << endl;
      return;
    }

  /*--------------------------------------------------------------------*
   * Create coordinate system factory and initialize
   *--------------------------------------------------------------------*/

  CoordSysCode csCode;

  ParmParse coordSysPP("coordsys");
  string coordSysString("UNDEFINED");
  coordSysPP.query("type", coordSysString);
  if (coordSysString == "cylinderspokes")
    csCode = CYLINDERSPOKES;
  else if (coordSysString == "cylinderequiangular")
    csCode = CYLINDEREQUIANGULAR;
  else if (coordSysString == "cylindertransition")
    csCode = CYLINDERTRANSITION;
  else if (coordSysString == "cubedsphereshellpanel")
    csCode = CUBEDSPHERE2D;
  else if (coordSysString == "cubedsphere2dflat")
    csCode = CUBEDSPHERE2DFLAT;
  else if (coordSysString == "doublecartesian")
    csCode = DOUBLECARTESIAN;
  else if (coordSysString == "triplecartesian")
    csCode = TRIPLECARTESIAN;
  else if (coordSysString == "ramp")
    csCode = RAMP;

  // Set up the coordinate system... factory
  //
  //  Real dx = domainLength / (Real(numCells[0]));
  Vector<int> ratios(refRatios);
  // Coordinate system for mapped grids
  MultiBlockCoordSysFactory* coordSysFactPtr = NULL;
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
        RealVect centralRectSize = RealVect::Unit;
        vector<Real> centralRectSizeVect(SpaceDim, 0.);
        if (coordSysPP.queryarr("central_rect_size", centralRectSizeVect, 0, SpaceDim))
          {
            for (int idir = 0; idir < SpaceDim; idir++)
              centralRectSize[idir] = centralRectSizeVect[idir];
          }
        Real outerRadius = 1.5;
        coordSysPP.query("outer_radius", outerRadius);
        RealVect centerPoint = RealVect::Zero;
        vector<Real> centerPointVect(SpaceDim, 0.);
        if (coordSysPP.queryarr("center_point", centerPointVect, 0, SpaceDim))
          {
            for (int idir = 0; idir < SpaceDim; idir++)
              centerPoint[idir] = centerPointVect[idir];
          }

        CylinderEquiangularCSFactory* cylinderCSFactPtr =
          new CylinderEquiangularCSFactory;
        cylinderCSFactPtr->setCenterPoint(centerPoint);
        cylinderCSFactPtr->setCentralRectSize(centralRectSize);
        cylinderCSFactPtr->setOuterRadius(outerRadius);
        coordSysFactPtr = cylinderCSFactPtr;
        break;
      }
    case CYLINDERTRANSITION:
      {
        MayDay::Error("CylinderTransition not implemented yet.");
        break;
      }
    case CUBEDSPHERE2D:
      {
        pout() << "cubed sphere shell coordinate system" << endl;
        CubedSphereShellCSFactory* cubedSphere2DCSFactPtr =
          new CubedSphereShellCSFactory();
        cubedSphere2DCSFactPtr->setHeight(1.);
        cubedSphere2DCSFactPtr->setRadius(1.);

        coordSysFactPtr = cubedSphere2DCSFactPtr;
        break;
      }
    case CUBEDSPHERE2DFLAT:
      {
        CubedSphere2DCSFactory* cubedSphere2DCSFactPtr =
          new CubedSphere2DCSFactory;
        coordSysFactPtr = cubedSphere2DCSFactPtr;
        ((CubedSphere2DCSFactory*) coordSysFactPtr)->setFlatMap(true);
        break;
      }
    case DOUBLECARTESIAN:
      {
        DoubleCartesianCSFactory* doubleCartesianCSFactPtr =
          new DoubleCartesianCSFactory;
        coordSysFactPtr = doubleCartesianCSFactPtr;
        break;
      }
    case TRIPLECARTESIAN:
      {
        TripleCartesianCSFactory* tripleCartesianCSFactPtr =
          new TripleCartesianCSFactory;
        coordSysFactPtr = tripleCartesianCSFactPtr;
        break;
      }
    case RAMP:
      {
        SchwarzChristoffelRampCSFactory* rampCSFactPtr =
          new SchwarzChristoffelRampCSFactory(numCells[0],
                                              alphaRad,
                                              leadLength,
                                              rampLength);
        coordSysFactPtr = new SingleBlockCSAdaptorFactory(rampCSFactPtr);
        break;
      }
    default: 
      {
        MayDay::Error("bogus code for coordinate system");
      }
    }


  // The AMRLevelFactory takes responsibility for destruction of coordSysFactPtr

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

      pout() << "slope_order = "
             << (useFourthOrderSlopes ? "4th" : "2nd") << endl;
      pout() << "use_primitive_slope_limiting = "
             << (usePrimLimiting ? "yes" : "no") << endl;
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

      pout() << "checkpoint_interval = " << checkpointInterval << endl;
      pout() << "plot_interval = " << plotInterval << endl;

      pout() << "CFL = " << cfl << endl;
      pout() << "initial_CFL = " << initialCFL << endl;

      pout() << "maximum_dt_growth = " << maxDtGrowth << endl;
      pout() << "dt_tolerance_factor = " << dtToleranceFactor << endl;
    }

  // Set up the physics for polytropic gas dynamics
  Real M0sq=-1.;
  ppPolytropic.query("fourth_order_artificial_viscosity_parameter",M0sq);
  MOLMappedPolytropicPhysics polytropicMappedPhysics(smallPressure);
  polytropicMappedPhysics.setPhysIBC(ibc);
  if (M0sq > 0.)
    {
      polytropicMappedPhysics.setFourthOrderArtificialViscosityParameter(M0sq);
    }
  if (verbosity >= 2)
    {
      pout() << "fourth-order artificial viscosity parameters = " << M0sq
             << endl;
    }
  // Cast to physics base class pointer for technical reasons.
  MOLPhysics* molPhysics = static_cast<MOLPhysics*> (&polytropicMappedPhysics);

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
  AMRLevelMappedConsFactory amrLevelFact(coordSysFactPtr,
                                         stabilityStrategy,
                                         taggingStrategy);
  amrLevelFact.plotPrefix(plotPrefix);
  amrLevelFact.spaceOrder(spaceOrder);
  amrLevelFact.limitFaceValues(usePrimLimiting); // WAS limitFaceValues
  amrLevelFact.highOrderLimiter(highOrderLimiter); // WAS limitFaceValues
  amrLevelFact.initialAverage(initialAverage);
  amrLevelFact.useFlattening(useFlattening);
  amrLevelFact.noPPM(noPPM);
  amrLevelFact.doDeconvolution(doDeconvolution);
  amrLevelFact.doFaceDeconvolution(doFaceDeconvolution);
  amrLevelFact.useArtificialViscosity(useArtificialViscosity);
  amrLevelFact.artificialViscosity(artificialViscosity);
  amrLevelFact.useArtVisc(useArtVisc);
  amrLevelFact.ratioArtVisc(ratioArtVisc);
  amrLevelFact.forwardEuler(forwardEuler);
  amrLevelFact.CFL(cfl);
  amrLevelFact.domainLength(domainLength);
  amrLevelFact.refinementThreshold(refineThresh);
  amrLevelFact.refinementIsScaled(refinementIsScaled);
  amrLevelFact.tagPressure(tagPressure);
  amrLevelFact.tagVorticity(tagVorticity);
  amrLevelFact.tagBufferSize(tagBufferSize);
  amrLevelFact.verbosity(verbosity);
  amrLevelFact.initialDtMultiplier(initialCFL);
  amrLevelFact.molPhysics(molPhysics);

  { // scope of AMR amr;
    AMR amr;

    // Set up the AMR object
    amr.define(maxLevel, refRatios, probDomain, &amrLevelFact);

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
            //      std::string gridFile;
            //      ppAmr.query("fixed_hierarchy",gridFile);

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
            GenFuncs::readBoxes(amrGrids, ppAmr, probDomain,
                                maxGridSize, blockFactor,
                                numLevels, refRatios, verbosity);
            amr.setupForFixedHierarchyRun(amrGrids,1);
          }
      }
    else // read from restart file
      {
        std::string restartFile;
        ppAmr.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
        HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
        // read from checkpoint file
        amr.setupForRestart(handle);
        handle.close();
#else
        MayDay::Error("amrRun restart only defined with hdf5");
#endif
        // This section added by petermc, 20 Jan 2010
        if (ppAmr.contains("fixed_hierarchy"))
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
