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
#include "AMRLevelAdvectMappedFactory.H"

#include "MOLPolytropicPhysics.H"
#include "MOLAdvectionPhysics.H"

#include "CartesianCS.H"
//#include "AffineCS.H"
#include "RThetaZCS.H"
#include "TwistedCS.H"
#include "WarpedCS.H"

#include "GaussianMappedIBC.H"
#include "GaussianAdvectMappedIBC.H"
#include "CircleAdvectMappedIBC.H"
//#include "SinusoidalIBC.H"
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

OldTimer Everything    ("gov Everything", 0);
OldTimer TimeReadInput ("gov Read Input",   Everything);
OldTimer TimeSetupAMR  ("gov Setup AMR",    Everything);
OldTimer TimeRun       ("gov Run",          Everything);
OldTimer TimeConclude  ("gov Conclude",     Everything);

// Possible pressure relationships for the initial condition
#define PRESSURE_ISENTROPIC 0
#define PRESSURE_CONSTANT   1

// amrMOL is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrMOL();

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

  OldTimer::TimerInit(rank);

  Everything.start();

  // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrMOL...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

#ifdef TRAP_FPE
  enableFpExceptions ();
#endif

  // Run amrMOL, i.e., do the computation
  amrMOL();

  Everything.stop();

#ifndef CH_NTIMER
  Real end_memory = get_memory_usage_from_OS();

  pout() << endl
         << "Everything completed --- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << end_memory
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << Everything.wc_time()
         << " sec (wall-clock)" << endl << endl;
#endif

#if !defined(CH_NTIMER) && defined(CH_MPI)
  Real avg_memory, min_memory, max_memory;
  gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

  // OldTimer::TimerSummary();

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

void amrMOL()
{
  // Start timing the reading of the input file
  TimeReadInput.start();

  // Read inputs that are prefixed with "godunov."
  ParmParse ppgodunov("godunov");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppgodunov.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Stop after this number of steps
  int nstop = 0;
  ppgodunov.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppgodunov.get("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppgodunov.get("domain_length",domainLength);

  // Set the location of the lower left corner
  std::vector<Real> x0a(SpaceDim, 0.0);
  ppgodunov.queryarr("x0", x0a, 0, SpaceDim);
  RealVect x0;
  for ( int idir=0 ; idir<SpaceDim ; ++idir ) x0[idir] = x0a[idir] ;

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i)
  {
    numCells[i] = 0;
  }
  ppgodunov.getarr("num_cells",numCells,0,SpaceDim);

  CH_assert(D_TERM(   (numCells[0] > 0),
                   && (numCells[1] > 0),
                   && (numCells[2] > 0)));
  CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                   && (numCells[1] % 2 == 0),
                   && (numCells[2] % 2 == 0)));

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

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppgodunov.get("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppgodunov.get("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppgodunov.get("fill_ratio",fillRatio);

  // The hyperbolic codes use a grid buffer of 1
  int gridBufferSize = 1;
  ppgodunov.query("grid_buffer_size", gridBufferSize);

  // Order of the normal predictor (CTU -> 0, PLM -> 1, PPM -> 2)
  // petermc, 10 Oct 2009, set default PPM
  std::string normalPred = "PPM";
  ppgodunov.query("normal_predictor",normalPred);
  int normalPredOrder;
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
  ppgodunov.query("use_fourth_order_slopes",inFourthOrderSlopes);
  useFourthOrderSlopes = (inFourthOrderSlopes == 1);

  // Arbitrary fiat by petermc, 17 June 2008
  useFourthOrderSlopes = true;

  // Do slope limiting:  default true
  int inPrimLimiting = 1;
  bool usePrimLimiting;
  ppgodunov.query("use_prim_limiting",inPrimLimiting);
  usePrimLimiting = (inPrimLimiting == 1);

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

  // Do slope limiting using characteristics.  NOT USED.
  //  int inCharLimiting = 0;
  //  bool useCharLimiting;
  //  ppgodunov.get("use_char_limiting",inCharLimiting);
  //  useCharLimiting = (inCharLimiting == 1);

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

  // Compute timestep from cell-centered velocities:  default false
  int inDtFromCells = 0;
  bool dtFromCells;
  ppgodunov.query("dt_from_cells", inDtFromCells);
  dtFromCells = (inDtFromCells == 1);

  // End timing the reading of the input file
  TimeReadInput.stop();

#ifndef CH_NTIMER
  pout() << "Input Read completed --- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage_from_OS()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeReadInput.wc_time()
         << " sec (wall-clock)" << endl;
#endif

  // Start timing AMR solver setup
  TimeSetupAMR.start();

  // Create and define IBC (initial and boundary condition) object
  PhysMappedIBC* ibc;

  // A minimum pressure needed to construct PolytropicPhysics - used in slope
  // flattening
  Real smallPressure;

  // Don't use source term by default
  bool useSourceTerm = false;

  // Source term multiplier
  Real sourceTermScaling = 0.0;

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  // Create coordinate system factory and initialize
  ParmParse coordSysPP("coordsys");
  string coordSysString("cartesian");
  coordSysPP.get("type", coordSysString);

  // Set up the coordinate system... factory
  //
  Real dx = domainLength / probDomain.domainBox().longside();
  RealVect dxVect( RealVect(D_DECL(dx,dx,dx)) );
  Vector<int> ratios(refRatios);
  NewCoordSysFactory* coordSysFact;
  if (coordSysString == "cylindrical")
  {
    coordSysFact= new RThetaZCSFactory();
  }
  else if (coordSysString == "affine")
  {
     RealVect transVect;
     std::vector<Real> b( SpaceDim, 0.0 );
     coordSysPP.getarr("translation_vector", b, 0, SpaceDim );
     D_TERM(transVect[0]=b[0];,
            transVect[1]=b[1];,
            transVect[2]=b[2];)

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
     CH_assert(2*radius<=domainLength);
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
     D_TERM(scale[0]=b[0];,
            scale[1]=b[1];,
            scale[2]=b[2];)

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
  else
    {
      RealVect stretch(RealVect::Unit);
      if (coordSysPP.contains("stretch"))
        {
          std::vector<Real> tempvect(SpaceDim, 1);
          coordSysPP.getarr("stretch", tempvect, 0, SpaceDim);
          for (int dir=0; dir<SpaceDim; dir++) stretch[dir] = tempvect[dir];
        }

      coordSysFact = new CartesianCSFactory(stretch);
    }

  // Cast to physics base class pointer for technical reasons

  bool isAdvection = false;

  MOLPhysics* molPhysics;

  // Determine the sample problem specified
  std::string problemString;
  if (ppgodunov.contains("problem"))
    {
      ppgodunov.query("problem",problemString);

      // Print some parameters
      if (verbosity >= 2)
        {
          pout() << "problem = " << problemString << endl;
        }

      if (problemString == "gaussiansmooth")
        {
          // For all gas dynamics
          Real gamma = 1.4;
          ppgodunov.get("gamma",gamma);
          pout() << "gamma = " << gamma << endl;

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
              pout() << "Invalid pressure, \""
                     << pressureString
                     << "\", specified in input file"
                     << endl;
              return;
            }

          vector<Real> centerpp(SpaceDim,0.5);
          RealVect center;
          ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
          for (int i = 0; i < SpaceDim; i++)
            {
              center[i] = centerpp[i];
            }

          // for smoothing:  width of physical domain
          vector<Real> widthpp(SpaceDim,1.);
          RealVect width;
          ppgodunov.queryarr("initial_width",widthpp,0,SpaceDim);
          for (int i = 0; i < SpaceDim; i++)
            {
              width[i] = widthpp[i];
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
              pout() << "initial_center = " << D_TERM(center[0] << "  " <<,
                                                      center[1] << "  " <<,
                                                      center[2] << ) endl;
              pout() << "initial_width = " << D_TERM(width[0] << "  " <<,
                                                     width[1] << "  " <<,
                                                     width[2] << ) endl;
              pout() << "initial_size = " << size << endl;
              pout() << "initial_velocity = " << D_TERM(velocity[0] << "  " <<,
                                                        velocity[1] << "  " <<,
                                                        velocity[2] << ) endl;
            }

          // Define IBC for gaussian problem
          GaussianMappedIBC* gaussianibc = new GaussianMappedIBC;
          gaussianibc->setFortranCommon(smallPressure,
                                        gamma,
                                        ambientDensity,
                                        deltaDensity,
                                        pressure,
                                        center,
                                        width,
                                        size,
                                        velocity,
                                        artificialViscosity);
          // Need this before calling GaussianMappedIBC.
          // Also need gaussibc->setCoordSys.
          gaussianibc->setTime(0.);

          ibc = gaussianibc;

          // Set up the physics for polytropic gas dynamics
          MOLPolytropicPhysics* polytropicPhysics =
            new MOLPolytropicPhysics(smallPressure);
          polytropicPhysics->setPhysIBC(ibc);
          molPhysics = static_cast<MOLPhysics*> (polytropicPhysics);
        }
      else if (problemString == "gaussianadvection")
        {
          isAdvection = true;
          GaussianAdvectMappedIBC* advectibc = new GaussianAdvectMappedIBC;
          ParmParse gaussianPP("gaussian");
          std::string velType;
          gaussianPP.get("velType", velType);
          if (velType == "uniform")
            {
              RealVect vel;
              Vector<Real> vela(SpaceDim);
              gaussianPP.getarr("vel", vela,0,SpaceDim);
              D_TERM(vel[0]=vela[0];,
                     vel[1]=vela[1];,
                     vel[2]=vela[2];)

                advectibc->setUniformVel(vel);
            }
          else if (velType == "solidBody")
            {
              Real omega;
              gaussianPP.get("omega", omega);
              RealVect rotationCenter;
              Vector<Real> rotCtra(SpaceDim);
              gaussianPP.getarr("rotationCenter", rotCtra,0,SpaceDim);
              D_TERM(rotationCenter[0]=rotCtra[0];,
                     rotationCenter[1]=rotCtra[1];,
                     rotationCenter[2]=rotCtra[2];)

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

              D_TERM(vel[0]=vela[0];,
                     vel[1]=vela[1];,
                     vel[2]=vela[2];)
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

          advectibc->setParams(r0,center,x0);

          ibc = advectibc;

          // Set up the physics for advection
          MOLAdvectionPhysics* advectionPhysics = new MOLAdvectionPhysics();
          advectionPhysics->setPhysIBC(ibc);
          molPhysics = static_cast<MOLPhysics*> (advectionPhysics);
        }
      else if (problemString == "circleadvection")
        {
          isAdvection = true;
          CircleAdvectMappedIBC* advectibc =  new CircleAdvectMappedIBC;
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

          ibc = advectibc;

          // Set up the physics for advection
          MOLAdvectionPhysics* advectionPhysics = new MOLAdvectionPhysics();
          advectionPhysics->setPhysIBC(ibc);
          molPhysics = static_cast<MOLPhysics*> (advectionPhysics);
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

  // Set up the AMRLevel... factory
  AMRLevelMappedConsFactory* amrLevelFactPtr;
  if (isAdvection)
    {
      amrLevelFactPtr = new AMRLevelAdvectMappedFactory();
    }
  else
    {
      amrLevelFactPtr = new AMRLevelMappedConsFactory();
    }
  amrLevelFactPtr->spaceOrder(spaceOrder);
  amrLevelFactPtr->limitFaceValues(usePrimLimiting); // WAS limitFaceValues
  amrLevelFactPtr->initialAverage(initialAverage);
  amrLevelFactPtr->useFlattening(useFlattening);
  amrLevelFactPtr->noPPM(noPPM);
  amrLevelFactPtr->doDeconvolution(doDeconvolution);
  amrLevelFactPtr->doFaceDeconvolution(doFaceDeconvolution);
  amrLevelFactPtr->useArtificialViscosity(useArtificialViscosity);
  amrLevelFactPtr->artificialViscosity(artificialViscosity);
  amrLevelFactPtr->useArtVisc(useArtVisc);
  amrLevelFactPtr->ratioArtVisc(ratioArtVisc);
  amrLevelFactPtr->forwardEuler(forwardEuler);
  // amrLevelFactPtr->enforceMinVal(enforceMinVal, minVal);
  amrLevelFactPtr->CFL(cfl);
  amrLevelFactPtr->domainLength(domainLength);
  amrLevelFactPtr->refinementThreshold(refineThresh);
  amrLevelFactPtr->tagBufferSize(tagBufferSize);
  amrLevelFactPtr->verbosity(verbosity);
  amrLevelFactPtr->initialDtMultiplier(initialCFL);
  // amrLevelFactPtr->IBC(ibc);
  amrLevelFactPtr->molPhysics(molPhysics);
  amrLevelFactPtr->dtFromCells(dtFromCells);
  amrLevelFactPtr->coordinateSystemFactory(coordSysFact);
  // molPhysics,
  // normalPredOrder,
  // useCharLimiting,
  // useFlattening,
  // useArtificialViscosity,
  // artificialViscosity,
  // useSourceTerm,
  // sourceTermScaling);

  // Set up output files -- do this with AMRLevelFactory
  // in order to be able to do mapped-grid output as well.
  std::string prefix;
  if (ppgodunov.contains("plot_prefix"))
    {
      ppgodunov.query("plot_prefix",prefix);
      amrLevelFactPtr->plotPrefix(prefix);
    }

  { // scope of AMR amr;
  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel, refRatios, probDomain, amrLevelFactPtr);

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
  if (ppgodunov.contains("plot_prefix"))
  {
    ppgodunov.query("plot_prefix",prefix);
    amr.plotPrefix(prefix);
  }

  if (ppgodunov.contains("chk_prefix"))
  {
    ppgodunov.query("chk_prefix",prefix);
    amr.checkpointPrefix(prefix);
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
      readBoxes(amrGrids, ppgodunov, probDomain,
                maxGridSize, blockFactor,
                numLevels, refRatios, verbosity);
      amr.setupForFixedHierarchyRun(amrGrids,1);
    }
  }
  else
  {
    std::string restartFile;
    ppgodunov.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrMOL restart only defined with hdf5");
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

#ifndef CH_TIMER
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
  }
  delete amrLevelFactPtr;
  delete molPhysics;
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
