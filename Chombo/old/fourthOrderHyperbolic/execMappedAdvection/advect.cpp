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
#include "AMRLevel.H"
#include "AMRLevelAdvectFactory.H"

#include "GaussianIBC.H"
#include "ConstantIBC.H"
#include "CompactSupportIBC.H"
#include "SinusoidalIBC.H"
#include "CircleIBC.H"
#include "TopHatIBC.H"
#include "SlottedCircleIBC.H"
#include "PolynomialPatchIBC.H"

#include "CartesianCoordSys.H"
//#include "AffineCoordSys.H"
#include "RThetaZCoordSys.H"
#include "TwistedCoordSys.H"
#include "WarpedCoordSys.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#ifdef USE_ARRAYVIEW

#include "UsingNamespace.H"

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

#ifdef CH_USE_TIMER
Chombo::Timer *all_timer ,*setup_timer ,*solve_timer ,*timestep_timer ,*shutdown_timer ;
#endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
  setChomboMPIErrorHandler();
#endif

#ifdef CH_USE_TIMER
  // timers
  all_timer = new Chombo::Timer("All",0) ;
  setup_timer = new Chombo::Timer("Setup",*all_timer) ;
  solve_timer = new Chombo::Timer("Solve",*all_timer,1) ;
  timestep_timer = new Chombo::Timer("TimeStep",*solve_timer) ;
  shutdown_timer = new Chombo::Timer("Shutdown",*all_timer) ;

  Chombo::Timer::TimerInit(0);
  all_timer->start();
  setup_timer->start();
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

#ifdef CH_USE_TIMER
  shutdown_timer->stop() ;
  all_timer->stop() ;
  Chombo::Timer::TimerSummary();
#endif

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

  // Stop after this number of steps
  int nstop = 0;
  ppadvect.query("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppadvect.query("max_time",stopTime);

  // order of spatial  accuracy to use (2 or 4)
  int spaceOrder = 4;
  ppadvect.query("spaceOrder", spaceOrder);

  // limit face values in advection
  bool limitFaceValues = false;
  ppadvect.query("limitFaceValues", limitFaceValues);

  // if true, enforce a min value on advected quantity
  bool enforceMinVal = false;
  ppadvect.query("enforceMinVal", enforceMinVal);

  Real minVal = 0.0;
  std::string lowOrderFluxScheme("CTU");
  if (enforceMinVal)
    {
      ppadvect.get("minVal", minVal);
      ppadvect.get("lowOrderFluxScheme", lowOrderFluxScheme);
    }

  bool redistributeNegativeVal(false);
  ppadvect.query("redistributeNegativeVal", redistributeNegativeVal);
  int maxRedistributionPasses(4);
  if (redistributeNegativeVal)
    {
      ppadvect.get("maxRedistributionPasses", maxRedistributionPasses);
      CH_assert( maxRedistributionPasses>0 );
    }

  // use hyperviscosity
  bool useHyperviscosity = false;
  ppadvect.query("useHyperviscosity", useHyperviscosity);

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

 CH_assert(D_TERM(   (numCells[0] > 0),
                && (numCells[1] > 0),
                && (numCells[2] > 0)));
 CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                && (numCells[1] % 2 == 0),
                && (numCells[2] % 2 == 0)));

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

  pout() << "limit face values = " << limitFaceValues << endl;

  pout() << "enforceMinVal = " << enforceMinVal;
  if (enforceMinVal)
    {
      pout() << "; minVal = " << minVal;
      pout() << "; lowOrderFluxScheme = " << lowOrderFluxScheme;
    }
  pout() << endl;

  pout() << "redistributeNegativeVal = " << redistributeNegativeVal;

  if (redistributeNegativeVal)
    {
      pout() << "; maxRedistributionPasses = " << maxRedistributionPasses;
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

  ppadvect.get("problem", problemString);

  BasicIBC* ibc = NULL;
  if (problemString == "gaussian")
    {
      GaussianIBC* advectibc = new GaussianIBC;
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
          MayDay::Error("GaussianIBC -- bad velType");
        }

      Real r0;
      gaussianPP.query("radius",r0);

      // Set the location of the center of the gaussian
      std::vector<Real> centera(SpaceDim,0.0);
      gaussianPP.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;



      advectibc->setParams(r0,center,x0);

      ibc  = advectibc;
    }
  else if (problemString == "circle")
    {
      CircleIBC* advectibc = new CircleIBC;
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
          MayDay::Error("CircleIBC -- bad velType");
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
  else if (problemString == "tophat")
    {
      TopHatIBC* advectibc = new TopHatIBC;
      ParmParse tophatPP("tophat");
      std::string velType;
      tophatPP.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          tophatPP.getarr("vel", vela,0,SpaceDim);
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)

          advectibc->setUniformVel(vel);
        }
      else if (velType == "solidBody")
        {
          Real omega;
          tophatPP.get("omega", omega);
          RealVect rotationCenter;
          Vector<Real> rotCtra(SpaceDim);
          tophatPP.getarr("rotationCenter", rotCtra,0,SpaceDim);
          D_TERM(rotationCenter[0]=rotCtra[0];,
                 rotationCenter[1]=rotCtra[1];,
                 rotationCenter[2]=rotCtra[2];)

         advectibc->setSolidBodyRotation(rotationCenter, omega);
        }
      else
        {
          MayDay::Error("TopHatIBC -- bad velType");
        }

      Real r0;
      tophatPP.query("radius",r0);

      // Set the location of the center of the top hat
      std::vector<Real> centera(SpaceDim,0.0);
      tophatPP.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;

      advectibc->setParams(r0,center,x0);

      ibc  = advectibc;
    }
  else if (problemString == "slottedCircle")
    {
      SlottedCircleIBC* advectibc = new SlottedCircleIBC;
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
          MayDay::Error("CircleIBC -- bad velType");
        }

      Real r0;
      circlePP.get("radius",r0);

      // Set the location of the center of the circle
      std::vector<Real> centera(SpaceDim,0.0);
      circlePP.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;

      Real slotWidth;
      circlePP.get("slotWidth", slotWidth);

      Real slotLength;
      circlePP.get("slotLength", slotLength);

      int slotDir = 0;
      circlePP.query("slotDir", slotDir);

      advectibc->setParams(r0,center,
                           slotWidth, slotLength, slotDir, x0);

      ibc  = advectibc;
    }
  else if (problemString == "constant")
    {
      ConstantIBC* advectibc = new ConstantIBC;
      ParmParse constantPP("constant");
      std::string velType;
      constantPP.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          constantPP.getarr("vel", vela,0,SpaceDim);
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)

          advectibc->setUniformVel(vel);
        }
      else if (velType == "solidBody")
        {
          Real omega;
          constantPP.get("omega", omega);
          RealVect rotationCenter;
          Vector<Real> rotCtra(SpaceDim);
          constantPP.getarr("rotationCenter", rotCtra,0,SpaceDim);
          D_TERM(rotationCenter[0]=rotCtra[0];,
                 rotationCenter[1]=rotCtra[1];,
                 rotationCenter[2]=rotCtra[2];)

         advectibc->setSolidBodyRotation(rotationCenter, omega);
        }
      else if (velType == "transosc")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          constantPP.getarr("vel", vela,0,SpaceDim);
          if (constantPP.contains("transvel"))
            {
              constantPP.getarr("transvel", vela,0,SpaceDim);
            }
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)
          Real oscAmp;
          constantPP.get("oscAmp", oscAmp);
             advectibc->setTranslatingOscillation(vel,oscAmp);
        }
      else
        {
          MayDay::Error("ConstantIBC -- bad velType");
        }

      Real mag;
      constantPP.query("magnitude",mag);

      advectibc->setParams(mag);

      ibc = advectibc;
    }
  else if (problemString == "compact")
    {
      CompactSupportIBC* advectibc = new CompactSupportIBC;
      ParmParse compactPP("compact");
      std::string velType;
      compactPP.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          compactPP.getarr("vel", vela,0,SpaceDim);
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)

          advectibc->setUniformVel(vel);
        }
      else if (velType == "solidBody")
        {
          Real omega;
          compactPP.get("omega", omega);
          RealVect rotationCenter;
          Vector<Real> rotCtra(SpaceDim);
          compactPP.getarr("rotationCenter", rotCtra,0,SpaceDim);
          D_TERM(rotationCenter[0]=rotCtra[0];,
                 rotationCenter[1]=rotCtra[1];,
                 rotationCenter[2]=rotCtra[2];)

         advectibc->setSolidBodyRotation(rotationCenter, omega);
        }
      else if (velType == "transosc")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          compactPP.getarr("vel", vela,0,SpaceDim);
          if (compactPP.contains("transvel"))
            {
              compactPP.getarr("transvel", vela,0,SpaceDim);
            }
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)
          Real oscAmp;
          compactPP.get("oscAmp", oscAmp);
             advectibc->setTranslatingOscillation(vel,oscAmp);
        }
      else
        {
          MayDay::Error("CompactSupportIBC -- bad velType");
        }

      Real mag;
      compactPP.query("magnitude",mag);

      Real width;
      compactPP.query("width",width);

      // Set the location of the center of the function
      std::vector<Real> centera(SpaceDim,0.0);
      compactPP.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;

      advectibc->setParams(mag,width,center);

      ibc = advectibc;
    }
  else if (problemString == "sinusoidal")
    {
      SinusoidalIBC* advectibc = new SinusoidalIBC;
      ParmParse sinusoidalPP("sinusoidal");
      std::string velType;
      sinusoidalPP.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          sinusoidalPP.getarr("vel", vela,0,SpaceDim);
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)

          advectibc->setUniformVel(vel);
        }
      else if (velType == "transosc")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          sinusoidalPP.getarr("vel", vela,0,SpaceDim);
          if (sinusoidalPP.contains("transvel"))
            {
              sinusoidalPP.getarr("transvel", vela,0,SpaceDim);
            }
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)
          Real oscAmp;
          sinusoidalPP.get("oscAmp", oscAmp);
          advectibc->setTranslatingOscillation(vel,oscAmp);
        }
      else
        {
          MayDay::Error("SinusoidalIBC -- bad velType");
        }

      Real mag;
      sinusoidalPP.query("magnitude",mag);

      // Set the offset of the sinusoidal
      std::vector<Real> offseta(SpaceDim,0.0);
      sinusoidalPP.queryarr("offset",offseta,0,SpaceDim);
      RealVect offset;
      for ( int d=0 ; d<SpaceDim ; ++d ) offset[d] = offseta[d] ;

      // Set the wavenumber of the sinusoidal
      std::vector<Real> wavenumbera(SpaceDim,0.0);
      sinusoidalPP.queryarr("wavenumber",wavenumbera,0,SpaceDim);
      RealVect wavenumber;
      for ( int d=0 ; d<SpaceDim ; ++d ) wavenumber[d] = wavenumbera[d] ;


      advectibc->setParams(mag,wavenumber,offset);

      ibc  = advectibc;
    }
  else if (problemString == "polynomialPatch")
    {
      PolynomialPatchIBC* advectibc = new PolynomialPatchIBC;
      ParmParse polynomialPatchPP("polynomialPatch");
      std::string velType;
      polynomialPatchPP.get("velType", velType);
      if (velType == "uniform")
        {
          RealVect vel;
          Vector<Real> vela(SpaceDim);
          polynomialPatchPP.getarr("vel", vela,0,SpaceDim);
          D_TERM(vel[0]=vela[0];,
                 vel[1]=vela[1];,
                 vel[2]=vela[2];)

          advectibc->setUniformVel(vel);
        }
      else
        {
          MayDay::Error("PolynomialPatchIBC -- bad velType");
        }

      Real rad, mag;
      polynomialPatchPP.query("radius",rad);
      polynomialPatchPP.query("magnitude",mag);

      // Set the location of the center of the polynomialPatch
      std::vector<Real> centera(SpaceDim,0.0);
      polynomialPatchPP.queryarr("center",centera,0,SpaceDim);
      RealVect center;
      for ( int d=0 ; d<SpaceDim ; ++d ) center[d] = centera[d] ;



      advectibc->setParams(rad,center,mag);

      ibc  = advectibc;
    }

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
  RealVect dxVect( D_DECL(dx,dx,dx) );
  Vector<int> ratios(refRatios);
  CoordSysFactory<FArrayBox,FluxBox>* coordSysFact;
  if (coordSysString == "cylindrical")
  {
     coordSysFact= new RThetaZCoordSysFactory(probDomain,
                                              ratios,
                                              dxVect);
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
     coordSysFact = new AffineCoordSysFactory(probDomain,
                                              ratios,
                                              dxVect,
                                              transVect,
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

     coordSysFact = new TwistedCoordSysFactory(probDomain,
                                               ratios,
                                               dxVect,
                                               radius,
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

     coordSysFact = new WarpedCoordSysFactory(probDomain,
                                              ratios,
                                              dxVect,
                                              scale,
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

      coordSysFact = new CartesianCoordSysFactory(probDomain,
                                                  ratios,
                                                  dxVect,
                                                  stretch);
      }

  // Set up the AMRLevel... factory
  //
  AMRLevelAdvectFactory amrAdvectFact;
  amrAdvectFact.spaceOrder(spaceOrder);
  amrAdvectFact.limitFaceValues(limitFaceValues);
  amrAdvectFact.enforceMinVal(enforceMinVal, minVal);
  if (enforceMinVal)
     amrAdvectFact.lowOrderFluxScheme(lowOrderFluxScheme);
  amrAdvectFact.redistributeNegativeVal(redistributeNegativeVal,
                                        maxRedistributionPasses);
  amrAdvectFact.useHyperviscosity(useHyperviscosity);
  if (hyperviscosity>=0)
     amrAdvectFact.hyperviscosity(hyperviscosity);
  amrAdvectFact.CFL(cfl);
  amrAdvectFact.domainLength(domainLength);
  amrAdvectFact.refinementThreshold(refineThresh);
  amrAdvectFact.tagBufferSize(tagBufferSize);
  amrAdvectFact.verbosity(verbosity);
  amrAdvectFact.initialDtMultiplier(initialCFL);
  amrAdvectFact.IBC(ibc);
  amrAdvectFact.coordinateSystemFactory(coordSysFact);


  // Set up output files -- do this with AMRLevelFactory
  // in order to be able to do mapped-grid output as well.
  std::string prefix;
  if (ppadvect.contains("plot_prefix"))
    {
      ppadvect.query("plot_prefix",prefix);
      amrAdvectFact.plotPrefix(prefix);
    }

  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrAdvectFact);

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
    // defaults to maxGridSize
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

#ifdef CH_USE_TIMER
  amr.timer(timestep_timer);
  setup_timer->stop();
#endif

  // Run the computation
#ifdef CH_USE_TIMER
  solve_timer->start();
#endif

  amr.run(stopTime,nstop);

#ifdef CH_USE_TIMER
  solve_timer->stop();
#endif

  // Output the last plot file and statistics
#ifdef CH_USE_TIMER
  shutdown_timer->start() ;
#endif

  amr.conclude();

  // clean up memory
  // as things stand now, coordSysFact is deleted by amrAdvectFactory's destructor
  //delete coordSysFact;
  delete ibc;

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
