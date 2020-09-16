#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


/*******************************************************************************
 *
 *  \file volumeAMRTest.cpp
 *
 *  \brief Testing of grid metrics in an AMR framework
 *
 *  The effect of AMR on the solution are evaluated for mapped grids.  The
 *  solutions is not updated.
 *
 *  There are two refinement patterns and two types of initial conditions.  The
 *  possible refinement patters are selected by the input file:
 *    volumeAMRTest.rotating.inputs - a rotating patch
 *    volumeAMRTest.center.inputs   - a centered patch appears and disappears
 *  If an input file is not specified, this test runs with the rotating patch.
 *  For each case, two types of IC can be specified:
 *    constant                      - tests freestream preservation
 *    linear                        - tests conservation
 *  If an input file is not specified, both tests are run.  Otherwise, the
 *  appropriate test is selected from the input file.
 *
 ******************************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include "ParmParse.H"
#include "parstream.H"
#include "AMR.H"
#include "AMRLevel.H"
#include "BoxIterator.H"
#include "NodeAMRIO.H"
#include "computeNorm.H"
#include "LoadBalance.H"

#include "FourthOrderUtil.H"

#include "AMRLevelVolTest.H"
#include "AMRLevelVolTestFactory.H"
#include "SingleBlockCSAdaptor.H"
#include "VolTestParameters.H"

#include "CartesianCS.H"
#include "RThetaZCS.H"
#include "RThetaPhiCS.H"
#include "TwistedCS.H"
#include "WarpedCS.H"
#include "AnalyticCS.H"

#include "UsingNamespace.H"

/*==============================================================================
 * Definition of global parameters
 *============================================================================*/

//--In VolTestParameters.H

int VT::g_runStatus = 0;
Real VT::g_conservedJU0 = 0.;
Real VT::g_conservedJU1 = 1.;

//--Global variables for handling output:

static const char *pgmname = "volumeAMRTest";
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;


/*==============================================================================
 * Prototypes
 *============================================================================*/

void
trial(const VT::ICType, const int a_scaleMesh, Real* a_errorNorms);

int
parseTestOptions(int argc ,char* argv[]);


/*==============================================================================
 *  Routines for comparing floating point numbers
 *  src/generalFuncs.H
 *============================================================================*/

namespace GenFuncs
{

// Comparison with limit tol^2 as x and y -> 0.
// Returns true if not equal
template <typename T>
inline bool compare(const T &x, const T &y, int prec)
{
  const T tol = std::pow(10., -std::abs(prec));
  return std::fabs(x - y) >
    (std::min(std::fabs(x), std::fabs(y)) + tol)*tol;
}

// Comparison with limit tol as x and y -> 0.
// Return true if not equal
template <typename T>
inline bool compare1(const T &x, const T &y, int prec)
{
  const T tol = std::pow(10., -std::abs(prec));
  return std::fabs(x - y) >
    std::min(std::fabs(x), std::fabs(y))*tol + tol;
}

}  // End of namespace GenFuncs


/*******************************************************************************
 ///////////////////////////////////////////////////////////////////////////////
 *******************************************************************************
 *
 * exe
 *
 ******************************************************************************/


/*==============================================================================
 *
 *  Entry
 *
 *============================================================================*/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif

  const int inFileArg = parseTestOptions(a_argc, a_argv);

  // Check for an input file
  const char* inFile = "volumeAMRTest.inputs";

  int numICType = 2;
  if (inFileArg < a_argc)
    {
      inFile = a_argv[inFileArg];
      if (std::strcmp(inFile, "volumeAMRTest.inputs") != 0)
        {
          numICType = 1;
        }
      // Otherwise this is probably a test run with 'make run'
    }

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  // Parse the command line and the input file
  ParmParse pp(a_argc - (inFileArg + 1), a_argv + inFileArg + 1, NULL, inFile);

  int eekflagsTotal = 0;
#if CH_SPACEDIM >= 2
  VT::ICType icType = VT::ICTypeConstant;
  for (int iICTypeTrial = 0; iICTypeTrial != numICType; ++iICTypeTrial)
    {
      if (numICType == 1)
        {
          std::string icTypeString;
          pp.get("ic_type", icTypeString);
          if (icTypeString == "constant")
            {
              icType = VT::ICTypeConstant;
            }
          else if (icTypeString == "linear")
            {
              icType = VT::ICTypeLinear;
            }
          else
            {
              MayDay::Error("Unknown argument to keyword \"ic_type\"");
            }
        }
      else
        {
          icType = static_cast<VT::ICType>(iICTypeTrial);
        }

      Real trialErrorNorms[3];
      trial(icType, 1, trialErrorNorms);

      // Test results from the trial
#ifdef CH_USE_DOUBLE
      const int consvTol = std::numeric_limits<Real>::digits10 - 1;
#else
      const int consvTol = std::numeric_limits<Real>::digits10 - 1;
#endif

      int eekflag = 0;
      switch (icType)
        {
        case VT::ICTypeConstant:
          eekflag = VT::g_runStatus;
          if (VT::g_runStatus == 0)
            {
              pout() << indent2 << pgmname
                     << ": freestream preservation passed." << endl;
            }
          else
            {
              pout() << indent2 << pgmname
                     << ": freestream preservation failed on "
                     << VT::g_runStatus << " cells." << endl;
            }
          break;
        case VT::ICTypeLinear:
        {
          if (GenFuncs::compare(VT::g_conservedJU0,
                                VT::g_conservedJU1,
                                consvTol))
            {
              eekflag = 1;
              const int savedPrec = pout().precision();
              pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
              pout().precision(std::numeric_limits<Real>::digits10);
              pout() << indent2 << pgmname
                     << ": conservation failed." << endl;
              pout() << indent2 << VT::g_conservedJU0 << " and "
                     << VT::g_conservedJU1 << std::endl;
              pout().setf(std::ios_base::fmtflags(0),
                          std::ios_base::floatfield);
              pout().precision(savedPrec);
            }
          else
            {
              pout() << indent2 << pgmname
                     << ": conservation passed." << endl;
            }
          break;
        }
        }
      if (eekflag)
        {
          ++eekflagsTotal;
        }
    }

#endif
  if (eekflagsTotal == 0)
    {
      pout() << indent << pgmname
           << ": test passed." << endl;
    }
  else
    {
      pout() << indent << pgmname
           << ": test FAILED with error code "
           << eekflagsTotal << endl;
    }
#ifdef CH_MPI
  // Exit MPI
  MPI_Finalize();
#endif
  return eekflagsTotal;
}


/*==============================================================================
 *
 * Run a trial
 *
 *============================================================================*/

void trial(const VT::ICType a_icType, const int a_scaleMesh, Real* a_errorNorms)
{

/*--------------------------------------------------------------------*
 * Input
 *--------------------------------------------------------------------*/

  // Default parser
  ParmParse pp;

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  pp.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

//--Read options that define the problem

  // Solution type given by argument a_icType

  // Refinement type
  std::string probTypeString;
  VT::RefType refType = VT::RefTypeRotating;
  pp.get("refine_type", probTypeString);
  if (probTypeString == "center")
    {
      refType = VT::RefTypeCenter;
    }
  else if (probTypeString == "rotating")
    {
      refType = VT::RefTypeRotating;
    }
  else
    {
      MayDay::Error("Unknown argument to keyword \"refine_type\"");
    }

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      numCells[dir] = 0;
    }
  pp.getarr("num_cells",numCells,0,SpaceDim);
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      numCells[dir] *= a_scaleMesh;
    }

//--Read options necessary to initialize the amr level factory

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  pp.get("tag_buffer_size",tagBufferSize);

//--Read options necessary to define the amr class

  // Maximum AMR level limit
  int maxLevel = 0;
  pp.get("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  pp.getarr("ref_ratio",refRatios,0,maxLevel+1);
  refRatios[maxLevel] = 1;

//--Read options necessary to setup the amr class

  // Get plot prefix
  std::string plotPrefix("plot");
  if (pp.contains("plot_prefix"))
    {
      pp.query("plot_prefix", plotPrefix);
    }

  // Get checkpoint prefix
  std::string checkpointPrefix("check");
  if (pp.contains("checkpoint_prefix"))
    {
      pp.query("checkpoint_prefix", checkpointPrefix);
    }

  // Set up plot file writing
  int plotInterval = 0;
  pp.query("plot_interval",plotInterval);

  // Set up checkpointing
  int checkpointInterval = 0;
  pp.query("checkpoint_interval",checkpointInterval);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  pp.get("max_grid_size",maxGridSize);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  pp.query("dt_tolerance_factor",dtToleranceFactor);

  //
  Real fillRatio = 0.75;
  pp.get("fill_ratio",fillRatio);

  // Minimum dimension of a grid
  int blockFactor = 1;
  pp.get("block_factor",blockFactor);

  // Requires (int)((numGhost + refRatio - 1)/refRatio) + 2
  int gridBufferSize = 3;
  pp.query("grid_buffer_size", gridBufferSize);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  pp.getarr("regrid_interval",regridIntervals,0,numReadLevels);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  pp.query("max_dt_growth",maxDtGrowth);

  // Determine if a fixed or variable time step will be used.
  Real fixedDt = 1./(100.*a_scaleMesh);  // Denominator is number of steps
  pp.query("fixed_dt",fixedDt);

  // Set up whether to use subcycling in time.
  int useSubcyclingInt = 1;
  pp.query("use_subcycling", useSubcyclingInt);
  bool useSubcycling = (useSubcyclingInt == 1);

//--Options necessary to run the amr class

  // Stop after this number of steps
  int nstop = 0;
  pp.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  pp.get("max_time",stopTime);

/*--------------------------------------------------------------------*
 * Setup the AMR classes
 *--------------------------------------------------------------------*/

//--Set up the problem

  bool isPeriodic[SpaceDim];
  {
    vector<int> isPeriodicI(SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir)
      {
        isPeriodicI[dir] = 0;
      }
    pp.queryarr("is_periodic", isPeriodicI, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir)
      {
        isPeriodic[dir] = (isPeriodicI[dir] == 1);
      }
  }
  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

//--Create coordinate system factory and initialize

  ParmParse coordSysPP("coordsys");
  string coordSysString("cartesian");
  coordSysPP.get("type", coordSysString);

  NewCoordSysFactory* coordSysFact = 0;
  if (coordSysString == "cylindrical")
    {
      RealVect stretch(RealVect::Unit);
      if (coordSysPP.contains("stretch"))
        {
          std::vector<Real> tempvect(SpaceDim, 1);
          coordSysPP.getarr("stretch", tempvect, 0, SpaceDim);
          for (int dir=0; dir<SpaceDim; dir++) stretch[dir] = tempvect[dir];
        }
      coordSysFact= new RThetaZCSFactory(stretch);
    }
  else if (coordSysString == "spherical")
    {
      Real stretch = 1.;
      coordSysPP.query("stretch", stretch);
      Real rMin = 1.;
      coordSysPP.query("min_radius", rMin);
      Real phiMin = Pi/4.;
      coordSysPP.query("min_inclination", phiMin);
      coordSysFact = new RThetaPhiCSFactory(stretch, rMin, phiMin);
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
      //CH_assert(2*radius<=domainLength[domainLength.maxDir(false)]);
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
      RealVect origin(RealVect::Zero);
      if (coordSysPP.contains("origin"))
        {
          std::vector<Real> tempvect(SpaceDim, 0);
          coordSysPP.getarr("origin", tempvect, 0, SpaceDim);
          for (int dir=0; dir<SpaceDim; dir++) origin[dir] = tempvect[dir];
        }
      RealVect stretch(RealVect::Unit);
      if (coordSysPP.contains("stretch"))
        {
          std::vector<Real> tempvect(SpaceDim, 1);
          coordSysPP.getarr("stretch", tempvect, 0, SpaceDim);
          for (int dir=0; dir<SpaceDim; dir++) stretch[dir] = tempvect[dir];
        }

      coordSysFact = new CartesianCSFactory(origin, stretch);
    }

//--Create the level factory and initialize

  AMRLevelVolTestFactory* amrLevelFactory =
    new AMRLevelVolTestFactory(coordSysFact, a_icType, refType, tagBufferSize,
                               plotPrefix);
  AMRLevel::verbosity(verbosity);

//--Create the AMR object and initialize

  AMR amr;
  amr.define(maxLevel, refRatios, probDomain, amrLevelFactory);

  amr.plotPrefix(plotPrefix);
  amr.checkpointPrefix(checkpointPrefix);
  amr.plotInterval(plotInterval);
  amr.checkpointInterval(checkpointInterval);
  amr.maxGridSize(maxGridSize);
  amr.maxBaseGridSize(maxGridSize);
  amr.dtToleranceFactor(dtToleranceFactor);
  amr.fillRatio(fillRatio);
  amr.blockFactor(blockFactor);
  amr.gridBufferSize(gridBufferSize);
  amr.verbosity(verbosity);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(maxDtGrowth);
//   amr.timeEps();
  amr.fixedDt(fixedDt);
//   amr.initialTime();
  amr.useSubcyclingInTime(useSubcycling);

/*--------------------------------------------------------------------*
 * Setup the grid and solution
 *--------------------------------------------------------------------*/

  // Set up input files
  if (!pp.contains("restart_file"))
  {
    if (!pp.contains("fixed_hierarchy"))
    {
      amr.setupForNewAMRRun();
    }
    else
    {
      MayDay::Error("Fixed hierarchy not available for test framework");
//       int numLevels = maxLevel+1;
//       Vector<Vector<Box> > amrGrids(numLevels);
//       GenFuncs::readBoxes(amrGrids, pp, probDomain,
//                           maxGridSize, blockFactor,
//                           numLevels, refRatios, verbosity,
//                           a_scaleMesh);
//       amr.setupForFixedHierarchyRun(amrGrids,1);
    }
  }
  else
  {
    std::string restartFile;
    pp.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
  }

/*--------------------------------------------------------------------*
 * Run the problem
 *--------------------------------------------------------------------*/

  // Run and time the computation
  amr.run(stopTime,nstop);

/*--------------------------------------------------------------------*
 * Conclude and cleanup
 *--------------------------------------------------------------------*/

  // Output the last plot file and statistics - time the process
  amr.conclude();

  delete amrLevelFactory;
}

/*--------------------------------------------------------------------*
 *  Parse the standard test options (-v -q) out of the command line.
 *--------------------------------------------------------------------*/

int
parseTestOptions( int argc ,char* argv[] )
{
  int i = 1;
  for (; i < argc ; ++i )
    {
      if ( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              // argv[i] = "" ;
            }
          else if ( strncmp( argv[i], "--help", 7 ) == 0)
            {
              pout() << "Usage: <executable name> [-v -q --help] "
                "[ [<inputfile>] <PP args>]" << endl;
              MayDay::Abort();
            }
        }
      else
        {
          return i;
        }
    }
  return i;
}


/*******************************************************************************
 ///////////////////////////////////////////////////////////////////////////////
 *******************************************************************************
 *
 * src/AMRLevelVolTestFactory.cpp
 *
 ******************************************************************************/

/*******************************************************************************
 *
 * \file AMRLevelVolTestFactory.cpp
 *
 * \brief Non-inline definitions for class AMRLevelVolTestFactory
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Default constructor
/*--------------------------------------------------------------------*/

AMRLevelVolTestFactory::AMRLevelVolTestFactory(
  NewCoordSysFactory *const a_coordSysFact,
  const VT::ICType          a_icType,
  const VT::RefType         a_refType,
  const int                 a_tagBufferSize,
  const std::string&        a_plotPrefix)
  :
  m_coordSysFact(new SingleBlockCSAdaptorFactory(a_coordSysFact)),
  m_icType(a_icType),
  m_refType(a_refType),
  m_tagBufferSize(a_tagBufferSize),
  m_plotPrefix(a_plotPrefix)
{
}

/*--------------------------------------------------------------------*/
///  Destructor
/**  Will delete coordinate system factory.
 *//*-----------------------------------------------------------------*/

AMRLevelVolTestFactory::~AMRLevelVolTestFactory()
{
  delete m_coordSysFact;
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Create a new AMR level
/*--------------------------------------------------------------------*/

AMRLevel*
AMRLevelVolTestFactory::new_amrlevel() const
{
  AMRLevelVolTest* amr = new AMRLevelVolTest(m_coordSysFact,
                                             m_icType,
                                             m_refType,
                                             m_tagBufferSize,
                                             m_plotPrefix);
//   amr->time();                 // default 0. or set by AMR.initialTime(Real)
//   amr->initialDtMultiplier();  // default 0.1
  return static_cast<AMRLevel*>(amr);
}


/*******************************************************************************
 ///////////////////////////////////////////////////////////////////////////////
 *******************************************************************************
 *
 * src/AMRLevelVolTest.cpp
 *
 ******************************************************************************/

/*******************************************************************************
 *
 * \file AMRLevelVolTest.cpp
 *
 * \brief Non-inline definitions for class AMRLevelVolTest
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Default constructor
/*--------------------------------------------------------------------*/

AMRLevelVolTest::AMRLevelVolTest()
  :
  m_numStates(1),
  m_levelStep(0),
  m_coordSysFact(NULL),
  m_icType(VT::ICTypeConstant),
  m_refType(VT::RefTypeRotating),
  m_levelGridMetrics(m_numStates, 4)
{
}

/*--------------------------------------------------------------------*/
///  Factory constructor
/*--------------------------------------------------------------------*/

AMRLevelVolTest::AMRLevelVolTest(
  MultiBlockCoordSysFactory *const a_coordSysFact,
  const VT::ICType          a_icType,
  const VT::RefType         a_refType,
  const int                 a_tagBufferSize,
  const std::string&        a_plotPrefix)
  :
  m_numStates(1),
  m_levelStep(0),
  m_coordSysFact(a_coordSysFact),
  m_icType(a_icType),
  m_refType(a_refType),
  m_tagBufferSize(a_tagBufferSize),
  m_plotPrefix(a_plotPrefix),
  m_levelGridMetrics(m_numStates, 4)
{
  // verbosity(3);
}

/*--------------------------------------------------------------------*/
///  Destructor
/*--------------------------------------------------------------------*/

AMRLevelVolTest::~AMRLevelVolTest()
{
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::define(AMRLevel*            a_coarser_level_ptr,
                        const ProblemDomain& a_problem_domain,
                        int                  a_level,
                        int                  a_ref_ratio)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::define " << a_level << endl;
    }

  // Call inherited define
  AMRLevel::define(a_coarser_level_ptr,
                   a_problem_domain,
                   a_level,
                   a_ref_ratio);

  // Miscellaneous terms
  m_stateNames.push_back("u");
  m_numGhost = 1;  // (switch to m_numGhostJU) sufficient for second-order gradients of <JU>
  m_ghostVect = m_numGhost*IntVect::Unit;
//  m_numGhostU = 5;  // Standard for 4th-order Cartesian as of 11/29/10
  m_numGhostU = 1;
  m_ghostVectU = m_numGhostU*IntVect::Unit;

  // Define the metrics
  LevelGridMetrics* coarserLGM = NULL;
  if (a_coarser_level_ptr)
    coarserLGM = &(static_cast<AMRLevelVolTest*>(a_coarser_level_ptr)->
                   m_levelGridMetrics);
  RealVect dxVect;  // Computational mesh spacing
  int maxSize = 0;
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      maxSize = std::max(maxSize, a_problem_domain.domainBox().size(dir));
    }
  const Real minDx = 1./maxSize;
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      dxVect[dir] = minDx;
    }

  if (s_verbosity >= 3)
    {
      pout() << "DX: " << dxVect << std::endl;
    }
  m_levelGridMetrics.define(this,
                            m_coordSysFact,
                            coarserLGM,
                            dxVect,
                            m_ghostVectU);
}

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

Real
AMRLevelVolTest::advance()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::advance " << m_level << " step " << m_dt
             << endl;
    }

//**
//   if (m_level == 0)
//     {  // DBG
//       Real norms[3];
//       computeNorm(0, norms);
//       const int savedPrec = pout().precision();
//       pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
//       pout().precision(std::numeric_limits<Real>::digits10);
//       pout() << "Conservation at begining of advance for step " << m_levelStep
//              << " on level " << m_level << ": " << norms[0] << endl;
//       pout().setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
//       pout().precision(savedPrec);
//     }
//**

//--Fill ghost cells

  fillGhostsAndComputeU();

  m_time += m_dt;
  ++m_levelStep;
  return m_dt;
}

/*--------------------------------------------------------------------*/
///  Post time-step
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::postTimeStep()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::postTimeStep " << m_level << endl;
    }

  if (m_hasFiner)
    {
      // Average from finer level data
      AMRLevelVolTest *volTestFinerPtr = getFinerLevel();
      volTestFinerPtr->m_coarseAverage.averageToCoarse(
        m_JU,
        volTestFinerPtr->m_JU);
      // Average metric terms (these only average if required)
      m_levelGridMetrics.postTimeStep();
    }
}

/*--------------------------------------------------------------------*/
///  Tag for regrid
/*--------------------------------------------------------------------*/

void AMRLevelVolTest::tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::tagCells " << m_level << " (step: "
             << m_levelStep << ")\n";
    }

  IntVectSet localTags;

  switch (m_refType)
    {
    case VT::RefTypeCenter:
    {
      const int iTime = static_cast<int>(time()/0.25);
      if (((iTime % 2) == 0))
        {
          const DisjointBoxLayout& grid = m_levelGridMetrics.getBoxes();
          DataIterator dit = grid.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const Box box = grid[dit];
              const Box refBox = grow(box, -3);
              BoxIterator bit(box);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  if (refBox.contains(bit()))
                    {
                      localTags |= bit();
                    }
                }
            }
        }
      break;
    }
    case VT::RefTypeRotating:
    {

// First find the new location of a point originating at 'locA' and rotating
// about 'axis'.  The new location is given by 'locB'.  Assumes a Cartesian
// coordinate system with origin at 0 and a domain from -1 to 1

      RealVect locB;   // Rotated position vector
      {
        RealVect locA;   // Original position vector (0.35355, 0.35355, 0.35355)
        D_EXPR(
          locA[0] = 0.4*sqrt(0.4),
          locA[1] = 0.4*sqrt(0.4),
          locA[2] = 0.4*sqrt(0.4)
          );
        RealVect RM[3];  // Rotation matrix

        // The rotation matrix is a function of angle 'theta' rotation about
        // unit vector 'axis' (right-hand rule)
        Real axis[3] =
        {
          0, 0, 1
        };
        // One of 3 axis that rotates to -loc at 180 deg (assuming
        // loc[0] = locA[1] = locA[2]).  The others are (0, 1, -1) and
        // (-1, 0, 1).
        if (SpaceDim == 3)
        {
          axis[0] =  1./std::sqrt(2.);
          axis[1] = -1./std::sqrt(2.);
          axis[2] =  0./std::sqrt(2.);
        }
        Real theta = 2*PI*time();
        Real c = std::cos(theta);
        Real s = std::sin(theta);
        Real c1 = 1 - c;

        D_EXPR(
          RM[0][0] = c1*std::pow(axis[0], 2) + c,
          RM[0][1] = c1*axis[0]*axis[1] - s*axis[2],
          RM[0][2] = c1*axis[0]*axis[2] + s*axis[1]
          );

        D_EXPR(
          RM[1][0] = c1*axis[0]*axis[1] + s*axis[2],
          RM[1][1] = c1*std::pow(axis[1], 2) + c,
          RM[1][2] = c1*axis[1]*axis[2] - s*axis[0]
          );

        D_EXPR(
          RM[2][0] = c1*axis[0]*axis[2] - s*axis[1],
          RM[2][1] = c1*axis[1]*axis[2] + s*axis[0],
          RM[2][2] = c1*std::pow(axis[2], 2) + c
          );

        // Compute the new location
        D_EXPR(
          locB[0] = RM[0].dotProduct(locA),
          locB[1] = RM[1].dotProduct(locA),
          locB[2] = RM[2].dotProduct(locA)
          );
      }

// Offset B into a domain from 0 to 1

      D_EXPR(
        locB[0] = 0.5*locB[0] + 0.5,
        locB[1] = 0.5*locB[1] + 0.5,
        locB[2] = 0.5*locB[2] + 0.5
        );

// Assume a circle/sphere with radius 'rad' and tag cells within the radius

      const Real rad = 0.1;
//  const Real rad = 0.6;  // To tag the upper quadrant for really coarse meshes

      // Map indices to a hypothetical domain extending from 0 to 1
      RealVect dx2(D_DECL(1./problemDomain().size(0),
                          1./problemDomain().size(1),
                          1./problemDomain().size(2)));

      const DisjointBoxLayout& grid = m_levelGridMetrics.getBoxes();
      DataIterator dit = grid.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
//       NewCoordSys *const coordSys = m_levelGridMetrics.getCoordSys(dit);
          const Box box = grid[dit];
          BoxIterator bit(box);
          for (bit.begin(); bit.ok(); ++bit)
            {
              RealVect Xi = bit();
              Xi += 0.5;
              Xi *= dx2;
//           Xi *= m_levelGridMetrics.dxVect();
//           RealVect X = coordSys->realCoord(Xi);
              RealVect X = Xi;
              X -= locB;
              if (X.vectorLength() <= rad)
                {
                  localTags |= bit();
                }
            }
        }
      break;
    }
    }  // End of switch

  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  if (m_icType == VT::ICTypeLinear)
    {
      const Box testDomainBox = grow(m_problem_domain.domainBox(), -1);
      if (!testDomainBox.contains(localTagsBox))
        {
          pout() << "Warning: AMRLevelVolTest::tagCells has tagged cells "
            "adjacent to the domain boundary.  Conservation tests will likely "
            "fail.\n";
        }
    }
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;

  a_tags = localTags;
}

/*--------------------------------------------------------------------*/
///  Tag cells at the start of a solution
/**  The grid is not refined at the start
 *//*-----------------------------------------------------------------*/

void
AMRLevelVolTest::tagCellsInit(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::tagCellsInit " << m_level << endl;
    }

  if (m_refType == VT::RefTypeRotating)
    {
      tagCells(a_tags);
    }
}

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::preRegrid(int                         a_base_level,
                           const Vector<Vector<Box> >& a_new_grids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::preRegrid " << a_base_level << "  "
             << m_level << endl;
    }

  // Note that U needs to be up-to-date here.  This requires doing what was done
  // in advance() again if the solution actually changed.  U only needs to be
  // maintained in the valid region.  U needs at least one
  // layer of ghosts in the valid region (do not worry about invalid ghosts).
  // Be careful because this routine is going from fine to coarse, but
  // fortunately, we don't need any ghost cells from U in the coarser level
  // Here we can just get away with an exchange.
  //**FIXME This only needs to exchange 1 layer of ghost cells so we can correct
  //**the valid cells.
  m_U.exchange(m_exchangeCopier1);

//   std::cout << "The new grid: " << a_new_grids[m_level][0] << std::endl; DBG
  DisjointBoxLayout grids = loadBalance(a_new_grids[m_level]);
  const LevelData<FArrayBox>* coarseUPtr = NULL;
  if (m_hasCoarser)
    {
      coarseUPtr = &getCoarserLevel()->m_U;
    }
  m_levelGridMetrics.preRegrid(a_base_level,
                               grids,
                               coarseUPtr,
                               m_U,
                               m_JU);

  // We need proper ghosts in m_JU now for interpolating to fine
  //**I don't think this is necessary since during regrid, we copy to a
  //**coarsened fine mesh which is only filled by valid cells.  If properly
  //**nested, this is all the information we need.
//   m_JU.exchange(m_exchangeCopier);
}

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

void AMRLevelVolTest::regrid(const Vector<Box>& a_new_grids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::regrid " << m_level << endl;
    }

  m_level_grids = a_new_grids;

  // Save the old data
  LevelData<FArrayBox> JUOld;
  JUOld.define(m_JU);

  // Note:  if grids are stored in AMRLevel*, they will need to be defined
  // for any fine levels that no longer exist.  This is because preRegrid, where
  // mapped grids are normally defined, is not called for these levels.
  //  However, regrid is called for these levels.  LevelGridMetrics, on the
  // other hand, properly allocates these with 0 boxes during the preRegrid of
  // the "new" finest level.  I.e., do this (since LGM has marked the metrics as
  // undefined on these levels):
//   if (!m_levelGridMetrics.metricsDefined()))
//     {
//       m_grids = loadBalance(a_newGrids);
//     }

  // Set up data structures
  levelSetup();

  // Reshape with the new grids
  m_U.define(m_levelGridMetrics.getBoxes(), m_numStates, m_ghostVectU);
  m_JU.define(m_levelGridMetrics.getBoxes(), m_numStates, m_ghostVect);

//--Fill <JU>

  // Interpolate from coarser level
  if (m_hasCoarser)
    {
      Interval dummyIntv;  // Only used for multiblock
      m_levelGridMetrics.regrid(m_JU,
                                getCoarserLevel()->m_U,
                                getCoarserLevel()->m_JU,
                                dummyIntv);
    }

  // Copy from old state
  JUOld.copyTo(m_JU.interval(),
               m_JU,
               m_JU.interval());

//--We need <U> at the interior if there is a finer level (use the
//--LevelGridMetrics flag because it knows if a level is used or not).  To
//--compute <U>, we need to first exchange <JU>.

  if (m_levelGridMetrics.hasFiner())
    {
      m_JU.exchange(m_exchangeCopier1);
      // Compute the new <U>
      m_levelGridMetrics.computeValidU(m_U, m_JU);
    }
}

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

void AMRLevelVolTest::postRegrid(int a_base_level)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::postRegrid " << m_level << endl;
    }
}

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

void AMRLevelVolTest::initialGrid(const Vector<Box>& a_new_grids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::initialGrid " << m_level << endl;
    }

  // Save grids
  m_level_grids = a_new_grids;

  // Load balance and store the DBL in metrics class
  DisjointBoxLayout grids = loadBalance(a_new_grids);
  m_levelGridMetrics.initialGrid(&grids);
  levelSetup();

  m_U.define(m_levelGridMetrics.getBoxes(), m_numStates, m_ghostVectU);
  m_JU.define(m_levelGridMetrics.getBoxes(), m_numStates, m_ghostVect);
}

/*--------------------------------------------------------------------*/
///  Compute the grid metrics.
/*--------------------------------------------------------------------*/

void  AMRLevelVolTest::postInitialGrid(const bool a_restart)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::postInitialGrid " << m_level << endl;
    }
  m_levelGridMetrics.postInitialGrid(NULL);
}

/*--------------------------------------------------------------------*/
///
/**
 *  \note
 *  <ul>
 *    <li> However <U> is determined, it is needed on
 *         \f$\mathcal(G)(\Omega^\ell_k,1)\f$ so use exchange if
 *         needed
 *    <li> In this example, we specify \<U\> and compute \<JU\> from
 *         \<J\> and \<U\>.  It is common to instead start with
 *         pointwise U in which case \<U\> must be computed.
 *  </ul>
 *//*-----------------------------------------------------------------*/

void AMRLevelVolTest::initialData()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::initialData " << m_level << endl;
    }
  if (m_levelGridMetrics.getBoxes().size() == 0)
    {
      return;
    }

//--Specify <U>

  switch (m_icType)
    {
    case VT::ICTypeConstant:
    {
      for (DataIterator dit = m_levelGridMetrics.getDataIterator(); dit.ok();
           ++dit)
        {
          m_U[dit].setVal(1.0);
        }
      break;
    }
    case VT::ICTypeLinear:
    {
      // Specify linear solution from 1 to 2
      // Map indices to a hypothetical domain extending from 0 to 1
      RealVect dx2(D_DECL(1./problemDomain().size(0),
                          1./problemDomain().size(1),
                          1./problemDomain().size(2)));
      for (DataIterator dit = m_levelGridMetrics.getDataIterator(); dit.ok();
           ++dit)
        {
          FArrayBox &fabU = m_U[dit];
          Box initBox = grow(m_levelGridMetrics.getBoxes()[dit], 1);
          for (BoxIterator bit(initBox); bit.ok(); ++bit)
            {
              RealVect Xi = bit();
              Xi += 0.5;
              Xi *= dx2;
//           RealVect X = coordSys->realCoord(Xi);
              RealVect X = Xi;
              fabU(bit()) = X.sum() + 1.;
            }
        }
      break;
    }
    }

//--Compute <JU> from <J> and <U>

  for (DataIterator dit = m_levelGridMetrics.getDataIterator(); dit.ok();
       ++dit)
    {
      const Box interiorBox = m_levelGridMetrics.getBoxes()[dit];
      FArrayBox& JUFab = m_JU[dit];
      const FArrayBox& UFab = m_U[dit];
      const FArrayBox& JFab = m_levelGridMetrics.m_J[dit];
      fourthOrderCellProd(JUFab, UFab, JFab, interiorBox, problemDomain());
    }
}

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

void AMRLevelVolTest::postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::postInitialize " << m_level << endl;
    }

  if (m_hasFiner)
    {
      // Volume weighted average from finer level data
      AMRLevelVolTest* finerLevel = getFinerLevel();
      finerLevel->m_coarseAverage.averageToCoarse(m_JU, finerLevel->m_JU);
    }

//--This is our conserved quantity.  Find L1 norm for testing conservation

  if (m_level == 0 && m_icType == VT::ICTypeLinear)
    {
      Real norms[3];
      computeNorm(0, norms);
      VT::g_conservedJU0 = norms[0];
    }
}

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

Real AMRLevelVolTest::computeDt()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::computeDt " << m_level << endl;
    }

  return m_dt;
}

/*--------------------------------------------------------------------*/
///
/*--------------------------------------------------------------------*/

Real AMRLevelVolTest::computeInitialDt()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::computeInitialDt " << m_level << endl;
    }

  return m_dt;
}

//--HDF5 Routines

#ifdef CH_USE_HDF5

/*--------------------------------------------------------------------*/
///  Write the checkpoint header
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::writeCheckpointHeader" << endl;
    }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = m_stateNames[comp];
    }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }
}

/*--------------------------------------------------------------------*/
///  Write the checkpoint data
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::writeCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_int ["tag_buffer_size"] = m_tagBufferSize;
  header.m_real["dx"]              = m_levelGridMetrics.dxVect()[0];  //***FIXME
  header.m_real["dt"]              = m_dt;
  header.m_real["time"]            = m_time;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM(if (m_problem_domain.isPeriodic(0))
           {
             header.m_int ["is_periodic_0"] = 1;
           }
         else
           {
             header.m_int ["is_periodic_0"] = 0;
           } ,

         if (m_problem_domain.isPeriodic(1))
           {
             header.m_int ["is_periodic_1"] = 1;
           }
         else
           {
             header.m_int ["is_periodic_1"] = 0;
           } ,

         if (m_problem_domain.isPeriodic(2))
           {
             header.m_int ["is_periodic_2"] = 1;
           }
         else
           {
             header.m_int ["is_periodic_2"] = 0;
           } );

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  write(a_handle,m_JU.boxLayout());
  write(a_handle,m_JU,"data");
}

/*--------------------------------------------------------------------*/
///  Read checkpoint header
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::readCheckpointHeader(HDF5Handle& a_handle)
{

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::readCheckpointHeader" << endl;
    }

  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << "hdf5 header data:" << endl;
      pout() << header << endl;
    }

  // Get the number of components
  if (header.m_int.find("num_components") == header.m_int.end())
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointHeader: checkpoint file "
                    "does not have num_components");
    }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates)
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointHeader: num_components in "
                    "checkpoint file does not match solver");
    }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      if (header.m_string.find(compStr) == header.m_string.end())
        {
          MayDay::Error("AMRLevelVolTest::readCheckpointHeader: checkpoint "
                        "file does not have enough component names");
        }

      stateName = header.m_string[compStr];
      if (stateName != m_stateNames[comp])
        {
          MayDay::Error("AMRLevelVolTest::readCheckpointHeader: state_name in "
                        "checkpoint does not match solver");
        }
    }
}

/*--------------------------------------------------------------------*/
///  Read checkpoint data for this level
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << "hdf5 header data:" << endl;
      pout() << header << endl;
    }

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointLevel: file does not "
                    "contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
    {
      pout() << "read ref_ratio = " << m_ref_ratio << endl;
    }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointLevel: file does not "
                    "contain tag_buffer_size");
    }
  m_tagBufferSize = header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
    {
      pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
    }

//   // Get dx
//   if (header.m_real.find("dx") == header.m_real.end())
//     {
//       MayDay::Error("AMRLevelVolTest::readCheckpointLevel: file does not "
//                     "contain dx");
//     }
//   m_dx = header.m_real["dx"];

//   if (s_verbosity >= 2)
//     {
//       pout() << "read dx = " << m_dx << endl;
//     }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointLevel: file does not "
                    "contain dt");
    }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
    {
      pout() << "read dt = " << m_dt << endl;
    }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointLevel: file does not "
                    "contain time");
    }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
    {
      pout() << "read time = " << m_time << endl;
    }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointLevel: file does not "
                    "contain prob_domain");
    }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
           {
             isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
           }
         else
           {
             isPeriodic[0] = false;
           } ,

         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
           {
             isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
           }
         else
           {
             isPeriodic[1] = false;
           } ,

         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
           {
             isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
           }
         else
           {
             isPeriodic[2] = false;
           } );

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> boxes;
  const int gridStatus = read(a_handle,boxes);

  if (gridStatus != 0)
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointLevel: file does not "
                    "contain a Vector<Box>");
    }

  // Create level domain
  DisjointBoxLayout grids = loadBalance(boxes);

  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = grids;

  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = constGrids[lit()];
      m_level_grids.push_back(b);
    }

  if (s_verbosity >= 4)
    {
      pout() << "read level domain: " << endl;
      LayoutIterator lit = grids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& b = grids[lit()];
          pout() << lit().intCode() << ": " << b << endl;
        }
      pout() << endl;
    }

  // Reshape state with new grids
  m_JU.define(grids,m_numStates);
  const int dataStatus = read<FArrayBox>(a_handle,
                                         m_JU,
                                         "data",
                                         grids);

  if (dataStatus != 0)
    {
      MayDay::Error("AMRLevelVolTest::readCheckpointLevel: file does not "
                    "contain state data");
    }

  // Set up data structures
  levelSetup();
  //**FIXME need to define m_U
}

/*--------------------------------------------------------------------*/
///  Write plotfile header
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::writePlotHeader" << endl;
    }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates + 2;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = m_stateNames[comp];
    }
  header.m_string["component_1"] = "Ju";
  header.m_string["component_2"] = "J";

  // Write the header
  header.writeToFile(a_handle);
  a_handle.setGroup("/Expressions");
//   HDF5HeaderData expressions;
//   m_level.getPhysicsPtrConst()->expressions(expressions);
//   expressions.writeToFile(a_handle);

  if (s_verbosity >= 4)
    {
      pout() << header << endl;
    }
}

/*--------------------------------------------------------------------*/
///  Write plotfile data for this level
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::writePlotLevel " << m_level << endl;
    }

  // Write out the mapped-grid geometry info (only effective at level 0)
  writeMappedPlotFile();

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_levelGridMetrics.dxVect()[0];  // ***FIXME
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 4)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  LevelData<FArrayBox> outData(m_levelGridMetrics.getBoxes(), 3, IntVect::Unit);
  m_U.copyTo(m_U.interval(), outData, Interval(0, 0));
  m_JU.copyTo(m_JU.interval(), outData, Interval(1, 1));
  m_levelGridMetrics.m_J.copyTo(m_levelGridMetrics.m_J.interval(),
                                outData,
                                Interval(2, 2));

  write(a_handle, m_levelGridMetrics.getBoxes());
  write(a_handle, outData, "data", IntVect::Unit);
}

/*--------------------------------------------------------------------*/
///  Write the mapping for all levels
/*   This is really annoying.  Can't we put this as new blocks in the
 *   HDF5 file?
 *--------------------------------------------------------------------*/

void AMRLevelVolTest::writeMappedPlotFile() const
{
  // only do this on level 0
  if (m_level == 0)
    {
      if (s_verbosity >= 3)
        {
          pout() << "AMRLevelVolTest::writeMappedPlotFile" << endl;
        }

      // gather AMR Levels and create node-centered dataset of
      // node locations of mapped grids

      Vector<AMRLevel*> vectAMRLevels;
      {
        // cast away const for this to call this function
        AMRLevelVolTest* nonConstThis = const_cast<AMRLevelVolTest*>(this);

        vectAMRLevels = nonConstThis->getAMRLevelHierarchy();
      }

      int numLevels = vectAMRLevels.size();

      Vector<int> vectRefRatio(numLevels,0);
      Vector<DisjointBoxLayout> vectGrids(numLevels);
      Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(numLevels, NULL);

      const AMRLevelVolTest* levelPtr = this;
      // loop over levels and set things up
      for (int lev=0; lev<numLevels; lev++)
        {
          const DisjointBoxLayout& levelGrids =
            levelPtr->m_levelGridMetrics.getBoxes();
          vectGrids[lev] = levelGrids;
          vectRefRatio[lev] = levelPtr->m_ref_ratio;
          vectNodeLoc[lev] = new LevelData<NodeFArrayBox>(levelGrids,
                                                          SpaceDim,
                                                          IntVect::Unit);

          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const NewCoordSys *const levelCS =
                levelPtr->m_levelGridMetrics.getCoordSys(levelGrids[dit]);
              NodeFArrayBox& thisNodeFAB = (*vectNodeLoc[lev])[dit];
              // node-centered FAB
              FArrayBox& thisFAB = thisNodeFAB.getFab();
              BoxIterator bit(thisFAB.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  // location in index space
                  RealVect nodeIndexLoc(iv);
                  nodeIndexLoc *= levelPtr->m_levelGridMetrics.dxVect();
                  // now convert to real space location
                  RealVect nodeRealLoc = levelCS->realCoord(nodeIndexLoc);
                  D_EXPR6(thisFAB(iv,0) = nodeRealLoc[0],
                          thisFAB(iv,1) = nodeRealLoc[1],
                          thisFAB(iv,2) = nodeRealLoc[2],
                          thisFAB(iv,3) = nodeRealLoc[3],
                          thisFAB(iv,4) = nodeRealLoc[4],
                          thisFAB(iv,5) = nodeRealLoc[5]);
                }
            } // end loop over grids

          // advance to next level
          levelPtr = levelPtr->getFinerLevel();
        } // end loop over levels

      // create names
      Vector<string> locationNames(SpaceDim);
      D_EXPR6(locationNames[0] = "x",
              locationNames[1] = "y",
              locationNames[2] = "z",
              locationNames[3] = "u",
              locationNames[4] = "v",
              locationNames[5] = "w");

      // create filename
      char iter_str[80];

      sprintf(iter_str,
              "%s%06d.%dd.map.hdf5",
              m_plotPrefix.c_str(), m_levelStep, SpaceDim);

      string fileName(iter_str);

      // Choose a dx
      Real dx = m_levelGridMetrics.dxVect()[
        m_levelGridMetrics.dxVect().maxDir(false)];
      // now call nodal WriteAMRHierarchy function...
      WriteAMRHierarchyHDF5(fileName,
                            vectGrids,
                            vectNodeLoc,
                            locationNames,
                            m_problem_domain.domainBox(),
                            dx,
                            m_dt,
                            m_time,
                            vectRefRatio,
                            numLevels);

      // now clean up our mess
      for (int lev=0; lev<vectNodeLoc.size(); lev++)
        {
          delete vectNodeLoc[lev];
        }
    }
}
#endif

/*--------------------------------------------------------------------*/
///  Conclude by computing error norms
/*
 *--------------------------------------------------------------------*/

void
AMRLevelVolTest::conclude(int a_step) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::conclude " << m_level  << endl;
    }
  if (m_level == 0)
    {
      switch (m_icType)
        {
        case VT::ICTypeConstant:
        {
          // Test that <U> is still 1.0 everywhere
          VT::g_runStatus = allCellsTest(1.0);
          break;
        }
        case VT::ICTypeLinear:
        {
          Real norms[3];
          // Test that <JU> was conserved
          computeNorm(0, norms);
          VT::g_conservedJU1 = norms[0];
          break;
        }
        }
    }
}

/*--------------------------------------------------------------------*/
///  Create a load-balanced DisjointBoxLayout from a collection of
///  boxes
/**  \param[in]  a_grids
 *                      Collection of boxes
 *   \return            The load balanced layout
 *//*-----------------------------------------------------------------*/

DisjointBoxLayout
AMRLevelVolTest::loadBalance(const Vector<Box>& a_grids) const
{
  // Load balance and create boxlayout
  Vector<int> procMap;

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
    {
      pout() << "AMRLevelAdvect::loadBalance: procesor map: " << endl;

      for (int igrid = 0; igrid < a_grids.size(); ++igrid)
        {
          pout() << igrid << ": " << procMap[igrid] << "  " << endl;
        }

      pout() << endl;
    }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}

/*--------------------------------------------------------------------*/
///  Set up data structures for the level
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::levelSetup()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelVolTest::levelSetup " << m_level << endl;
    }

  m_exchangeCopier1.exchangeDefine(m_levelGridMetrics.getBoxes(),
                                   IntVect::Unit);
  m_exchangeCopier.exchangeDefine(m_levelGridMetrics.getBoxes(),
                                  m_ghostVectU);

  AMRLevelVolTest* volTestCoarserPtr = getCoarserLevel();
  AMRLevelVolTest* volTestFinerPtr   = getFinerLevel();

  m_hasCoarser = (volTestCoarserPtr != NULL);
  m_hasFiner   = (volTestFinerPtr   != NULL);

  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();

      m_coarseAverage.define(m_levelGridMetrics.getBoxes(),
                             m_numStates,
                             nRefCrse);
    }
}

/*--------------------------------------------------------------------*/
///  Fill ghost cells
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::fillGhostsAndComputeU()
{

//--From neighbouring boxes on this level

  m_JU.exchange(m_exchangeCopier1);

//--From a coarser level (also sets <U> in the invalid ghost cells)

  if (m_hasCoarser)
    {
      Interval dummyIntv;  // Only used for multiblock
      getCoarserLevel()->m_levelGridMetrics.presetThis2FnInterpolatorCrFnLevU(
        getCoarserLevel()->m_U, dummyIntv);
      m_levelGridMetrics.fillFineGhostCells(m_U, m_JU);
    }

//--Compute the new <U> at the interior.

  m_levelGridMetrics.computeValidU(m_U, m_JU);

//--A full exchange of <U> to set the valid ghost cells

  m_U.exchange(m_exchangeCopier);
}

/*--------------------------------------------------------------------*/
///  Get the next coarser level
/*--------------------------------------------------------------------*/

AMRLevelVolTest*
AMRLevelVolTest::getCoarserLevel() const
{
  AMRLevelVolTest* volTestCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
    {
      volTestCoarserPtr = dynamic_cast<AMRLevelVolTest*>(m_coarser_level_ptr);
      if (volTestCoarserPtr == NULL)
        {
          MayDay::Error("AMRLevelVolTest::getCoarserLevel: dynamic cast "
                        "failed");
        }
    }

  return volTestCoarserPtr;
}

/*--------------------------------------------------------------------*/
///  Get the next finer level
/*--------------------------------------------------------------------*/

AMRLevelVolTest*
AMRLevelVolTest::getFinerLevel() const
{
  AMRLevelVolTest* volTestFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
    {
      volTestFinerPtr = dynamic_cast<AMRLevelVolTest*>(m_finer_level_ptr);
      if (volTestFinerPtr == NULL)
        {
          MayDay::Error("AMRLevelVolTest::getFinerLevel: dynamic cast failed");
        }
    }

  return volTestFinerPtr;
}

/*--------------------------------------------------------------------*/
///  Compute a norm
/*--------------------------------------------------------------------*/

void
AMRLevelVolTest::computeNorm(const int a_type, Real a_norm[3]) const
{
  CH_assert(m_level == 0);

  // Gather AMR Levels
  Vector<AMRLevel*> vAMRLevels;
  {
    // cast away const for this to call this function
    AMRLevelVolTest* nonConstThis = const_cast<AMRLevelVolTest*>(this);
    vAMRLevels = nonConstThis->getAMRLevelHierarchy();
  }
  int numLevels = vAMRLevels.size();
  // Gather ref ratios
  Vector<int> vRefRatios(numLevels - 1);
  for (int iLev = 0; iLev < numLevels - 1; ++iLev)
    {
      vRefRatios[iLev] =
        dynamic_cast<AMRLevelVolTest*>(vAMRLevels[iLev])->refRatio();
    }

  switch (a_type)
    {

//--L1 norm of <JU> for testing conservation

    case 0:


    {
      Vector<LevelData<FArrayBox>*> vPhi(numLevels);
      for (int iLev = 0; iLev != numLevels; ++iLev)
        {
          vPhi[iLev] =
            &(dynamic_cast<AMRLevelVolTest*>(vAMRLevels[iLev])->m_JU);
        }
      a_norm[0] = ::computeNorm(vPhi,
                                vRefRatios,
                                1.,
                                Interval(0, 0),
                                1,
                                0);
      break;
    }

    case 1:

//--L1, L2, and LInfinity norm of <U> - <Uexact> for convergence

    {
      MayDay::Error("No error testing is allowed!");
      break;
    }
    }
}

/*--------------------------------------------------------------------*/
///  Test the value in all cells
/*--------------------------------------------------------------------*/

int
AMRLevelVolTest::allCellsTest(Real a_testVal) const
{
#ifdef CH_USE_DOUBLE
  const int tol = std::numeric_limits<Real>::digits10 - 2;
#else
  const int tol = std::numeric_limits<Real>::digits10 - 1;
#endif

  CH_assert(m_level == 0);

  // Gather AMR Levels
  Vector<AMRLevel*> vAMRLevels;
  {
    // cast away const for this to call this function
    AMRLevelVolTest* nonConstThis = const_cast<AMRLevelVolTest*>(this);
    vAMRLevels = nonConstThis->getAMRLevelHierarchy();
  }
  int numLevels = vAMRLevels.size();

  pout() << "      " << "Testing <U> in all cells:\n";
  int status = 0;
  for (int iLev = 0; iLev != numLevels; ++iLev)
    {
      const AMRLevelVolTest* volTestLevelPtr =
        dynamic_cast<const AMRLevelVolTest*>(vAMRLevels[iLev]);
      for (DataIterator dit =
             volTestLevelPtr->m_levelGridMetrics.getDataIterator();
           dit.ok(); ++dit)
        {
          const Box box = volTestLevelPtr->m_levelGridMetrics.getBoxes()[dit];
          const FArrayBox& fab = volTestLevelPtr->m_U[dit];
          for (BoxIterator bit(box); bit.ok(); ++bit)
            {
              if (GenFuncs::compare(fab(bit()), a_testVal, tol))
                {
                  pout() << "      " << "   Failure (" << fab(bit())
                         << ") in cell " << bit() << " on level " << iLev
                         << std::endl;
                  ++status;
                }
            }
        }
    }
  if (status == 0)
    {
      pout() << "      " << "   Success\n";
    }
  return status;
}
