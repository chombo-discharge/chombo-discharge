#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>
#include "parstream.H"

#include "LayoutIterator.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "CH_Timer.H"

#include "AMRLevelSelfGravity.H"
#include "AMRIO.H"
#include "Ancillae.H"
#include "PhysIBC.H"
#include "FABView.H"

#include "AMRMultiGrid.H"
#include "AMRPoissonOp.H"
#include "BiCGStabSolver.H"
#include "RelaxSolver.H"
#include "BCFunc.H"

#include "LoHiCenter.H"
#include "LoHiSide.H"
#include "computeSum.H"

#include "InterpF_F.H"

#define SOLVER_RESET_PHI false
#define SOLVER_NORM_TYPE (1)
#define SOLVER_MAX_ITER  (20)
#define SOLVER_MIN_ITER  (2)

#ifdef CH_USE_FLOAT
#  define SOLVER_TOLERANCE (1.0e-4)
#else
#  define SOLVER_TOLERANCE (1.0e-10)
#endif

#undef  USE_RELAX_SOLVER
#define SOLVER_NUM_SMOOTH 8
#define SOLVER_NUM_MG     1
#define SOLVER_HANG       1.e-14
#define SOLVER_NORM_THRES 1.e-10

//this is a dummy function for periodic BCs
void NoOpBc(FArrayBox& a_state,
            const Box& a_box,
            const ProblemDomain& a_domain,
            Real a_dx,
            bool a_homogeneous)
{
}

// Constructor
AMRLevelSelfGravity::AMRLevelSelfGravity()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity default constructor" << endl;
  }
  m_initial_dt_multiplier = 0.1;
  m_gdnvPhysics    = NULL;
  m_refCellTagger  = NULL;
  m_gradient       = NULL;

  m_paramsDefined = false;
}

// Destructor
AMRLevelSelfGravity::~AMRLevelSelfGravity()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity destructor" << endl;
  }

  if (m_gdnvPhysics != NULL)
  {
    delete m_gdnvPhysics;
  }

  if (m_refCellTagger != NULL)
  {
    delete m_refCellTagger;
  }

  if (m_gradient != NULL)
  {
    delete m_gradient;
  }
}

// This instance should never get called - historical
void AMRLevelSelfGravity::define(AMRLevel*  a_coarserLevelPtr,
                           const Box& a_problemDomain,
                           int        a_level,
                           int        a_refRatio)
{
  ProblemDomain physdomain(a_problemDomain);

  MayDay::Error("AMRLevelSelfGravity::define -\n\tShould never be called with a Box for a problem domain");
}

// Define new AMR level
void AMRLevelSelfGravity::define(AMRLevel*            a_coarserLevelPtr,
                           const ProblemDomain& a_problemDomain,
                           int                  a_level,
                           int                  a_refRatio)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::define " << a_level << endl;
  }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelSelfGravity* amrGodPtr = dynamic_cast<AMRLevelSelfGravity*>(a_coarserLevelPtr);

    // FM:: why do we need these call if those params have been
    // already defined ? Why not all the other parameters then ?
    if (amrGodPtr != NULL)
    {
      m_cfl = amrGodPtr->m_cfl;
      m_domainLength = amrGodPtr->m_domainLength;
      m_tagBufferSize= amrGodPtr->m_tagBufferSize;
    }
    else
    {
      MayDay::Error("AMRLevelSelfGravity::define: a_coarserLevelPtr is not castable to AMRLevelSelfGravity*");
    }
  }

  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().longside();

  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more
  m_numGhost = 1;
  m_numRhsGhost = 1;
  m_numForceGhost = 2;

  // Set the number of ghost cells appropriately
  if (m_useFourthOrderSlopes) m_numForceGhost += 2;

  CH_assert(allDefined());
  // define patch integrator
  m_gdnvPhysics->define(m_problem_domain,m_dx);

  // Get additional information from the patch integrator
  m_numStates  = m_gdnvPhysics->numConserved();
  m_stateNames = m_gdnvPhysics->stateNames();
  m_plotNames  = m_gdnvPhysics->plotNames();

  // define poisson operator
  m_plotNames.push_back("potential");


  // define refCellTagger operator
  m_refCellTagger->define(m_problem_domain,m_dx,a_refRatio);
}


void AMRLevelSelfGravity::defineParams(const Real&             a_cfl,
                                       const Real&             a_domainLength,
                                       const int&              a_verbosity,
                                       const int&              a_tagBufferSize,
                                       const int&             a_maxInitRefLevel,
                                       const Real&            a_initialDtMultiplier,
                                       const GodunovPhysics* const a_godunovPhysics,
                                       const int&             a_normalPredOrder,
                                       const bool&            a_useFourthOrderSlopes,
                                       const bool&            a_usePrimLimiting,
                                       const bool&            a_useCharLimiting,
                                       const bool&            a_useFlattening,
                                       const bool&          a_useArtificialViscosity,
                                       const Real&            a_artificialViscosity,
                                       const RefCellTagger* const  a_refCellTagger,
                                       const bool&                 a_useDeltaPhiCorr,
                                       const StencilType&          a_stencil)
{
  // Set the CFL number
  m_cfl = a_cfl;

  // Set the physical dimension of the longest side of the domain
  m_domainLength = a_domainLength;

  verbosity(a_verbosity);

  // Set the tag buffer size
  m_tagBufferSize = a_tagBufferSize;

  initialDtMultiplier(a_initialDtMultiplier);

  if (m_gdnvPhysics != NULL)
  {
    delete m_gdnvPhysics;
    m_gdnvPhysics = NULL;
  }
  m_gdnvPhysics = dynamic_cast<SelfGravityPhysics*>(a_godunovPhysics->new_godunovPhysics());

  m_normalPredOrder = a_normalPredOrder;

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_usePrimLimiting      = a_usePrimLimiting;
  m_useCharLimiting      = a_useCharLimiting;
  m_useFlattening        = a_useFlattening;

  // Artificial viscosity coefficient must be greater than zero
  CH_assert(!a_useArtificialViscosity || (a_artificialViscosity >= 0.0));

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = a_useArtificialViscosity;
  m_artificialViscosity    = a_artificialViscosity;

  // create RefCellTagger object for tagging cells to be refined
  m_refCellTagger = a_refCellTagger->new_refCellTagger();

  m_maxInitRefLevel = a_maxInitRefLevel;

  // Set the stencil for differentiating the grav. potential:
  if (m_gradient != NULL)
  {
    delete m_gradient;
    m_gradient = NULL;
  }
  switch(a_stencil)
    {
    case TwoPts:
      m_gradient = new TwoPtsGradient();
      m_numPhiGhost = 1;
      break;
    case FourPts:
      m_gradient = new FourPtsGradient();
      m_numPhiGhost = 2;
      break;
    case TenPts:
      m_gradient = new TenPtsGradient();
      m_numPhiGhost = 1;
      break;
    default:
      MayDay::Error("AMRLevelSelfGravity::defineparams: a_stencil NOT defined");
    }

  // whether or not the deltaPhi correction should be applie
  m_useDeltaPhiCorr = a_useDeltaPhiCorr;

  m_paramsDefined = true;
}


// Advance by one timestep
Real AMRLevelSelfGravity::advance()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::advance level " << m_level << " at time " << m_time << endl;
  }

  Real newDt      = 1.0e9;
  Real newDt_hydr = 1.0e9;

  // Copy the new to the old
  for (DataIterator di=m_UNew.dataIterator(); di.ok(); ++di)
  {
    m_UOld[di].copy(m_UNew[di]);
  }

  // Set up arguments to LevelGodunov::step based on whether there are
  // coarser and finer levels

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  const LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  const LevelData<FArrayBox>* coarserDataOld = &dummyData;
  const LevelData<FArrayBox>* coarserDataNew = &dummyData;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelSelfGravity* coarserPtr = getCoarserLevel();

    // Recall that my flux register goes between my level and the next
    // finer level
    coarserFR = &coarserPtr->m_fluxRegister;

    coarserDataOld = &coarserPtr->m_UOld;
    coarserDataNew = &coarserPtr->m_UNew;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }

  // A finer level exists
  if (m_hasFiner)
  {
    // Recall that my flux register goes between my level and the next
    // finer level
    finerFR = &m_fluxRegister;
  }

  LevelData<FArrayBox> source;
  setupSourceTerm(source,m_UNew,m_forceNew,m_time,m_dt);

  // we don't need the flux in the simple hyperbolic case...
  LevelData<FArrayBox> flux[SpaceDim];
  newDt_hydr = m_levelGodunov.step(m_UNew,
                                   flux,
                                   *finerFR,
                                   *coarserFR,
                                   source,
                                   *coarserDataOld,
                                   tCoarserOld,
                                   *coarserDataNew,
                                   tCoarserNew,
                                   m_time,
                                   m_dt);

  // synchronize pressure and entropy
  for (DataIterator di=m_grids.dataIterator(); di.ok(); ++di)
  {
    m_gdnvPhysics->synchronize(m_UNew[di],m_UOld[di],m_grids.get(di));
  }

  // Update the time and store the new timestep
  m_time += m_dt;
  newDt  = newDt_hydr;


#ifdef GRAVITY
  newDt_hydr = hundred;
  for (DataIterator di=m_grids.dataIterator(); di.ok(); ++di)
  {
    const Box& bx = m_grids.get(di);

    Real newDt_srce =
    m_gdnvPhysics->applySource(m_UNew[di],m_UOld[di],source[di],m_time,m_dt,bx);
    newDt_hydr = Min(newDt_hydr,newDt_srce);
  }
#endif

  newDt_hydr *= m_cfl;
  if (s_verbosity >=3) pout() << " Dt/r Hydro " << newDt_hydr << endl;


#ifdef GRAVITY
  // save new Phi/Force -> old Phi/Force
  for (DataIterator di=m_phiNew.dataIterator(); di.ok(); ++di)
  {
    m_phiOld[di].copy(m_phiNew[di]);
    m_forceOld[di].copy(m_forceNew[di]);
  }

  if (m_hasFiner && !m_isThisFinestLev)
    {
      ellipticSolver(m_level,true);

      computeForce(m_forceNew,m_phiNew,m_time);
    }
#endif


  Real returnDt = globalMin(newDt);

  m_dtOld = m_dt;
  m_dtNew = returnDt;

  return returnDt;
}


// Things to do after a timestep
void AMRLevelSelfGravity::postTimeStep()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::postTimeStep " << m_level << endl;
  }

  if (m_hasFiner)
  {
    // uold, U before refluxing, needed in synchronize when Pc is used
    LevelData<FArrayBox> UOld;
    UOld.define(m_UNew);

    // Reflux
    Real scale = -1.0/m_dx;
    m_fluxRegister.reflux(m_UNew,scale);

    // synchronize pressure and entropy and more...
    for (DataIterator di=m_grids.dataIterator(); di.ok(); ++di)
    {
      m_gdnvPhysics->synchronize(m_UNew[di],UOld[di],m_grids.get(di));
    }

    // Average from finer level data
    AMRLevelSelfGravity* amrGodFinerPtr = getFinerLevel();

    amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                    amrGodFinerPtr->m_UNew);

    // in high mach number flows, preserve smooth pressure
    // distribution predating interpolation
    for (DataIterator di=m_grids.dataIterator(); di.ok(); ++di)
    {
      m_gdnvPhysics->setPressureToEntropy(m_UNew[di],m_grids.get(di));
    }
  }

#ifdef GRAVITY
  // if there is a force term, first we solve for the potential
  // (Poisson eq.) on all levels of the grid hierarchy that are
  // aligned in time. This is done by gravity() which also calls, in
  // the proper order, functions that apply 2nd order corrections and
  // compute the force field taking into account interlevel and domain
  // boundary conditions. After all this has been taken care of, if
  // necessary, particles are reassigned to the a new grid patch and
  // hierarchy level by manageParticles(). If no force term is present
  // only manageParticles() needs to be called.
  const Real eps = 5.e-2;
  Real crseTime = -one;
  if (m_level> 0) crseTime = m_coarser_level_ptr->time();
  const bool stepsLeft = abs(crseTime - m_time) > (eps*m_dt);

  if (m_level==0 || stepsLeft) gravity(m_level);
#endif

  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 2 && m_level == 0)
  {
    int nRefFine = 1;

    if (s_verbosity<3) pout() << "AMRLevelSelfGravity::postTimeStep:" << endl;
    pout() << "  Sums:" << endl;
    for (int comp = 0; comp < m_numStates; comp++)
    {
      Interval curComp(comp,comp);
      Real integral = computeSum(m_UNew,NULL,nRefFine,m_dx,curComp);

      pout() << "    " << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << integral
             << " --- " << m_stateNames[comp];

      if (comp == 0 && !first)
      {
        pout() << " (" << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << (integral-last_integral)
               << " " << setw(23)
                      << setprecision(16)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific)
                      << (integral-orig_integral)
               << ")";
      }

      pout() << endl;

      if (comp == 0)
      {
        if (first)
        {
          orig_integral = integral;
          first = false;
        }

        last_integral = integral;
      }
    }
  }

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::postTimeStep done " << endl;
  }
}


// Create tags for regridding
void AMRLevelSelfGravity::tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::tagCells " << m_level << endl;
  }

  // Create tags based on various criteria
  IntVectSet localTags;

  if (m_refCellTagger->needGhostCells())
  {
    // If there is a coarser level interpolate undefined ghost cells
    if (m_hasCoarser)
    {
      const AMRLevelSelfGravity* amrGodCoarserPtr = getCoarserLevel();
      PiecewiseLinearFillPatch pwl(m_UNew.disjointBoxLayout(),
                                 amrGodCoarserPtr->m_UNew.disjointBoxLayout(),
                                 m_numStates,
                                 amrGodCoarserPtr->m_problem_domain,
                                 amrGodCoarserPtr->m_ref_ratio,
                                 1);

      pwl.fillInterp(m_UNew,
                     amrGodCoarserPtr->m_UNew,
                     amrGodCoarserPtr->m_UNew,
                     1.0,
                     0,
                     0,
                     m_numStates);
    }
    m_UNew.exchange(Interval(0,SpaceDim+2));
  }

  if (m_refCellTagger->refineShocks())
  {
    localTags |= m_refCellTagger->tagShocks(m_UNew);
  }
  if (m_refCellTagger->refineVorticity())
  {
    localTags |= m_refCellTagger->tagVorticity(m_UNew);
  }
  if (m_refCellTagger->refineGradient("gas"))
  {
    localTags |= m_refCellTagger->tagGradient(m_UNew);
  }
  if (m_refCellTagger->refineJeans())
  {
    localTags |= m_refCellTagger->tagJeans(m_UNew);
  }

#ifdef GRAVITY
  if (m_refCellTagger->refineOverdense())
  {
    Real frhs = twothird;
    // get Poisson RHS, then convert it to mass density
    for (DataIterator di = m_rhs.dataIterator(); di.ok(); ++di)
    {
      m_rhs[di] *= frhs;
      if (m_problem_domain.isPeriodic()) m_rhs[di] += one;
    }
    localTags |= m_refCellTagger->tagOverdense(m_rhs);
  }
  if (m_refCellTagger->refineGradient("phi"))
  {
    // phi's ghost cells must have be filled when the force was computed
    localTags |= m_refCellTagger->tagGradient(m_phiNew);
  }
#endif

  // localTags may or may not include the previous tags depending on how the
  // refine-mode has been defined.
  if (m_refCellTagger->refineRegion())
  {
    localTags = m_refCellTagger->tagRegion(localTags,m_level,m_grids);
  }


  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();

#ifdef INFINITE_DOMAIN_BC
  if (m_level==0)
  {
    localTagsBox &= grow(m_problem_domain,-4);
  }
  else
#endif
  {
    localTagsBox &= m_problem_domain;
  }
  localTags &= localTagsBox;

  a_tags = localTags;
}


// Create tags at initialization
void AMRLevelSelfGravity::tagCellsInit(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::tagCellsInit " << m_level << endl;
  }
  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.

  if (m_level<m_maxInitRefLevel) tagCells(a_tags);
}

// Set up data on this level after regridding
void AMRLevelSelfGravity::regrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::regrid " << m_level << endl;
  }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  if (s_verbosity >= 4)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;

    pout() << "new grids: " << endl;
    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      pout() << constGrids[lit()] << endl;
    }
  }

  const IntVect ivGhost = m_numGhost*IntVect::Unit;
  // Save data for later
  for (DataIterator di=m_UNew.dataIterator(); di.ok(); ++di)
  {
    m_UOld[di].copy(m_UNew[di]);
  }

  // Reshape state with new grids
  m_UNew.define(m_grids,m_numStates,ivGhost);
  resetToZero(m_UNew);

#ifdef GRAVITY
  for (DataIterator di=m_phiNew.dataIterator(); di.ok(); ++di)
  {
    m_phiOld[di].copy(m_phiNew[di]);
  }
  m_phiNew.define(m_grids,1,m_numPhiGhost*IntVect::Unit);
  resetToZero(m_phiNew);

  if (m_useDeltaPhiCorr)
  {
    m_deltaPhi.define(m_grids,1,m_numPhiGhost*IntVect::Unit);
    resetToZero(m_deltaPhi);
  }

  m_forceNew.define(m_grids,SpaceDim,m_numForceGhost*IntVect::Unit);
  m_forceOld.define(m_grids,SpaceDim,m_numForceGhost*IntVect::Unit);

  m_rhs.define(m_grids,1,m_numRhsGhost*IntVect::Unit);
  resetToZero(m_rhs);
#endif

  // Set up data structures
  levelSetup();

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    AMRLevelSelfGravity* amrCrsePtr = getCoarserLevel();
    m_fineInterp.interpToFine(m_UNew,amrCrsePtr->m_UNew);

    // preserve smooth pressure distribution predating interpolation
    for (DataIterator di=m_grids.dataIterator(); di.ok(); ++di)
    {
      m_gdnvPhysics->setPressureToEntropy(m_UNew[di],m_grids.get(di));
    }

#ifdef GRAVITY
    m_fineInterpPhi.interpToFine(m_phiNew,amrCrsePtr->m_phiNew);
#endif
  }

  // Copy from old state
  m_UOld.copyTo(m_UOld.interval(),m_UNew,m_UNew.interval());
  m_UOld.define(m_grids,m_numStates,ivGhost);

#ifdef GRAVITY
  m_phiOld.copyTo(m_phiOld.interval(),m_phiNew,m_phiNew.interval());
  m_phiOld.define(m_grids,1,m_numPhiGhost*IntVect::Unit);
#endif
}


void AMRLevelSelfGravity::postRegrid(int a_base_level)
{
  if (s_verbosity >= 3)
  {
    pout() << " AMRLevelSelfGravity::postRegrid " << m_level << " base level " << a_base_level << endl;
  }

  if (m_hasFiner)
  {
    // define finest level
    AMRLevelSelfGravity* amrFinerPtr = getFinerLevel();
    m_isThisFinestLev = (amrFinerPtr->m_grids.size()<=0);
  }

  if (m_level == a_base_level+1)
  {
    // define the new finest level
    AMRLevelSelfGravity* amrCoarserPtr = getCoarserLevel();
    amrCoarserPtr->isThisFinestLev(!m_grids.size()>0);

#ifdef GRAVITY
    // in AMRC::timestep() we make sure that regrid() is called only
    // after postTimeStep() has been executed for all synchronized
    // levels. Thus this regrid() was called from the coarsest among
    // the synchronized levels, which, however, is not necessarily the
    // same as a_base_level ! In any case, this means that at this
    // point all regridding operations have been carried out.  Since
    // the grid hierarchy has changed, we now need to call gravity()
    // to compute the new multilevel gravitational potential. That
    // call must be made from the coarsest among the synchronized
    // levels: so first task is to climb down the ranks and find out
    // that level. But note that if m_level==1, the level we look for
    // can only be a_base_level=0.

    // we now we are in a_base_level+1;
    int crsestSynchLev = m_level-1;
    AMRLevelSelfGravity* amrSynchLevPtr = amrCoarserPtr;

    if (m_level>1)
    {
      // look for the right level
      const Real eps = 5.e-2;
      const Real smallT = eps * m_dt;
      AMRLevelSelfGravity* amrCoarserPtr = amrSynchLevPtr->getCoarserLevel();

      while (crsestSynchLev>0)
      {
        const Real dTimeLevs= amrSynchLevPtr->time()-amrCoarserPtr->time();
        if (abs(dTimeLevs) > smallT) break;

        // need to update allFinerPcls for all synchronized levels
        --crsestSynchLev;

        amrSynchLevPtr = amrCoarserPtr;
        amrCoarserPtr = amrSynchLevPtr->getCoarserLevel();
      }
    }

    // kind of need to call gravity() from the synch lev because of
    // reference to data of the base level there.
    const bool srceCorr= false;
    amrSynchLevPtr->gravity(crsestSynchLev,srceCorr);
#endif
  }

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::postRegrid done " << endl;
  }
}


// Initialize grids
void AMRLevelSelfGravity::initialGrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::initialGrid " << m_level << endl;
  }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  if (s_verbosity >= 4)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;

    pout() << "new grids: " << endl;
    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      pout() << constGrids[lit()] << endl;
    }
  }

  // Define old and new state data structures
  const IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  m_UOld.define(m_grids,m_numStates,ivGhost);
  resetToZero(m_UNew);
  resetToZero(m_UOld);

#ifdef GRAVITY
  m_phiNew.define(m_grids,1,m_numPhiGhost*IntVect::Unit);
  m_phiOld.define(m_grids,1,m_numPhiGhost*IntVect::Unit);
  resetToZero(m_phiNew);

  if (m_useDeltaPhiCorr)
  {
    m_deltaPhi.define(m_grids,1,ivGhost);
    resetToZero(m_deltaPhi);
  }

  m_forceNew.define(m_grids,SpaceDim,m_numForceGhost*IntVect::Unit);
  m_forceOld.define(m_grids,SpaceDim,m_numForceGhost*IntVect::Unit);

  const IntVect ivRhsGhost = m_numRhsGhost*IntVect::Unit;
  m_rhs.define(m_grids,1,ivRhsGhost);
  resetToZero(m_rhs);
#endif

  // Set up data structures
  levelSetup();
}

// Initialize data
void AMRLevelSelfGravity::initialData()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::initialData " << m_level << endl;
  }

  if (m_level>m_maxInitRefLevel) return;

  PhysIBC* physIBCPtr = m_gdnvPhysics->getPhysIBC();
  physIBCPtr->initialize(m_UNew);
}

// Things to do after initialization
void AMRLevelSelfGravity::postInitialize()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::postInitialize " << m_level << endl;
  }

  if (m_hasFiner)
  {
    AMRLevelSelfGravity* amrFinerPtr   = getFinerLevel();
    m_isThisFinestLev = (amrFinerPtr->m_grids.size()<=0);
  }

  // by the time the execution has reached this point, the finest level
  // should be the maxInitRefLevel, unless the requirements for refinement
  // were not met (or something else did not work)...
  if (!m_isThisFinestLev)
  {
    AMRLevelSelfGravity* amrFinerPtr = getFinerLevel();

    amrFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                 amrFinerPtr->m_UNew);

    // preserve smooth pressure distribution predating interpolation
    for (DataIterator di=m_grids.dataIterator(); di.ok(); ++di)
    {
      m_gdnvPhysics->setPressureToEntropy(m_UNew[di],m_grids.get(di));
    }
  }

#ifdef GRAVITY
  const bool srceCorr = false;
  if (m_level==0) gravity(m_level,srceCorr);
#endif
}


#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelSelfGravity::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::writeCheckpointHeader" << endl;
  }

  //
  a_handle.setGroup("/");

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

// Write checkpoint data for this level
void AMRLevelSelfGravity::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::writeCheckpointLevel" << endl;
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
  header.m_real["dx"]              = m_dx;
  header.m_real["dt"]              = m_dt;
  header.m_real["time"]            = m_time;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM(
         if (m_problem_domain.isPeriodic(0))
           header.m_int ["is_periodic_0"] = 1;
         else
           header.m_int ["is_periodic_0"] = 0; ,

         if (m_problem_domain.isPeriodic(1))
           header.m_int ["is_periodic_1"] = 1;
         else
           header.m_int ["is_periodic_1"] = 0; ,

         if (m_problem_domain.isPeriodic(2))
           header.m_int ["is_periodic_2"] = 1;
         else
           header.m_int ["is_periodic_2"] = 0; );

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  // write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_grids);
  write(a_handle,m_UNew,"data");
}

// Read checkpoint header
void AMRLevelSelfGravity::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::readCheckpointHeader" << endl;
  }

  //
  a_handle.setGroup("/");

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
    MayDay::Error("AMRLevelSelfGravity::readCheckpointHeader: checkpoint file does not have num_components");
  }
  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates)
  {
    MayDay::Error("AMRLevelSelfGravity::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
    {
      MayDay::Error("AMRLevelSelfGravity::readCheckpointHeader: checkpoint file does not have enough component names");
    }

    stateName = header.m_string[compStr];
    if (stateName != m_stateNames[comp])
    {
      MayDay::Error("AMRLevelSelfGravity::readCheckpointHeader: state_name in checkpoint does not match solver");
    }
  }
}

// Read checkpoint data for this level
void AMRLevelSelfGravity::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::readCheckpointLevel" << endl;
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
    MayDay::Error("AMRLevelSelfGravity::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelSelfGravity::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelSelfGravity::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelSelfGravity::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain time");
  }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
  {
    pout() << "read time = " << m_time << endl;
  }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelSelfGravity::readCheckpointLevel: file does not contain prob_domain");
  }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
           isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
         else
           isPeriodic[0] = false; ,

         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
           isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
         else
           isPeriodic[1] = false; ,

         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
           isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
         else
           isPeriodic[2] = false;);

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);

  if (gridStatus != 0)
  {
    MayDay::Error("AMRLevelSelfGravity::readCheckpointLevel: file does not contain a Vector<Box>");
  }

  // Create level domain
  m_grids = loadBalance(grids);

  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = m_grids;

  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
  {
    const Box& b = constGrids[lit()];
    m_level_grids.push_back(b);
  }

  if (s_verbosity >= 4)
  {
    pout() << "read level domain: " << endl;
    LayoutIterator lit = m_grids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = m_grids[lit()];
      pout() << lit().intCode() << ": " << b << endl;
    }
    pout() << endl;
  }

  // Reshape state with new grids
  m_UNew.define(m_grids,m_numStates);
  resetToZero(m_UNew);
  const int dataStatus = read<FArrayBox>(a_handle,
                                         m_UNew,
                                         "data",
                                         m_grids);

  if (dataStatus != 0)
  {
    MayDay::Error("AMRLevelSelfGravity::readCheckpointLevel: file does not contain state data");
  }
  // m_UOld.define(m_grids,m_numStates); // ,ivGhost);
  m_UOld.define(m_UNew);

  // Set up data structures
  levelSetup();

  // define finest level
  if (m_level>0 && m_hasCoarser)
  {
    AMRLevelSelfGravity* amrCrsePtr = getCoarserLevel();

    // the coarser level is not the finest level
    amrCrsePtr->isThisFinestLev(m_grids.size()<=0);
  }

#ifdef GRAVITY
  const IntVect ivPhiGhost = m_numPhiGhost*IntVect::Unit;
  m_phiNew.define(m_grids,1,ivPhiGhost);
  m_phiOld.define(m_grids,1,ivPhiGhost);
  resetToZero(m_phiNew);

  if (m_useDeltaPhiCorr)
  {
    m_deltaPhi.define(m_grids,1,ivPhiGhost);
    resetToZero(m_deltaPhi);
  }

  m_forceNew.define(m_grids,SpaceDim,m_numForceGhost*IntVect::Unit);
  m_forceOld.define(m_grids,SpaceDim,m_numForceGhost*IntVect::Unit);
  const IntVect ivRhsGhost = m_numRhsGhost*IntVect::Unit;
  m_rhs.define(m_grids,1,ivRhsGhost);
  resetToZero(m_rhs);

  // now I need to setup the gravitational field potential. I need to
  // call gravity() from the finest level to do so, but cannot rely on
  // "isThisFinestLev" till all levels have been read in. So I will
  // read this information directly from the header

  // Read the header for this level
  a_handle.setGroup("/");
  HDF5HeaderData tmpHeader;
  tmpHeader.readFromFile(a_handle);
  const int numLevels   = tmpHeader.m_int["num_levels"];
  const int finestLevel = numLevels-1;
  pout() << " N levs " << numLevels << endl;

  // not necessary but safer to call from the base lev
  if (m_level == finestLevel)
  {
    AMRLevelSelfGravity* baseLevPtr = this;
    while (baseLevPtr->getCoarserLevel() != NULL)
    {
      baseLevPtr = baseLevPtr->getCoarserLevel();
    }
    baseLevPtr->gravity(0,false);
  }

#endif
}

// Write plotfile header
void AMRLevelSelfGravity::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::writePlotHeader" << endl;
  }

  // reposition the hadle to root
  a_handle.setGroup("/");

  int numPlotNames = m_plotNames.size();

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = numPlotNames;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < numPlotNames; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_plotNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write plotfile data for this level
void AMRLevelSelfGravity::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::writePlotLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx;
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  const int numPlots = m_plotNames.size();
  LevelData<FArrayBox> UPlot(m_grids,numPlots,m_numGhost*IntVect::Unit);

#ifdef GRAVITY
  Interval phiInterval(numPlots-1,numPlots-1); // one before last
  m_phiNew.copyTo(m_phiNew.interval(),UPlot,phiInterval);
#endif

  //
  computePlotVars(UPlot,m_UNew);

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  // write(a_handle,m_UNew.boxLayout());
  // write(a_handle,m_UNew,"data");

  // Write the data for this level
  write(a_handle,UPlot.boxLayout());
  write(a_handle,UPlot,"data");
}
#endif

// Returns the dt computed earlier for this level
Real AMRLevelSelfGravity::computeDt()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::computeDt " << m_level << endl;
  }

  Real newDt;
  newDt = m_dtNew;

  return newDt;
}

// Compute dt using initial data
Real AMRLevelSelfGravity::computeInitialDt()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::computeInitialDt " << m_level << endl;
  }

#if 0 // this is the old way to do this
  Real newDTLocal = 1.e-2; // = 1.e-2;

  Real newDTHydr = m_initial_dt_multiplier * m_dx
                 / m_levelGodunov.getMaxWaveSpeed(m_UNew);
  if (s_verbosity >= 3) pout() << " newDTHydr = " << newDTHydr << endl;

  newDTLocal =  Min(newDTLocal, newDTHydr);
#endif

  LevelData<FArrayBox> source;
  // m_dt isn't used in the source term computation (although
  // it is required by the function), so it's safe to pass in the
  // undefined dt.
  setupSourceTerm(source, m_UNew, m_forceNew, m_time, m_dt);
  Real maxSpeed = 0;
  DataIterator dit = source.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = m_grids.get(dit());
      const FArrayBox& thisU = m_UNew[dit()];
      const FArrayBox& thisSrc = source[dit()];
      Real localMax = m_gdnvPhysics->getMaxWaveSpeedWithSource(thisU,
                                                               thisSrc,
                                                               m_dx,
                                                               gridBox);
      if (abs(localMax) > maxSpeed) maxSpeed = abs(localMax);
    } // end loop over grids

  // ensure that maxSpeed is not zero
  if (maxSpeed == 0.0)
  {
    //pout() << "Warning! maxSpeed is zero in AMRLevelSelfGravity::computeInitialDt!" << endl;
    maxSpeed = 1.0e-8;
  }
  Real newDtLocal = m_initial_dt_multiplier*m_dx/maxSpeed;

  Real newDT = globalMin(newDtLocal);

  if (s_verbosity >= 3)
  {
    pout() << " newDT = " << newDT << endl;
  }

  return newDT;
}


///
const LevelData<FArrayBox>& AMRLevelSelfGravity::getStateNew() const
{
  CH_assert(allDefined());
  return m_UNew;
}

///
const LevelData<FArrayBox>& AMRLevelSelfGravity::getStateOld() const
{
  CH_assert(allDefined());
  return m_UOld;
}

// set boolean stating whether or not this is the finest lev.
void AMRLevelSelfGravity::isThisFinestLev(const bool a_isThisFinestLev)
{
  m_isThisFinestLev = a_isThisFinestLev;
}

// return boolean stating whether or not this is the finest lev.
bool AMRLevelSelfGravity::isThisFinestLev() const
{
  return m_isThisFinestLev;
}

bool AMRLevelSelfGravity::allDefined() const
{
  return isDefined() && m_paramsDefined;
}

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelSelfGravity::loadBalance(const Vector<Box>& a_grids)
{
  CH_assert(allDefined());

  // Load balance and create boxlayout
  Vector<int> procMap;

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelSelfGravity::loadBalance: procesor map: " << endl;
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

// Setup menagerie of data structures
void AMRLevelSelfGravity::levelSetup()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << " AMRLevelSelfGravity::levelSetup " << m_level << endl;
  }

  AMRLevelSelfGravity* amrGodCoarserPtr = getCoarserLevel();
  AMRLevelSelfGravity* amrGodFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrGodCoarserPtr != NULL);
  m_hasFiner   = (amrGodFinerPtr   != NULL);

  // this will be finally set in postInitialize or postRegrid because
  // I can't really tell till all levels have been defined
  m_isThisFinestLev=true;

  if (m_hasCoarser)
  {
    int nRefCrse = m_coarser_level_ptr->refRatio();

    const DisjointBoxLayout& crseGrids = amrGodCoarserPtr->m_grids;

    m_coarseAverage.define(m_grids,
                           m_numStates,
                           nRefCrse);

    m_fineInterp.define(m_grids,
                        m_numStates,
                        nRefCrse,
                        m_problem_domain);

    // Maintain levelGodunov
    m_levelGodunov.define(m_grids,
                          crseGrids,
                          m_problem_domain,
                          nRefCrse,
                          m_dx,
                          m_gdnvPhysics,
                          m_normalPredOrder,
                          m_useFourthOrderSlopes,
                          m_usePrimLimiting,
                          m_useCharLimiting,
                          m_useFlattening,
                          m_useArtificialViscosity,
                          m_artificialViscosity,
                          m_hasCoarser,
                          m_hasFiner);

#ifdef GRAVITY
    const int numComp = 1;
    m_coarseAveragePhi.define(m_grids,
                              numComp,
                              nRefCrse);

    m_fineInterpPhi.define(m_grids,
                           numComp,
                           nRefCrse,
                           m_problem_domain);

    // Maintain QuadCFInterp
    m_quadCFInterp.define(m_grids, &crseGrids, m_dx,
                          nRefCrse, numComp, m_problem_domain);

    m_forcePatcher.define(m_grids, crseGrids, SpaceDim,
                          amrGodCoarserPtr->m_problem_domain,
                          amrGodCoarserPtr->m_ref_ratio,
                          m_numForceGhost);
#endif

    // This may look twisted but you have to do this this way because
    // the coarser levels get setup before the finer levels so, since
    // a flux register lives between this level and the next FINER
    // level, the finer level has to do the setup because it is the
    // only one with the information at the time of construction.

    // Maintain flux registers
    amrGodCoarserPtr->m_fluxRegister.define(m_grids,
                                            amrGodCoarserPtr->m_grids,
                                            m_problem_domain,
                                            amrGodCoarserPtr->m_ref_ratio,
                                            m_numStates);

    amrGodCoarserPtr->m_fluxRegister.setToZero();
  }
  else
  {
    m_levelGodunov.define(m_grids,
                          DisjointBoxLayout(),
                          m_problem_domain,
                          m_ref_ratio,
                          m_dx,
                          m_gdnvPhysics,
                          m_normalPredOrder,
                          m_useFourthOrderSlopes,
                          m_usePrimLimiting,
                          m_useCharLimiting,
                          m_useFlattening,
                          m_useArtificialViscosity,
                          m_artificialViscosity,
                          m_hasCoarser,
                          m_hasFiner);
  }
}

// Get the next coarser level
AMRLevelSelfGravity* AMRLevelSelfGravity::getCoarserLevel() const
{
  CH_assert(allDefined());

  AMRLevelSelfGravity* amrGodCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrGodCoarserPtr = dynamic_cast<AMRLevelSelfGravity*>(m_coarser_level_ptr);

    if (amrGodCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelSelfGravity::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrGodCoarserPtr;
}

// Get the next finer level
AMRLevelSelfGravity* AMRLevelSelfGravity::getFinerLevel() const
{
  CH_assert(allDefined());

  AMRLevelSelfGravity* amrGodFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrGodFinerPtr = dynamic_cast<AMRLevelSelfGravity*>(m_finer_level_ptr);

    if (amrGodFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelSelfGravity::getFinerLevel: dynamic cast failed");
    }
  }

  return amrGodFinerPtr;
}

// return pointer to
LevelData<FArrayBox>* AMRLevelSelfGravity::getPhi(const Real& a_time)
{
  if (s_verbosity >= 3)
    {
      pout() << " AMRLevelSelfGravity::getPhi " << m_level << " a_time= " << a_time << " m_time " << m_time << endl;
    }
  LevelData<FArrayBox>* phiPtr = NULL;

  const Real eps = 0.01*m_dt;
  if (Abs(a_time-m_time) <= eps) // case alpha=1; new synchronization point
    {
      phiPtr = &m_phiNew;
    }
  else if (Abs(a_time-(m_time-m_dtOld)) < eps) // case alpha=0; old synch point
    {
      phiPtr = &m_phiOld;
    }
  else
    {
      // define phiInt
      m_phiInt.define(m_grids,1,m_numPhiGhost*IntVect::Unit);
      //m_phiInt.define(m_grids,1);

      // need time interpolation
      Real alpha = (a_time-(m_time-m_dtOld))/m_dtOld;

      interpolateInTime(m_phiInt,m_phiOld,m_phiNew,
                        a_time,m_time,m_time-m_dtOld,m_dt,
                        Interval(0,0),Interval(0,0));

      // add deltaPhi
      if (m_useDeltaPhiCorr && !m_isThisFinestLev)
        {
          for (DataIterator di = m_phiInt.dataIterator(); di.ok(); ++di)
            {
              m_phiInt[di].plus(m_deltaPhi[di],alpha);
            }
        }

      phiPtr = &m_phiInt;
    }

  if (phiPtr == NULL)
    {
      MayDay::Error("AMRLevelSelfGravity::getPhi: something failed");
    }

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelSelfGravity:: getPhi done " << endl;
    }
  return phiPtr;
}


LevelData<FArrayBox>* AMRLevelSelfGravity::getPoissonRhs(const Real& a_time)
{
  if (s_verbosity >= 3)
  {
    pout() << " AMRLevelSelfGravity:: getPoissonRhs " << m_level << endl;
  }
  LevelData<FArrayBox>* poissonRhs = NULL;

  makePoissonRhs(m_rhs,m_UNew,a_time);

  poissonRhs = &m_rhs;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::getPoissonRhs done " << endl;
  }
  return poissonRhs;
}


//
void AMRLevelSelfGravity::gravity(const int a_baseLevel, const bool a_srceCorr)
{
  CH_assert(allDefined());

#ifdef GRAVITY
  if (s_verbosity >= 3)
    {
      pout() << " AMRLevelSelfGravity::gravity start: baseLevel = " << a_baseLevel << endl;
    }

  ellipticSolver(a_baseLevel,false);

  // need to do the two calls in this order, to ensure that boundary
  // conditions can be properly applied
  if (a_srceCorr)
    {
      // applies 2nd order corrections and compute the new force
      secondOrderCorrection();

      // same thing for finer levels
      AMRLevelSelfGravity* thisSelfGravityPtr = this;
      while (!thisSelfGravityPtr->isThisFinestLev())
        {
          thisSelfGravityPtr = thisSelfGravityPtr->getFinerLevel();

          thisSelfGravityPtr->secondOrderCorrection();
        }
    }
  else
    {
      // compute the new force
      computeForce(m_forceNew,m_phiNew,m_time);

      // same thing for finer levels
      AMRLevelSelfGravity* thisSelfGravityPtr = this;
      while (!thisSelfGravityPtr->isThisFinestLev())
        {
          thisSelfGravityPtr = thisSelfGravityPtr->getFinerLevel();

          thisSelfGravityPtr->computeForce(thisSelfGravityPtr->m_forceNew,
                                     thisSelfGravityPtr->m_phiNew,
                                     m_time);
        }
    }

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelSelfGravity::gravity done " << endl << endl;
    }
#endif   // gravity
}


//
void AMRLevelSelfGravity::ellipticSolver(const int  a_baseLevel,
                                         const bool a_isLevelSolve)
{
  CH_assert(allDefined());

#ifdef GRAVITY
  if (s_verbosity >= 3)
    {
      pout() << " AMRLevelSelfGravity::ellipticSolver start: baseLevel = " << a_baseLevel;
    }

  AMRLevelSelfGravity* thisSelfGravityPtr = this;

  // if this is a levelSolve then consider only one (this) level

  if (!a_isLevelSolve)
    {
      while (thisSelfGravityPtr->m_hasFiner && !thisSelfGravityPtr->isThisFinestLev())
        {
          thisSelfGravityPtr = thisSelfGravityPtr->getFinerLevel();
        }
    }

  const int finestLevel = thisSelfGravityPtr->m_level;

  if (s_verbosity >= 3)
    {
      pout() << "... and finestLevel = " << finestLevel << endl;
    }

  // up to finest level
  Vector<AMRLevelSelfGravity*>      amrLevel(finestLevel+1,NULL);
  Vector<LevelData<FArrayBox>*> amrPhi(finestLevel+1,NULL);

  const int bndryLevel = Max(0,a_baseLevel-1);

  // setup level pointers, up to finest level
  for (int lev = finestLevel; lev >= 0; --lev)
    {
      amrLevel[lev] = thisSelfGravityPtr;
      thisSelfGravityPtr = thisSelfGravityPtr->getCoarserLevel();
    }

  // set up potential, up to finest level
  for (int lev = bndryLevel; lev <= finestLevel; ++lev)
    {
      amrPhi[lev] = amrLevel[lev]->getPhi(m_time);
    }

  // up to max gravity level
  Vector<LevelData<FArrayBox>*> amrRhs(finestLevel+1,NULL);
  Vector<DisjointBoxLayout>   amrBoxes(finestLevel+1);
  Vector<ProblemDomain>     amrProbDom(finestLevel+1);
  Vector<Real>                   amrDx(finestLevel+1);
  Vector<int>                amrRefRat(finestLevel+1,1);

  // setup data structure pointers, up to max gravity level
  for (int lev = finestLevel; lev >= 0; --lev)
    {
      amrDx[lev]      = amrLevel[lev]->m_dx;
      amrBoxes[lev]   = amrLevel[lev]->m_grids;
      amrRefRat[lev]  = amrLevel[lev]->refRatio();
      amrProbDom[lev] = amrLevel[lev]->problemDomain();
    }


  // setup Poisson eq. right hand side, up to max gravity level
  for (int lev = a_baseLevel; lev <= finestLevel; ++lev)
    {
      amrRhs[lev]=amrLevel[lev]->getPoissonRhs(m_time);
    }

  bool isDomainCovered = (a_baseLevel==0);
  if (m_problem_domain.isPeriodic())
    {
      // should be: sum(rhs)=0; if not, then offset residual
      // and synchronize
      if (!isDomainCovered)
        {
          long numPtsDomain = amrProbDom[a_baseLevel].domainBox().numPts();

          // count number of cells on this level
          long numPtsLevel = 0;
          const DisjointBoxLayout& baseGrids = amrBoxes[a_baseLevel];
          for (LayoutIterator lit = baseGrids.layoutIterator(); lit.ok(); ++lit)
            {
              numPtsLevel += amrBoxes[a_baseLevel][lit].numPts();
            }

          isDomainCovered = (numPtsDomain==numPtsLevel);
        }

      // compute sum(rhs)
      if (!a_isLevelSolve)
        {
          if (isDomainCovered)
            {
              m_rhsOffset = computeSum(amrRhs, amrRefRat,
                                       amrDx[a_baseLevel],
                                       Interval(0,0),
                                       a_baseLevel);

              // divide offset by domain volume
              Real domainVol = pow(m_domainLength, SpaceDim);
              m_rhsOffset /= domainVol;

              if (s_verbosity >= 3)
                {
                  pout() << " gravity::rhs_resid: " << m_rhsOffset
                         << " level " << a_baseLevel << endl;
                }

              for (int lev = finestLevel; lev > a_baseLevel; --lev)
                {
                  amrLevel[lev]->m_rhsOffset = m_rhsOffset;
                }
            }
          else
            {
              // we needed this if the levels were created after
              // rhsOffset was computed
              for (int lev = finestLevel; lev > bndryLevel; --lev)
                {
                  amrLevel[lev]->m_rhsOffset = amrLevel[bndryLevel]->m_rhsOffset;
                }
            }

          // enforce sum(rhs)=0
          for (int lev = a_baseLevel; lev <= finestLevel; ++lev)
            {
              offset(*amrRhs[lev], m_rhsOffset);
            }

          if (s_verbosity >= 3)
            {
              if (isDomainCovered)
                {
                  // compute the new sum(rhs)
                  pout() << " gravity: rhs_resid: lev= : " << a_baseLevel
                         << " " << computeSum(amrRhs, amrRefRat,
                                              amrDx[a_baseLevel],
                                              Interval(0,0),
                                              a_baseLevel) << endl;
                }
            }
        }
      else if (isDomainCovered)
        {
          Real rhsOffset = computeSum(*amrRhs[a_baseLevel],
                                      NULL,1,
                                      amrDx[a_baseLevel]);
          // divide offset by domain volume
          Real domainVol = pow(m_domainLength, SpaceDim);
          rhsOffset /= domainVol;
          offset(*amrRhs[a_baseLevel], rhsOffset);
        }
    }


  bool reset = (a_baseLevel==0 ? true : SOLVER_RESET_PHI);

  AMRPoissonOpFactory opFactory;
  opFactory.define(amrProbDom[0],
                   amrBoxes,
                   amrRefRat,
                   amrDx[0],
                   &NoOpBc);

#ifdef USE_RELAX_SOLVER
  RelaxSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_imax = 4*SOLVER_NUM_SMOOTH;
  bottomSolver.m_eps = 1.0e-3;
  bottomSolver.m_verbosity = s_verbosity;
#else
  BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = s_verbosity;
  bottomSolver.m_normType  = SOLVER_NORM_TYPE;
#endif

  //
  AMRMultiGrid<LevelData<FArrayBox> > amrSolver;

  amrSolver.define(amrProbDom[0],opFactory,&bottomSolver,finestLevel+1);

  amrSolver.setSolverParameters(SOLVER_NUM_SMOOTH,SOLVER_NUM_SMOOTH,
                                SOLVER_NUM_SMOOTH,SOLVER_NUM_MG,
                                SOLVER_MAX_ITER,SOLVER_TOLERANCE,
                                SOLVER_HANG,SOLVER_NORM_THRES);

  amrSolver.solve(amrPhi,amrRhs,finestLevel,a_baseLevel,reset);

  if (!a_isLevelSolve)
    {
      // do postelliptic kind of operations
      for (int lev = finestLevel; lev > a_baseLevel; --lev)
        {
          // Average Phi from finer level data
          amrLevel[lev]->m_coarseAveragePhi.averageToCoarse(*amrPhi[lev-1],
                                                            *amrPhi[lev]);
        }

      if (m_useDeltaPhiCorr)
        {
          //
          for (int lev = finestLevel-1; lev >= a_baseLevel; --lev)
            {
              //
              LevelData<FArrayBox>& phi  = (*amrPhi[lev]);
              LevelData<FArrayBox>& dPhi = (amrLevel[lev]->m_deltaPhi);

              // save a copy of phi in dPhi and use phi with the
              // solver which needs ghost data
              for (DataIterator di=phi.dataIterator(); di.ok(); ++di)
                {
                  dPhi[di].copy(phi[di]);
                }

              // point amrPhi[lev] to deltaPhi; phi should still point
              // to m_phiNew of level lev
              //                  amrPhi[lev] = &(dPhi);

              // undo the previous offsetting, simlar to levelSolve case
              if (m_problem_domain.isPeriodic())
                {
                  // if (m_level==0)
                  if (lev==a_baseLevel && isDomainCovered)
                    {
                      offset((*amrRhs[lev]),
                             computeSum((*amrRhs[lev]),NULL,1,amrDx[lev]));
                    }
                  else
                    {
                      offset((*amrRhs[lev]),-m_rhsOffset);
                    }
                }

              // solve again, but on a single level
              amrSolver.solve(amrPhi,amrRhs,lev,lev,reset);

              //
              Real normDeltaPhi=zero;
              const DisjointBoxLayout& grids = amrBoxes[lev];
              for (DataIterator di = grids.dataIterator(); di.ok(); ++di)
                {
                  // temp to allow data swapping
                  FArrayBox temp(grids.get(di),1);

                  // pass copy of phi to temp
                  temp.copy(dPhi[di]);

                  // define dPhi = phi^comp-phi^level
                  dPhi[di] -= phi[di];

                  // restore phi
                  phi[di].copy(temp);

                  normDeltaPhi += dPhi[di].norm(grids.get(di));
                }

              //
              if (s_verbosity >= 3)
                {
                  pout() << " normDeltaPhi = " << normDeltaPhi << endl;
                }
            }
        }
    } // if !levelSolve

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelSelfGravity::ellipticSolver done " << endl << endl;
    }
#endif  // elliptic solver
}

// Allow for second order source term corrections
void AMRLevelSelfGravity::secondOrderCorrection()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::secondOrderCorrection: " << m_level << endl;
  }

  IntVect ivGhost = m_numForceGhost*IntVect::Unit;
  LevelData<FArrayBox> dForce(m_grids,SpaceDim,ivGhost);

  // first set to old force
  DataIterator di = m_grids.dataIterator();
  for (di.begin(); di.ok(); ++di)
  {
    dForce[di].copy(m_forceOld[di]);
  }

  // then compute new force
  computeForce(m_forceNew,m_phiNew,m_time);

  for (di.begin(); di.ok(); ++di)
  {
    dForce[di] -= m_forceNew[di];
    dForce[di] *= (-half);
  }

  for (di.begin(); di.ok(); ++di)
  {
    m_gdnvPhysics->forceCorrection(m_UNew[di],dForce[di],m_time,m_dt,
                                   m_grids.get(di));
  }

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::secondOrderCorrection: done " << endl;
  }
}

//
void AMRLevelSelfGravity::makePoissonRhs(LevelData<FArrayBox>&       a_rhs,
                                         const LevelData<FArrayBox>& a_U,
                                         const Real&                 a_time)
{
#ifdef GRAVITY
  // I am assuming they are all defined on the same grids; assert() this along the way

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelSelfGravity::makePoissonRhs: " << m_level << endl;
  }

  const DisjointBoxLayout& grids = a_rhs.getBoxes();
  DataIterator di = grids.dataIterator();

  // 0. set rhs = 0
  resetToZero(a_rhs);

  // 1. add gas
  if (a_U.isDefined())
  {
    CH_assert(a_rhs.getBoxes()==a_U.getBoxes());
    for (di.begin(); di.ok(); ++di)
    {
      a_rhs[di].plus(a_U[di],grids.get(di),0,0);
    }
  }

  int numCells = 0;
  for (di.begin(); di.ok(); ++di)
  {
    const Box& curBox = grids.get(di);
    numCells += curBox.numPts();

    if (m_problem_domain.isPeriodic()) a_rhs[di].plus(-one,curBox);
    a_rhs[di].mult(threehalf,curBox);
  }
  if (s_verbosity > 1)
  {
    const bool notPerVol = false;
    const Real totalRhs = globalAverage(a_rhs,0,notPerVol);

    pout() << "  " << endl;
    pout() << " level  " << m_level  << endl;
    pout() << " Total   RHS   = " << totalRhs << endl;
    pout() << " # cells       = " << numCells << endl;
    pout() << "  " << endl;
  }

  if (s_verbosity >= 3)
  {
    pout() << " AMRLevelSelfGravity::makePoissonRhs: done  " << endl;
  }
#endif  // gravity
}


// add gravitational source term correction.
void AMRLevelSelfGravity::setupSourceTerm(LevelData<FArrayBox>&       a_source,
                                          const LevelData<FArrayBox>& a_U,
                                          const LevelData<FArrayBox>&  a_force,
                                          const Real&                 a_time,
                                          const Real&                 a_dt)
{
#ifdef GRAVITY
  if (s_verbosity >= 4)
  {
    pout() << " AMRLevelSelfGravity::computeSourceTerm " << m_level << endl;
  }

  a_source.define(m_grids,m_gdnvPhysics->numPrimitives(),2*IntVect::Unit);

  for (DataIterator di = m_grids.dataIterator(); di.ok(); ++di)
  {
    // reset to zero
    a_source[di].setVal(zero);

    // box where you need the source term defined
    Box sbox = grow(m_grids.get(di),2);

    CH_assert(a_force[di].box().contains(sbox));
    m_gdnvPhysics->setupSourceTerm(a_source[di],a_force[di],a_U[di],
                                   a_time,a_dt,sbox);
  }

  if (s_verbosity >= 4)
  {
    pout() << " AMRLevelSelfGravity::setupSourceTerm done " << endl;
  }
#endif
}

// overloaded version that can be called from gravity()
void AMRLevelSelfGravity::computeForce()
{
  computeForce(m_forceNew,m_phiNew,m_time);
}

//FM 6.7.05: a_phi=>phi: temporary hack around solver misbehavior when phi ghosts>1
// function that computes F=-grad(phi)
void AMRLevelSelfGravity::computeForce(LevelData<FArrayBox>&   a_force,
                                       LevelData<FArrayBox>&  a_phi,
                                       const Real&            a_time)
{
  CH_assert(allDefined());

#ifdef GRAVITY
  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelSelfGravity::computeForce " << m_level << endl;
  }

  //FM 6.7.05: temporary hack around solver misbehavior when phi ghosts>1
  int numPhiGhost = 2;
  //LevelData<FArrayBox> a_phi(m_grids,1,numPhiGhost*IntVect::Unit);
  //phi.copyTo(a_phi.interval(),a_phi,a_phi.interval());

  a_phi.exchange(a_phi.interval());

  Real alpha;
  AMRLevelSelfGravity* amrCoarserPtr = NULL;
  if (m_hasCoarser)
  {
    amrCoarserPtr = getCoarserLevel();
    const Real tNew = amrCoarserPtr->m_time;
    const Real tOld = tNew - amrCoarserPtr->m_dt;

    const Real eps = 0.01*m_dt;
    if (Abs(a_time-tNew) <= eps) // case alpha=1; synchronization point;
    {
      alpha = one;
      m_quadCFInterp.coarseFineInterp(a_phi,amrCoarserPtr->m_phiNew);
      // ,numPhiGhost);
    }
    else if (Abs(a_time-tOld) < eps)  // case alpha=0
    {
      alpha = zero;
      pout() << " Warning:: alpha=0 should never happen ";
      pout() << " time " << a_time << " tOld " << tOld << " tNew " << tNew;
      pout() << " dt " << m_dt << endl;
      m_quadCFInterp.coarseFineInterp(a_phi,amrCoarserPtr->m_phiOld);
      // ,numPhiGhost);
    }
    else
    {
      CH_assert((tNew-tOld)>eps);
      alpha= (a_time-tOld)/(tNew-tOld);
      LevelData<FArrayBox> phiCrse(amrCoarserPtr->m_grids,1,a_phi.ghostVect());
      interpolateInTime(phiCrse,amrCoarserPtr->m_phiOld,
                        amrCoarserPtr->m_phiNew, alpha,one,zero,m_dt);
      m_quadCFInterp.coarseFineInterp(a_phi,phiCrse); //,numPhiGhost);
    }
  }

  const DisjointBoxLayout& grids = a_phi.getBoxes();
  for (DataIterator di = grids.dataIterator(); di.ok(); ++di)
  {
    // The current box
    const Box& box = grids.get(di);

    // f = -grad(phi), hence pass -dx
    m_gradient->gradient(a_force[di],a_phi[di],m_problem_domain,(-m_dx),box);
  }

  if (m_hasCoarser)
  {
    m_forcePatcher.fillInterp(a_force,
                   amrCoarserPtr->m_forceOld,
                   amrCoarserPtr->m_forceNew,
                   alpha,0,0,SpaceDim);

    if (s_verbosity >= 4)
    {
      pout() << "AMRLevelSelfGravity:: computeForce CFInterp OK!" << endl;
    }
  }

  // exchange here to get force on one zone outside the valid region:
  // needed for the predictor step and the source term
  a_force.exchange(a_force.interval());

  /*
  if (s_verbosity >= 3 && s_verbosity < 4)
  {
    for (DataIterator di = grids.dataIterator(); di.ok(); ++di)
    {
      const FArrayBox& fx = a_force[di].getFlux(0);
      const Box& box = fx.box();

      for (BoxIterator bi(box); bi.ok(); ++bi)
      {
        const Real& f = fx(bi(),0);
        pout() << " IV : ";
        for (int i=0;i<SpaceDim;i++) pout() << bi()[i] << " ";
        pout() << "    fx : " << f << endl;
      }
    }
  }
  */

  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelSelfGravity::computeForce done "  << endl;
  }
#endif
}
