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

#include <string>
#include "parstream.H"

#include "AMRLevelMappedCons.H"

#include "FArrayBox.H"
#include "NodeFArrayBox.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "computeSum.H"
#include "FourthOrderFillPatch.H"
#include "computeNorm.H"
#include "CellToEdge.H"
#include "AMRIO.H"
#include "NodeAMRIO.H"
#include "TimeInterpolatorRK4.H"
#include "PatchGodunovF_F.H"
#include "WaveVelocityF_F.H"
// #include "Divergence.H"
#include "LevelRK4.H"
#include "AdvectOpF_F.H"
#include "PolytropicPhysicsF_F.H"
#include "AMRLevel.H"
#include "FourthOrderUtil.H"
#include "computeMappedDt.H"
// #include "DebugInclude.H"
// for kludge
#include "PhysMappedIBC.H"

//////////////////////////////////////////////////////////////////////////////

// Constructor
AMRLevelMappedCons::AMRLevelMappedCons()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons default constructor" << endl;
  }
  m_levelStep = 0;
  m_cfl = 0.8;
  m_spaceOrder = 4;
  m_limitFaceValues = true;
  m_initialAverage = false;
  m_useFlattening = false;
  m_useArtVisc = false;
  m_ratioArtVisc = 0.;
  m_forwardEuler = false;
  m_enforceMinVal = true;
  m_noPPM = false;
  m_doDeconvolution = true;
  m_doFaceDeconvolution = true;
  m_useArtificialViscosity = false;
  m_minVal = -100000.0;
  m_domainLength = 1.0;
  m_refineThresh = 0.2;
  m_initial_dt_multiplier = 0.1;
  m_gdnvPhysics = NULL;
  m_coordSysFactPtr = NULL;
  m_dtFromCells = false;
}

//////////////////////////////////////////////////////////////////////////////

// Destructor
AMRLevelMappedCons::~AMRLevelMappedCons()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons destructor" << endl;
    }

  if (m_gdnvPhysics != NULL)
    {
      delete m_gdnvPhysics;
      m_gdnvPhysics = NULL;
    }
  if (m_coordSysPtr != NULL)
    {
      delete m_coordSysPtr;
      m_coordSysPtr = NULL;
    }
}

//////////////////////////////////////////////////////////////////////////////

// Define new AMR level
void AMRLevelMappedCons::define(AMRLevel*            a_coarserLevelPtr,
                          const ProblemDomain& a_problemDomain,
                          int                  a_level,
                          int                  a_refRatio)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::define " << a_level << endl;
    }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().longside();

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelMappedCons* amrConsPtr = dynamic_cast<AMRLevelMappedCons*>(a_coarserLevelPtr);

    if (amrConsPtr != NULL)
    {
      m_cfl = amrConsPtr->m_cfl;
      m_spaceOrder = amrConsPtr->m_spaceOrder;
      m_limitFaceValues = amrConsPtr->m_limitFaceValues;
      m_useFlattening = amrConsPtr->m_useFlattening;
      m_useArtVisc = amrConsPtr->m_useArtVisc;
      m_ratioArtVisc = amrConsPtr->m_ratioArtVisc;
      m_enforceMinVal = amrConsPtr->m_enforceMinVal;
      m_noPPM = amrConsPtr->m_noPPM;
      m_doDeconvolution = amrConsPtr->m_doDeconvolution;
      m_doFaceDeconvolution = amrConsPtr->m_doFaceDeconvolution;
      m_useArtificialViscosity = amrConsPtr->m_useArtificialViscosity;
      m_artificialViscosity = amrConsPtr->m_artificialViscosity;
      m_minVal = amrConsPtr->m_minVal;
      m_domainLength = amrConsPtr->m_domainLength;
      m_refineThresh = amrConsPtr->m_refineThresh;
      m_tagBufferSize = amrConsPtr->m_tagBufferSize;
    }
    else
    {
      MayDay::Error("AMRLevelMappedCons::define: a_coarserLevelPtr is not castable to AMRLevelMappedCons*");
    }
  }

  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more

  // Begin application-dependent code - dfm

  // if we're limiting, need 3 ghost cells in U. otherwise, 2 will do
  // NEW, petermc, 6 July 2008:  from 3 to 4
  // changed from 4 to 5 on 20 Aug 2009
  m_numGhost = 5;
  m_ghostVect = m_numGhost * IntVect::Unit;
  // NEW, petermc, 25 June 2008
  //  if (m_useArtVisc) m_numGhost += 2;

  //   m_numStates = 1;
  //   m_stateNames.push_back("U_0");
  //m_stateNames.push_back("U_1");

  // NEW, petermc, 11 June 2008:  what is m_gdnvPhysics ?
  // It is GodunovPhysics, and defined in the executable where
  // it is set to PolytropicPhysics.
  m_gdnvPhysics->define(m_problem_domain, m_dx);
  // Number and names of conserved states
  m_numStates  = m_gdnvPhysics->numConserved();
  m_stateNames = m_gdnvPhysics->stateNames();
  // For writing checkpoint files, include state names with OLD in front.
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      char stateNameChars[60];
      sprintf(stateNameChars, "OLD%s", m_stateNames[comp].c_str());
      m_stateNames.push_back(stateNameChars);
    }

  // End application-dependent code - dfm
}

//////////////////////////////////////////////////////////////////////////////

// Advance by one timestep
Real AMRLevelMappedCons::advance()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::advance level " << m_level
           << " to time " << m_time << endl;
  }

  // Begin application-dependent code - dfm

  // Copy the new to the old
  //  m_Unew.copyTo(m_Unew.interval(),
  //                m_Uold,
  //                m_Uold.interval());
  // petermc, 13 Jan 2009:  Copy ghost cells too.
  m_levelConsOperator.copySolnData(m_Uold, m_Unew);

  // End application-dependent code - dfm

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  LevelData<FArrayBox>* coarserUold = &dummyData;
  LevelData<FArrayBox>* coarserUnew = &dummyData;
  // A coarser level exists
  if (m_hasCoarser)
    {
      AMRLevelMappedCons* coarserPtr = getCoarserLevel();

      // Recall that my flux register goes between my level and the next
      // finer level
      coarserFR = &coarserPtr->m_fluxRegister;

      coarserUold = &coarserPtr->m_Uold;
      coarserUnew = &coarserPtr->m_Unew;

      tCoarserNew = coarserPtr->m_time;
      tCoarserOld = tCoarserNew - coarserPtr->m_dt;
    }

  // A finer level exists
  if (m_hasFiner)
    {
      // Recall that my flux register goes between my level and the next
      // finer level
      finerFR = &m_fluxRegister;
      finerFR->setToZero();
    }

  // Advance conservation law by one time step using 4th-order Runge-Kutta.

  // Set m_levelConsOperator.m_fluxes = 0.
  m_levelConsOperator.setCoordSys(m_coordSysPtr);
  m_levelConsOperator.resetFluxes();

  // We need 1 layer of ghost cells of UoldCopy.
  LevelData<FArrayBox> UoldCopy(m_grids, m_numStates, m_ghostVect);
  m_levelConsOperator.copySolnData(UoldCopy, m_Uold);
  Real oldTime = m_time;
  if (m_forwardEuler)
    {
      LevelData<FArrayBox> RHSTmp;
      m_levelConsOperator.defineRHSData(RHSTmp, m_Uold);
      // start by copying old solution into new solution
      m_levelConsOperator.copySolnData(m_Unew, m_Uold);
      m_levelConsOperator.evalCountMax(1); // to make 1 call to evalRHS
      m_levelConsOperator.resetEvalCount();
      m_levelConsOperator.evalRHS(RHSTmp, m_Uold,
                                  *finerFR, *coarserFR,
                                  *coarserUold, tCoarserOld,
                                  *coarserUnew, tCoarserNew,
                                  m_time, m_dt);
      m_levelConsOperator.updateODE(m_Unew, RHSTmp, m_dt);
    }
  else
    {
      m_levelConsOperator.evalCountMax(4); // to make 4 calls to evalRHS
      m_levelConsOperator.resetEvalCount();
      if (m_hasFiner)
        {
          // Fill in the interpolator from this level to next finer level:
          // we'll use it later when we calculate on next finer level.
          TimeInterpolatorRK4& timeInterpolator =
            getFinerLevel()->m_levelConsOperator.getTimeInterpolator();

          RK4LevelAdvance<LevelData<FArrayBox>,
            TimeInterpolatorRK4,
            LevelFluxRegister,
            LevelMappedConsOperator>(m_Unew, m_Uold, timeInterpolator,
                                     *coarserUold, tCoarserOld,
                                     *coarserUnew, tCoarserNew,
                                     *coarserFR, *finerFR,
                                     m_time, m_dt, m_levelConsOperator);
        }
      else
        {
          RK4LevelAdvance<LevelData<FArrayBox>,
            LevelFluxRegister,
            LevelMappedConsOperator>(m_Unew, m_Uold,
                                     *coarserUold, tCoarserOld,
                                     *coarserUnew, tCoarserNew,
                                     *coarserFR, *finerFR,
                                     m_time, m_dt, m_levelConsOperator);
        }
    }

  if (m_useArtificialViscosity)
    {
      // We need 1 layer of ghost cells of UoldCopy.
      m_levelConsOperator.fillGhosts(UoldCopy, oldTime,
                                     tCoarserOld, tCoarserNew);

      // find flux, fluxArtificialViscosity, from artificial viscosity;
      // then m_fluxes += fluxArtificialViscosity;
      // and m_Unew -= div(fluxArtificialViscosity) * m_dt.
      m_levelConsOperator.addArtificialViscosity(m_Unew, UoldCopy,
                                                 *finerFR, *coarserFR,
                                                 oldTime, m_dt);
    }

  if (m_useArtVisc)
    {
      addDissipation();
    }

  // Update the time and store the new timestep
  ++m_levelStep;
  m_time += m_dt;
  // 27 Feb 2009:  Update timestep.
  // m_dtNew = m_dt;
  //  Real newDt = m_dx / getMaxWaveSpeed(m_Unew);
  //  m_dtNew = m_cfl * newDt;
  Real returnDt = computeNewDt();
  if (m_dtFromCells)
    { // the old method
      returnDt *= m_cfl;
    }
  m_dtNew = returnDt;
  return returnDt;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::addDissipation()
{
  // THIS NEEDS TO BE REWRITTEN FOR MAPPED COORDINATES

  // m_Unew += (m_ratioArtVisc / SpaceDim^3 * h^6) *
  //           laplacian ( laplacian (m_Uold) )
  Real mu = m_ratioArtVisc;
  for (int ilap = 0; ilap < 3; ilap++)
    { // apply Laplacian 3 times
      mu *= m_dx*m_dx / (SpaceDim * 1.);
    }

  // const Interval& intvlU = m_Uold.interval();
  const IntVect& ghostU = m_Uold.ghostVect();
  // if ghostU >= 2*IntVect::Unit then this is IntVect::Unit
  IntVect ghostLap = min(ghostU, IntVect::Unit);
  // if ghostU >= 2*IntVect::Unit then this is IntVect::Zero
  IntVect ghostLapLap = ghostLap - IntVect::Unit;
  int numU = m_Uold.nComp();

  const Box& domainBox = m_problem_domain.domainBox();
  const IntVect& domainLo = domainBox.smallEnd();
  const IntVect& domainHi = domainBox.bigEnd();
  DataIterator dit = m_Unew.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& boxBase = m_grids[dit];
      // petermc, 20 Jan 2009:
      // These box settings assume that UoldFab is defined everywhere.
      // But it isn't necessarily.
      //      Box boxLapU = grow(boxBase, ghostU - IntVect::Unit);
      //      Box boxLapLapU = grow(boxBase, ghostU - 2*IntVect::Unit);
      // Instead, go the reverse way:  set laplacian on as much of
      // baseBox as possible.
      Box boxLapU = grow(boxBase, ghostLap);
      Box boxLapLapU = grow(boxBase, ghostLapLap);
      for (int dir=0; dir<SpaceDim; dir++)
        if (!m_problem_domain.isPeriodic(dir))
          {
            if (boxBase.smallEnd(dir) == domainLo[dir])
              { // shrink boxLapU and boxLapLapU on the low end
                boxLapU.growLo(dir, -1);
                boxLapLapU.growLo(dir, -2);
              }
            if (boxBase.bigEnd(dir) == domainHi[dir])
              { // shrink boxLapU and boxLapLapU on the high end
                boxLapU.growHi(dir, -1);
                boxLapLapU.growHi(dir, -2);
              }
          }
      FArrayBox lapUfab(boxLapU, numU);
      lapUfab.setVal(0.);
      const FArrayBox& UoldFab = m_Uold[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_CCLAPLACIAN(CHF_FRA(lapUfab),
                           CHF_CONST_FRA(UoldFab),
                           CHF_BOX(boxLapU),
                           CHF_CONST_INT(dir),
                           CHF_CONST_REAL(m_dx));
        }
      FArrayBox lapLapUfab(boxLapLapU, numU);
      lapLapUfab.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_CCLAPLACIAN(CHF_FRA(lapLapUfab),
                           CHF_CONST_FRA(lapUfab),
                           CHF_BOX(boxLapLapU),
                           CHF_CONST_INT(dir),
                           CHF_CONST_REAL(m_dx));
        }
      FArrayBox& UnewFab = m_Unew[dit];
      UnewFab.plus(lapLapUfab,
                   boxLapLapU, boxLapLapU, // source and dest boxes
                   mu, // scaling
                   0, 0, // source and dest components
                   numU); // number of components
    }
}

//////////////////////////////////////////////////////////////////////////////

// Things to do after a timestep
void AMRLevelMappedCons::postTimeStep()
{
  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::postTimeStep " << m_level << endl;
    }

  // Begin application-dependent code - PC.

  if (m_hasFiner)
    {
      // Reflux
      Real scale = -1.0/m_dx;
      m_fluxRegister.reflux(m_Unew, scale);

      // Average from finer level data
      AMRLevelMappedCons* amrConsFinerPtr = getFinerLevel();

      // Set m_Unew at cell c to the mean of amrConsFinerPtr->m_Unew
      // over all finer cells within c.
      amrConsFinerPtr->
        m_coarseAverage.averageToCoarse(m_Unew, amrConsFinerPtr->m_Unew);
    }
  // End application-dependent code - PC.

  if (s_verbosity >= 2 && m_level == 0)
    {
      int nRefFine = 1;

      pout() << "AMRLevelMappedCons::postTimeStep:" << endl;
      pout() << "  Sums:" << endl;
      for (int comp = 0; comp < m_numStates; comp++)
        {
          Interval curComp(comp,comp);
          // Begin application-dependent code - PC.
          Real integral = computeSum(m_Unew, NULL, nRefFine, m_dx, curComp);
          // End application-dependent code - PC.

          pout() << "t = " << m_time
                 << "    " << setw(23)
                 << setprecision(16)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << integral
                 << " --- " << m_stateNames[comp];

          if (comp == 0 )
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
              if (first)
                {
                  orig_integral = integral;
                  first = false;
                }
              last_integral = integral;
            }
          // also restore formatting
          pout() << setw(13) << setprecision(6) << endl;
        }
    }
  if (s_verbosity >= 1 && m_level == 0)
    {
      // petermc, 25 June 2008:  we do not have exact solution to compare
      //      writeErrorNorms();
      //Real energy = computeEnergy();
      //pout() << " integrated energy = " << energy << endl ;
    }

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::postTimeStep " << m_level
             << " finished" << endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

// Create tags for regridding
void AMRLevelMappedCons::tagCells(IntVectSet& a_tags)
{
  // CH_assert(allDefined());
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::tagCells " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::tagCellsInit " << m_level << endl;
    }

  // Begin application-dependent code - PC.

  // Create tags based on undivided gradient of density
  const DisjointBoxLayout& levelDomain = m_Unew.disjointBoxLayout();
  IntVectSet localTags;

  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
  {
    const AMRLevelMappedCons* amrGodCoarserPtr = getCoarserLevel();
    const LevelData<FArrayBox>& UnewCoarser = amrGodCoarserPtr->m_Unew;
    FourthOrderFillPatch filler(levelDomain,
                                UnewCoarser.disjointBoxLayout(),
                                m_numStates,
                                amrGodCoarserPtr->m_problem_domain,
                                m_ref_ratio,
                                1, // to fill in 1 layer of ghost cells at this level, meaning that you need 1 layer of COARSE ghost cells
                                true); // fixed time
    // No interpolation in time.
    filler.fillInterp(m_Unew, UnewCoarser, 0, 0, m_numStates);
  }
  m_Unew.exchange(Interval(0,m_numStates-1));

  // Compute relative gradient
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b = levelDomain[dit()];
    FArrayBox gradFab(b,SpaceDim);
    const FArrayBox& UFab = m_Unew[dit()];

    for (int dir = 0; dir < SpaceDim; ++dir)
    {
      const Box bCenter = b & grow(m_problem_domain,-BASISV(dir));

      const Box bLo     = b & adjCellLo(bCenter,dir);
      const int hasLo = ! bLo.isEmpty();

      const Box bHi     = b & adjCellHi(bCenter,dir);
      const int hasHi = ! bHi.isEmpty();

      FORT_GETRELGRADF(CHF_FRA1(gradFab,dir),
                       CHF_CONST_FRA1(UFab,0),
                       CHF_CONST_INT(dir),
                       CHF_BOX(bLo),
                       CHF_CONST_INT(hasLo),
                       CHF_BOX(bHi),
                       CHF_CONST_INT(hasHi),
                       CHF_BOX(bCenter));
    }

    FArrayBox gradMagFab(b,1);
    FORT_MAGNITUDEF(CHF_FRA1(gradMagFab,0),
                    CHF_CONST_FRA(gradFab),
                    CHF_BOX(b));

    // Tag where gradient exceeds threshold
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      if (gradMagFab(iv) >= m_refineThresh)
      {
        localTags |= iv;
      }
    }
  }

  // End application-dependent code - PC.

  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;
  a_tags = localTags;
  /*
  Vector<IntVectSet> allTags;

  const int destProc = uniqueProc(SerialTask::compute);

  gather(allTags,localTags,destProc);

  if (procID() == uniqueProc(SerialTask::compute))
  {
    for (int i = 0; i < allTags.size(); ++i)
    {
      a_tags |= allTags[i];
    }
  }
  */
}

//////////////////////////////////////////////////////////////////////////////

// Create tags at initialization
void AMRLevelMappedCons::tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

//////////////////////////////////////////////////////////////////////////////

// Set up data on this level after regridding
void AMRLevelMappedCons::regrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::regrid " << m_level << endl;
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
          pout() << constGrids[lit] << endl;
        }
    }

  // Save data for later
  // Begin application-dependent code - PC.

  LevelData<FArrayBox> Uold;
  Uold.define(m_Unew);
  m_Unew.copyTo(m_Unew.interval(),
                Uold,
                Uold.interval());

  // Reshape state with new grids
  m_Unew.define(m_grids, m_numStates, m_ghostVect);
  m_Uold.define(m_grids, m_numStates, m_ghostVect);

  // Set up data structures
  levelSetup();

  // Re-create the coordinate system
  m_coordSysPtr->regrid(m_grids);

  // Interpolate from coarser level
  if (m_hasCoarser)
    {
      AMRLevelMappedCons* amrConsCoarserPtr = getCoarserLevel();
      m_fineInterp.interpToFine(m_Unew,
                                amrConsCoarserPtr->m_Unew);

      // Begin application-dependent code - PC.
    }

  // Copy from old state
  Uold.copyTo(Uold.interval(),
              m_Unew,
              m_Unew.interval());
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////

// Initialize grids
void AMRLevelMappedCons::initialGrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::initialGrid " << m_level << endl;
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
          pout() << constGrids[lit] << endl;
        }
    }

  // Define old and new state data structures
  // Begin application-dependent code - PC.

  m_Unew.define(m_grids, m_numStates, m_ghostVect);
  m_Uold.define(m_grids, m_numStates, m_ghostVect);
  // End application-dependent code - PC.

  // Set up data structures
  levelSetup();
}

//////////////////////////////////////////////////////////////////////////////

// Initialize data
void AMRLevelMappedCons::initialData()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::initialData " << m_level << endl;
    }
  // Begin application-dependent code - PC.

  PhysMappedIBC* physIBCPtr = (PhysMappedIBC*) m_gdnvPhysics->getPhysIBC();
  // Kludge!  Need this for mapped IBC.
  // ((GaussianMappedIBC*) physIBCPtr)->setCoordSys(m_coordSysPtr);
  physIBCPtr->setCoordSys(m_coordSysPtr);
  physIBCPtr->setTime(m_time);
  physIBCPtr->initialize(m_Unew);

  m_Unew.exchange();
  if (m_initialAverage)
    { // call to new function, petermc, 19 Dec 2008
      // FOR MAPPED, THIS MAY BE A NEW FUNCTION.
      fourthOrderAverage(m_Unew, m_problem_domain);
    }

  // End application-dependent code - PC.
}

//////////////////////////////////////////////////////////////////////////////

// Things to do after initialization
void AMRLevelMappedCons::postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::postInitialize " << m_level << endl;
    }

  if (s_verbosity >= 1 && m_level == 0)
    {
      // petermc, 25 June 2008:  we do not have exact solution to compare
      //      writeErrorNorms();
    }
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelMappedCons::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::writeCheckpointHeader" << endl;
    }

  // We write out all components of m_Unew and all components of m_Uold.

  // Set up the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates*2;

  // Set up the component names:  These already include the OLD ones.
  char compStr[30];
  for (int comp = 0; comp < m_numStates*2; ++comp)
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

//////////////////////////////////////////////////////////////////////////////

// Write checkpoint data for this level
void AMRLevelMappedCons::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::writeCheckpointLevel" << endl;
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

  // Write the data for this level:  m_Unew and then m_Uold.
  // What about ghosts?
  LevelData<FArrayBox> outData(m_Unew.getBoxes(),2*m_numStates);
  Interval interval0(0,m_numStates-1);
  Interval interval1(m_numStates,2*m_numStates-1);
  m_Unew.copyTo(interval0, outData, interval0);
  m_Uold.copyTo(interval0, outData, interval1);
  write(a_handle,outData.boxLayout());
  write(a_handle,outData,"data");
}

//////////////////////////////////////////////////////////////////////////////

// Read checkpoint header
void AMRLevelMappedCons::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::readCheckpointHeader" << endl;
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
      MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: checkpoint file does not have num_components");
    }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates*2)
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: num_components in checkpoint file does not match solver");
    }

  if (header.m_int.find("iteration") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain iteration");
    }
  m_levelStep = header.m_int ["iteration"];

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates*2; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      if (header.m_string.find(compStr) == header.m_string.end())
        {
          MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: checkpoint file does not have enough component names");
        }

      stateName = header.m_string[compStr];
      if (stateName != m_stateNames[comp])
        {
          MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: state_name in checkpoint does not match solver");
        }
    }
}


//////////////////////////////////////////////////////////////////////////////

// Read checkpoint data for this level
void AMRLevelMappedCons::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::readCheckpointLevel" << endl;
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
    MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
    {
      pout() << "read dt = " << m_dt << endl;
    }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not contain time");
    }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
    {
      pout() << "read time = " << m_time << endl;
    }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not contain prob_domain");
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
  const int gridStatus = read(a_handle, grids);

  if (gridStatus != 0)
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not contain a Vector<Box>");
    }

  // Create level domain
  m_grids = loadBalance(grids);

  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = m_grids;

  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = constGrids[lit];
      m_level_grids.push_back(b);
    }

  if (s_verbosity >= 4)
    {
      pout() << "read level domain: " << endl;
      LayoutIterator lit = constGrids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& b = constGrids[lit];
          pout() << lit().intCode() << ": " << b << endl;
        }
      pout() << endl;
    }

  // Reshape state with new grids
  LevelData<FArrayBox> inData ;
  inData.define(m_grids, m_numStates*2, m_ghostVect);
  const int dataStatus = read<FArrayBox>(a_handle, inData, "data", m_grids);
  if (dataStatus != 0)
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not contain all state data");
    }
  m_Unew.define(m_grids, m_numStates, m_ghostVect);
  m_Uold.define(m_grids, m_numStates, m_ghostVect);

  Interval interval0(0, m_numStates-1);
  Interval interval1(m_numStates, 2*m_numStates-1);
  inData.copyTo(interval0, m_Unew, interval0);
  inData.copyTo(interval1, m_Uold, interval0);

  // Set up data structures
  levelSetup();
}

//////////////////////////////////////////////////////////////////////////////

// Write plotfile header
void AMRLevelMappedCons::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = 2*m_numStates + 2; // or 3*m_numStates+2 if Ucen written out

  // Set up the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
    { // <UJ> averaged over cell
      sprintf(compStr,"component_%d",comp);
      char stateNameChars[60];
      sprintf(stateNameChars, "%sJ", m_stateNames[comp].c_str());
      header.m_string[compStr] = stateNameChars;
    }
  for (int comp = 0; comp < m_numStates; ++comp)
    { // <U> averaged over cell
      sprintf(compStr,"component_%d",comp + m_numStates);
      header.m_string[compStr] = m_stateNames[comp];
    }
  //  for (int comp = 0; comp < m_numStates; ++comp)
  //    { // U at cell center
  //      sprintf(compStr,"component_%d",comp + 2*m_numStates);
  //      char stateNameChars[60];
  //      sprintf(stateNameChars, "%sCen", m_stateNames[comp].c_str());
  //      header.m_string[compStr] = stateNameChars;
  //    }
  { // J
    int comp = 2*m_numStates; // or 3*m_numStates if including Ucen
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = "J";

    // volume
    comp++;
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = "volume";
  }

  // Write the header
  header.writeToFile(a_handle);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  m_gdnvPhysics->expressions(expressions);
  expressions.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

// Write plotfile data for this level
void AMRLevelMappedCons::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::writePlotLevel" << endl;
  }

  // if we're on level 0, write out the mapped-grid geometry info
  writeMappedPlotFile();

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

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  /////////////////////////////////////////////////////////////////////////
  // Begin application-dependent code - PC.
  /////////////////////////////////////////////////////////////////////////

  // Create state intervals
  Interval baseSrcInterval(0,m_numStates-1);

  int count = 0;
  Interval phiJDstInterval(count,count+m_numStates-1);
  count += m_numStates;

  Interval phiDstInterval(count,count+m_numStates-1);
  count += m_numStates;

  //  Interval phiCenInterval(count,count+m_numStates-1);
  //  count += m_numStates;

  Interval jDstInterval(count,count);
  count += 1;

  Interval volDstInterval(count,count);
  count += 1;

  //  Interval velDstInterval(count,count+SpaceDim-1);
  //  count += SpaceDim;
  //
  //  Interval locDstInterval(count,count+SpaceDim-1);
  //  count += SpaceDim;
  LevelData<FArrayBox>& UnewCopy = (LevelData<FArrayBox>&) m_Unew;
  UnewCopy.exchange();
  // Set up arguments to LevelGodunov::step based on whether there are
  // coarser and finer levels
  const DisjointBoxLayout& layout = m_Unew.getBoxes();
  LevelData<FArrayBox> phiDummy(layout, m_numStates, m_ghostVect);
  // Make copy of Unew, which is <UJ>.
  LevelData<FArrayBox> outData(layout, count);
  // from baseSrcInterval to phiJDstInterval:  <UJ>
  m_Unew.copyTo(baseSrcInterval, outData, phiJDstInterval);

  // from baseSrcInterval to baseSrcInterval
  //  m_Unew.copyTo(baseSrcInterval, phiDummy, baseSrcInterval);
  const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  // phiDummy = m_Unew / cellAvgJ
  m_levelConsOperator.cellUJToCellU(phiDummy, m_Unew, cellAvgJ);
  // phiDstInterval to hold <U>
  phiDummy.copyTo(baseSrcInterval, outData, phiDstInterval);

  //  // Now fill phiDummy with cell-centered U.
  //  fourthOrderAverage(phiDummy, m_problem_domain, -1);
  //  // phiCenInterval to hold U at cell centers
  //  phiDummy.copyTo(baseSrcInterval, outData, phiCenInterval);

  // Jacobian LevelData here!
  cellAvgJ.copyTo(cellAvgJ.interval(), outData, jDstInterval);

  // Volume LevelData here!
  const LevelData<FArrayBox>& volume = m_coordSysPtr->getCellVolumes();
  volume.copyTo(volume.interval(), outData, volDstInterval);

  // write the BoxLayout and the data
  write(a_handle, layout);
  write(a_handle, outData, "data");
}


//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::writeMappedPlotFile() const
{
  // only do this on level 0
  if (m_level == 0)
    {
      // gather AMR Levels and create node-centered dataset of
      // node locations of mapped grids

      Vector<AMRLevel*> vectAMRLevels;
      {
        // cast away const for this to call this function
        AMRLevelMappedCons* nonConstThis = const_cast<AMRLevelMappedCons*>(this);

        vectAMRLevels = nonConstThis->getAMRLevelHierarchy();
      }

      int numLevels = vectAMRLevels.size();

      Vector<int> vectRefRatio(numLevels,0);
      Vector<DisjointBoxLayout> vectGrids(numLevels);
      Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(numLevels, NULL);

      const AMRLevelMappedCons* levelPtr = this;
      // loop over levels and set things up
      for (int lev=0; lev<numLevels; lev++)
        {
          const DisjointBoxLayout& levelGrids = levelPtr->m_grids;
          vectGrids[lev] = levelPtr->m_grids;
          vectRefRatio[lev] = levelPtr->m_ref_ratio;

          vectNodeLoc[lev] = new LevelData<NodeFArrayBox>(levelGrids,
                                                          SpaceDim);

          const CoordSys<FArrayBox,FluxBox>* levelCS = levelPtr->m_coordSysPtr;
          Real levelDx = levelPtr->m_dx;

          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              NodeFArrayBox& thisNodeFAB = (*vectNodeLoc[lev])[dit];
              // node-centered FAB
              FArrayBox& thisFAB = thisNodeFAB.getFab();
              BoxIterator bit(thisFAB.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  // location in index space
                  RealVect nodeIndexLoc(iv);
                  nodeIndexLoc *= levelDx;
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
              m_plotfile_prefix.c_str(), m_levelStep, SpaceDim);

      string fileName(iter_str);

      // now call nodal WriteAMRHierarchy function...
      WriteAMRHierarchyHDF5(fileName,
                            vectGrids,
                            vectNodeLoc,
                            locationNames,
                            m_problem_domain.domainBox(),
                            m_dx,
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

//////////////////////////////////////////////////////////////////////////////

// Returns the dt computed earlier for this level
Real AMRLevelMappedCons::computeDt()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::computeDt " << m_level << endl;
    }

  return m_dtNew;
}

//////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelMappedCons::computeInitialDt()
{
  //  Real newDT = m_initial_dt_multiplier * m_dx / getMaxWaveSpeed(m_Unew);
  Real newDT = m_initial_dt_multiplier * computeNewDt();

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::computeInitialDt on level " << m_level << " = " << newDT << endl;
    }

  return newDT;
}

//////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelMappedCons::computeNewDt()
{
  Real newDT;
  if (m_dtFromCells)
    { // the old method
      newDT = m_dx / getMaxWaveSpeed(m_Unew);
    }
  else
    {
      LevelData<FArrayBox> cellAvgU(m_grids, m_numStates, m_ghostVect);
      const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
      // cellAvgU = m_Unew / cellAvgJ
      m_levelConsOperator.cellUJToCellU(cellAvgU, m_Unew, cellAvgJ);

      LevelData<FArrayBox> cellAvgW(m_grids, m_numStates, m_ghostVect);
      DataIterator dit = m_grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& UcellFab = cellAvgU[dit];
          FArrayBox& WcellFab = cellAvgW[dit];
          Box bx = WcellFab.box();
          m_gdnvPhysics->consToPrim(WcellFab, UcellFab, bx);
        }

      // cellVel:  cell-averaged velocities
      LevelData<FArrayBox> cellVel;
      Interval velInt = m_gdnvPhysics->velocityInterval();
      aliasLevelData(cellVel, &cellAvgW, velInt);

      // cellWaveVel:  velocity +/- speed of sound
      LevelData<FArrayBox> cellWaveVel(m_grids, SpaceDim, m_ghostVect);
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& UcellFab = cellAvgU[dit];
          Box bx = UcellFab.box();
          FArrayBox soundSpeedFab(bx, 1);
          m_gdnvPhysics->soundSpeed(soundSpeedFab, UcellFab, bx);
          FArrayBox& cellVelFab = cellVel[dit];
          FArrayBox& cellWaveVelFab = cellWaveVel[dit];
          FORT_WAVEVELOCITY(CHF_FRA(cellWaveVelFab),
                            CHF_CONST_FRA(cellVelFab),
                            CHF_CONST_FRA1(soundSpeedFab, 0),
                            CHF_BOX(bx));
        }

      // faceVel:  face-averaged wave velocities
      LevelData<FluxBox> faceWaveVel(m_grids, SpaceDim, m_ghostVect);
      // second-order is sufficient
      CellToEdge(cellWaveVel, faceWaveVel);
      // seems I need m_dx factor here
      newDT = computeMappedDt(faceWaveVel, m_coordSysPtr, m_cfl);
    }
  return newDT;
}

  //////////////////////////////////////////////////////////////////////////////

Real AMRLevelMappedCons::getMaxWaveSpeed(const LevelData<FArrayBox>& a_U)
{
  // FOR MAPPED, THIS NEEDS TO BE REWRITTEN, I THINK.
  CH_TIME("AMRLevelMappedCons::getMaxWaveSpeed");
  const DisjointBoxLayout& layout = a_U.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();
  // Initial maximum wave speed
  Real speed = 0.0;
  // Loop over all grids to get the maximum wave speed
  for (dit.begin(); dit.ok(); ++dit)
    {
      // Get maximum wave speed for this grid
      const Box& bx = layout[dit];
      const FArrayBox& Ufab = a_U[dit];
      Real speedOverBox = m_gdnvPhysics->getMaxWaveSpeed(Ufab, bx);
      // Compute a running maximum
      speed = Max(speed, speedOverBox);
    }
  // Gather maximum wave speeds and broadcast the maximum over these
  Vector<Real> allSpeeds;
  gather(allSpeeds, speed, uniqueProc(SerialTask::compute));
  if (procID() == uniqueProc(SerialTask::compute))
    {
      speed = allSpeeds[0];
      for (int i = 1; i < allSpeeds.size (); ++i)
        {
          speed = Max(speed, allSpeeds[i]);
        }
    }
  broadcast(speed, uniqueProc(SerialTask::compute));
  // Return the maximum wave speed
  return speed;
}


//////////////////////////////////////////////////////////////////////////////

// Set the CFL number
void AMRLevelMappedCons::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

//////////////////////////////////////////////////////////////////////////////

// Set the spatial order of accuracy
void AMRLevelMappedCons::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}

//////////////////////////////////////////////////////////////////////////////

// sets whether to limit face values in levelOperator
void AMRLevelMappedCons::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::initialAverage(bool a_initialAverage)
{
  m_initialAverage = a_initialAverage;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
}

//////////////////////////////////////////////////////////////////////////////
void AMRLevelMappedCons::useArtVisc(bool a_useArtVisc)
{
  m_useArtVisc = a_useArtVisc;
}

//////////////////////////////////////////////////////////////////////////////
void AMRLevelMappedCons::ratioArtVisc(Real a_ratioArtVisc)
{
  m_ratioArtVisc = a_ratioArtVisc;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::forwardEuler(bool a_forwardEuler)
{
  m_forwardEuler = a_forwardEuler;
}

//////////////////////////////////////////////////////////////////////////////

// sets whether to enforce a min value for advected quantity
void AMRLevelMappedCons::enforceMinVal(bool a_enforceMinVal, Real a_minVal)
{
  m_enforceMinVal = a_enforceMinVal;
  m_minVal = a_minVal;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

// Set the physical dimension of the longest side of the domain
void AMRLevelMappedCons::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
}

//////////////////////////////////////////////////////////////////////////////

// Set the refinement threshold
void AMRLevelMappedCons::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}

//////////////////////////////////////////////////////////////////////////////

// Set the tag buffer size
void AMRLevelMappedCons::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::godunovPhysics(const GodunovPhysics* const a_gdnvPhysics)
{
  m_gdnvPhysics = a_gdnvPhysics->new_godunovPhysics();
}

////////////////////////////////////////////////////////////////////////////////
void AMRLevelMappedCons::coordinateSystem(CoordSysFactory<FArrayBox,FluxBox>* a_coordSysFact)
{
  m_coordSysFactPtr = a_coordSysFact;
}

//////////////////////////////////////////////////////////////////////////////

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelMappedCons::loadBalance(const Vector<Box>& a_grids)
{
  // Load balance and create boxlayout
  Vector<int> procMap;
  //if (procID() == uniqueProc(SerialTask::compute))
  //{
  //  LoadBalance(procMap,a_grids);
  //}
  //broadcast(procMap,uniqueProc(SerialTask::compute));

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelMappedCons::loadBalance: procesor map: " << endl;
    for (int igrid = 0; igrid < a_grids.size(); ++igrid)
    {
      pout() << igrid << ": " << procMap[igrid] << "  " << a_grids[igrid].volume() << endl;
    }
    pout() << endl;
  }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}

//////////////////////////////////////////////////////////////////////////////

// Setup menagerie of data structures
void AMRLevelMappedCons::levelSetup()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::levelSetup " << m_level << endl;
  }

  AMRLevelMappedCons* amrConsCoarserPtr = getCoarserLevel();
  AMRLevelMappedCons* amrConsFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrConsCoarserPtr != NULL);
  m_hasFiner   = (amrConsFinerPtr   != NULL);

  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();

      m_coarseAverage.define(m_grids, m_numStates, nRefCrse);

      m_fineInterp.define(m_grids, m_numStates, nRefCrse, m_problem_domain);

      const DisjointBoxLayout& coarserLevelLayout = amrConsCoarserPtr->m_grids;

      // Maintain levelConsOperator
      m_levelConsOperator.define(m_grids, coarserLevelLayout,
                                 m_problem_domain,
                                 nRefCrse,
                                 m_dx,
                                 m_gdnvPhysics,
                                 m_numStates,
                                 m_hasCoarser, m_hasFiner);

      // This may look twisted but you have to do this this way because the
      // coarser levels get setup before the finer levels so, since a flux
      // register lives between this level and the next FINER level, the finer
      // level has to do the setup because it is the only one with the
      // information at the time of construction.

      // Maintain flux registers
      amrConsCoarserPtr->m_fluxRegister.define(m_grids,
                                               amrConsCoarserPtr->m_grids,
                                               m_problem_domain,
                                               amrConsCoarserPtr->m_ref_ratio,
                                               m_numStates);
      amrConsCoarserPtr->m_fluxRegister.setToZero();
    }
  else
    {
      m_levelConsOperator.define(m_grids, DisjointBoxLayout(),
                                 m_problem_domain,
                                 m_ref_ratio,
                                 m_dx,
                                 m_gdnvPhysics,
                                 m_numStates,
                                 m_hasCoarser, m_hasFiner);
    }
  // this should happen whether or not there's a coarser level
  m_levelConsOperator.spaceOrder(m_spaceOrder);
  m_levelConsOperator.limitFaceValues(m_limitFaceValues);
  m_levelConsOperator.useFlattening(m_useFlattening);
  m_levelConsOperator.noPPM(m_noPPM);
  m_levelConsOperator.doDeconvolution(m_doDeconvolution);
  m_levelConsOperator.doFaceDeconvolution(m_doFaceDeconvolution);
  m_levelConsOperator.useArtificialViscosity(m_useArtificialViscosity);
  m_levelConsOperator.artificialViscosity(m_artificialViscosity);

  // Create the coordinate system
  m_coordSysPtr = (FourthOrderCoordSys*)
    m_coordSysFactPtr->getCoordSys(m_grids,
                                   m_problem_domain,
                                   m_ghostVect);
}

//////////////////////////////////////////////////////////////////////////////

// Get the next coarser level
AMRLevelMappedCons* AMRLevelMappedCons::getCoarserLevel() const
{
  AMRLevelMappedCons* amrConsCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrConsCoarserPtr = dynamic_cast<AMRLevelMappedCons*>(m_coarser_level_ptr);

    if (amrConsCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelMappedCons::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrConsCoarserPtr;
}

//////////////////////////////////////////////////////////////////////////////

// Get the next finer level
AMRLevelMappedCons* AMRLevelMappedCons::getFinerLevel() const
{
  AMRLevelMappedCons* amrConsFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrConsFinerPtr = dynamic_cast<AMRLevelMappedCons*>(m_finer_level_ptr);

    if (amrConsFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelMappedCons::getFinerLevel: dynamic cast failed");
    }
  }

  return amrConsFinerPtr;
}
