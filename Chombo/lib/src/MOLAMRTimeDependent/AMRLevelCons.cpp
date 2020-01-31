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

#include "AMRLevelCons.H"
#include "FArrayBox.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "computeSum.H"
#include "FourthOrderFillPatch.H"
#include "AMRIO.H"
#include "TimeInterpolatorRK4.H"
#include "LevelRK4.H"
#include "AdvectOpF_F.H"
#include "PolytropicPhysicsF_F.H"
#include "VorticityF_F.H"
#include "AMRLevel.H"
#include "FourthOrderUtil.H"
#include "PhysIBC.H"
// #include "DebugInclude.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////

// Constructor
AMRLevelCons::AMRLevelCons()
{
  CH_TIME("AMRLevelCons::AMRLevelCons");
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelCons default constructor" << endl;
  }
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////

// Destructor
AMRLevelCons::~AMRLevelCons()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons destructor" << endl;
    }

  delete m_levelConsOperatorPtr;
  delete m_molPhysics;
  delete m_coarseAveragePtr;
}

//////////////////////////////////////////////////////////////////////////////

// Define new AMR level
void AMRLevelCons::define(AMRLevel*            a_coarserLevelPtr,
                          const ProblemDomain& a_problemDomain,
                          int                  a_level,
                          int                  a_refRatio)
{
  CH_TIME("AMRLevelCons::define");
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::define " << a_level << endl;
    }

  // Create our operator
  m_levelConsOperatorPtr = new LevelConsOperator();

  // Create our averager
  m_coarseAveragePtr = new CoarseAverage();


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
    AMRLevelCons* amrConsPtr = dynamic_cast<AMRLevelCons*>(a_coarserLevelPtr);

    if (amrConsPtr != NULL)
      {
        transferSettingsFromAMRLevel(amrConsPtr);
      }
    else
      {
        MayDay::Error("AMRLevelCons::define: a_coarserLevelPtr is not castable to AMRLevelCons*");
      }
  }

  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more

  // Begin application-dependent code - dfm
  m_numGhost = 5;
  m_ghostVect = m_numGhost * IntVect::Unit;

  // MOLPhysics* m_molPhysics is defined in the executable where
  // it is set to PolytropicPhysics or whatever.
  m_molPhysics->define(m_problem_domain, m_dx);
  // Number and names of conserved states
  m_numStates  = m_molPhysics->numConserved();
  m_stateNames = m_molPhysics->stateNames();
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

void
AMRLevelCons::setDefaultValues()
{
  m_cfl = 0.8;
  m_spaceOrder = 4;
  m_limitFaceValues = true;
  m_highOrderLimiter = true;
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
  m_refinementIsScaled = false;
  m_tagPressure = false;
  m_tagVorticity = false;
  m_initial_dt_multiplier = 0.1;
  m_molPhysics = NULL;
  m_coarseAveragePtr = NULL;
}

//////////////////////////////////////////////////////////////////////////////

// Advance by one timestep
Real AMRLevelCons::advance()
{
  CH_TIME("AMRLevelCons::advance");
  advanceU();
  Real returnDt = computeNewDt();
  return returnDt;
}


//////////////////////////////////////////////////////////////////////////////

// Advance by one timestep
void AMRLevelCons::advanceU()
{
  CH_TIME("AMRLevelCons::advanceU");
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelCons::advance level " << m_level
           << " to time " << m_time << endl;
  }

  // Begin application-dependent code - dfm

  // Copy the new to the old
  //  m_Unew.copyTo(m_Unew.interval(),
  //                m_Uold,
  //                m_Uold.interval());
  // petermc, 13 Jan 2009:  Copy ghost cells too.
  m_levelConsOperatorPtr->copySolnData(m_Uold, m_Unew);

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
      AMRLevelCons* coarserPtr = getCoarserLevel();

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

  // Set m_levelConsOperatorPtr->m_fluxes = 0.
  m_levelConsOperatorPtr->resetFluxes();

  // We need 1 layer of ghost cells of UoldCopy.
  LevelData<FArrayBox> UoldCopy(m_grids, m_numStates, m_ghostVect);
  m_levelConsOperatorPtr->copySolnData(UoldCopy, m_Uold);
  Real oldTime = m_time;
  if (m_forwardEuler)
    {
      LevelData<FArrayBox> RHSTmp;
      m_levelConsOperatorPtr->defineRHSData(RHSTmp, m_Uold);
      // start by copying old solution into new solution
      m_levelConsOperatorPtr->copySolnData(m_Unew, m_Uold);
      m_levelConsOperatorPtr->evalCountMax(1); // to make 1 call to evalRHS
      m_levelConsOperatorPtr->resetEvalCount();
      m_levelConsOperatorPtr->evalRHS(RHSTmp, m_Uold,
                                  *finerFR, *coarserFR,
                                  *coarserUold, tCoarserOld,
                                  *coarserUnew, tCoarserNew,
                                  m_time, m_dt);
      m_levelConsOperatorPtr->updateODE(m_Unew, RHSTmp, m_dt);
    }
  else
    {
      m_levelConsOperatorPtr->evalCountMax(4); // to make 4 calls to evalRHS
      m_levelConsOperatorPtr->resetEvalCount();
      if (m_hasFiner)
        {
          // Fill in the interpolator from this level to next finer level:
          // we'll use it later when we calculate on next finer level.
          TimeInterpolatorRK4& timeInterpolator =
            getFinerLevel()->m_levelConsOperatorPtr->getTimeInterpolator();

          RK4LevelAdvance<LevelData<FArrayBox>,
            TimeInterpolatorRK4,
            LevelFluxRegister,
            LevelConsOperator>(m_Unew, m_Uold, timeInterpolator,
                               *coarserUold, tCoarserOld,
                               *coarserUnew, tCoarserNew,
                               *coarserFR, *finerFR,
                               m_time, m_dt, *m_levelConsOperatorPtr);
        }
      else
        {
          RK4LevelAdvance<LevelData<FArrayBox>,
            LevelFluxRegister,
            LevelConsOperator>(m_Unew, m_Uold,
                               *coarserUold, tCoarserOld,
                               *coarserUnew, tCoarserNew,
                               *coarserFR, *finerFR,
                               m_time, m_dt, *m_levelConsOperatorPtr);
        }
    }

  if (m_useArtificialViscosity)
    {
      // We need 1 layer of ghost cells of UoldCopy.
      m_levelConsOperatorPtr->fillGhosts(UoldCopy, oldTime,
                                         tCoarserOld, tCoarserNew);

      // find flux, fluxArtificialViscosity, from artificial viscosity;
      // then m_fluxes += fluxArtificialViscosity;
      // and m_Unew -= div(fluxArtificialViscosity) * m_dt.
      m_levelConsOperatorPtr->addArtificialViscosity(m_Unew, UoldCopy,
                                                 *finerFR, *coarserFR,
                                                 oldTime, m_dt);
    }

  if (m_useArtVisc)
    {
      addDissipation();
    }

  // Update the time and store the new timestep
  m_time += m_dt;
}


//////////////////////////////////////////////////////////////////////////////

Real AMRLevelCons::computeNewDt()
{
  CH_TIME("AMRLevelCons::computeNewDt");
  Real returnDt = m_dx / getMaxWaveSpeed(m_Unew);
  returnDt *= m_cfl;
  m_dtNew = returnDt;
  return returnDt;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::addDissipation()
{
  CH_TIME("AMRLevelCons::addDissipation");
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
      for (int idir=0; idir<SpaceDim; idir++)
        if (!m_problem_domain.isPeriodic(idir))
          {
            if (boxBase.smallEnd(idir) == domainLo[idir])
              { // shrink boxLapU and boxLapLapU on the low end
                boxLapU.growLo(idir, -1);
                boxLapLapU.growLo(idir, -2);
              }
            if (boxBase.bigEnd(idir) == domainHi[idir])
              { // shrink boxLapU and boxLapLapU on the high end
                boxLapU.growHi(idir, -1);
                boxLapLapU.growHi(idir, -2);
              }
          }
      FArrayBox lapUfab(boxLapU, numU);
      lapUfab.setVal(0.);
      const FArrayBox& UoldFab = m_Uold[dit];
      for (int idir=0; idir<SpaceDim; idir++)
        {
          FORT_CCLAPLACIAN(CHF_FRA(lapUfab),
                           CHF_CONST_FRA(UoldFab),
                           CHF_BOX(boxLapU),
                           CHF_CONST_INT(idir),
                           CHF_CONST_REAL(m_dx));
        }
      FArrayBox lapLapUfab(boxLapLapU, numU);
      lapLapUfab.setVal(0.);
      for (int idir=0; idir<SpaceDim; idir++)
        {
          FORT_CCLAPLACIAN(CHF_FRA(lapLapUfab),
                           CHF_CONST_FRA(lapUfab),
                           CHF_BOX(boxLapLapU),
                           CHF_CONST_INT(idir),
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
void AMRLevelCons::postTimeStep()
{
  CH_TIME("AMRLevelCons::postTimeStep");
  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::postTimeStep " << m_level << endl;
    }

  // Begin application-dependent code - PC.

  if (m_hasFiner)
    {
      // Reflux
      Real scale = -1.0/m_dx;
      m_fluxRegister.reflux(m_Unew, scale);

      // Average from finer level data
      AMRLevelCons* amrConsFinerPtr = getFinerLevel();

      // Set m_Unew at cell c to the mean of amrConsFinerPtr->m_Unew
      // over all finer cells within c.
      amrConsFinerPtr->
        m_coarseAveragePtr->averageToCoarse(m_Unew, amrConsFinerPtr->m_Unew);
    }
  // End application-dependent code - PC.

  if (s_verbosity >= 2 && m_level == 0)
    {
      int nRefFine = 1;

      pout() << "AMRLevelCons::postTimeStep:" << endl;
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
      pout() << "AMRLevelCons::postTimeStep " << m_level
             << " finished" << endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

// Create tags for regridding
void AMRLevelCons::tagCells(IntVectSet& a_tags)
{
  CH_TIME("AMRLevelCons::tagCells");
  // CH_assert(allDefined());
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::tagCells " << m_level << endl;
    }

  // Begin application-dependent code - PC.

  // Create tags based on undivided gradient of density
  const DisjointBoxLayout& levelDomain = m_Unew.disjointBoxLayout();
  IntVectSet localTags;

  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
  {
    const AMRLevelCons* amrLevelCoarserPtr = getCoarserLevel();
    const LevelData<FArrayBox>& UnewCoarser = amrLevelCoarserPtr->m_Unew;
    FourthOrderFillPatch filler(levelDomain,
                                UnewCoarser.disjointBoxLayout(),
                                m_numStates,
                                amrLevelCoarserPtr->m_problem_domain,
                                m_ref_ratio,
                                m_numGhost,
                                true); // fixed time
    // No interpolation in time.
    filler.fillInterp(m_Unew, UnewCoarser, 0, 0, m_numStates);
  }
  // m_Unew.exchange();
  m_levelConsOperatorPtr->exchangeGhosts(m_Unew);

  int numW = m_molPhysics->numPrimitives();
  // index into primitive variables
  int tagIndex = indexForTagging();

  // for vorticity, use m_refineThresh = 20.
  Real threshold = m_refineThresh;
  LevelData<FArrayBox> vecMag(levelDomain, 1);
  if (m_tagVorticity)
    {
      // also fills in vecMag
      threshold = maxVorticity(vecMag);
      threshold *= m_refineThresh;
    }
  else
    {
      // petermc, 15 July 2009:  possibility for m_refinementIsScaled
      if (m_refinementIsScaled) threshold *= m_dx;
    }

  // Compute relative gradient
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bx = levelDomain[dit];
      FArrayBox& vecMagFab = vecMag[dit];
      if (!m_tagVorticity)
        {
          // If tagging on pressure or density, vecFab is gradient.
          int numVecComps = SpaceDim;

          FArrayBox vecFab(bx, numVecComps);
          const FArrayBox& UFab = m_Unew[dit];

          Box bx1 = grow(bx, 1); // need this for gradient
          FArrayBox WFab(bx1, numW);
          m_molPhysics->consToPrim(WFab, UFab, bx1);

          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              const Box bxCenter = bx & grow(m_problem_domain,-BASISV(idir));

              const Box bxLo     = bx & adjCellLo(bxCenter,idir);
              const int hasLo = ! bxLo.isEmpty();

              const Box bxHi     = bx & adjCellHi(bxCenter,idir);
              const int hasHi = ! bxHi.isEmpty();

              FORT_GETRELGRADF(CHF_FRA1(vecFab, idir),
                               CHF_CONST_FRA1(WFab, tagIndex),
                               CHF_CONST_INT(idir),
                               CHF_BOX(bxLo),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(bxHi),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(bxCenter));
            }
          FORT_MAGNITUDEF(CHF_FRA1(vecMagFab,0),
                          CHF_CONST_FRA(vecFab),
                          CHF_BOX(bx));
        }
      // Tag where vector magnitude exceeds threshold
      BoxIterator bit(bx);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          if (vecMagFab(iv) >= threshold)
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
void AMRLevelCons::tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

//////////////////////////////////////////////////////////////////////////////

int AMRLevelCons::indexForTagging()
{
  int ind = (m_tagPressure) ?
    m_molPhysics->pressureIndex() :
    m_molPhysics->densityIndex();
  return ind;
}

//////////////////////////////////////////////////////////////////////////////

Real AMRLevelCons::maxVorticity(LevelData<FArrayBox>&  a_vecMag)
{
  CH_TIME("AMRLevelCons::maxVorticity");
  int numW = m_molPhysics->numPrimitives();
  const Interval& velInt = m_molPhysics->velocityInterval();
  int numVecComps = SpaceDim;
  if (SpaceDim == 2) numVecComps = 1;
  Real vorticityMaxLocal = 0.;
  const DisjointBoxLayout& levelDomain = m_Unew.disjointBoxLayout();
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bx = levelDomain[dit];
      FArrayBox vecFab(bx, numVecComps);
      const FArrayBox& UFab = m_Unew[dit];

      Box bx1 = grow(bx, 1); // need this for gradient
      FArrayBox WFab(bx1, numW);
      FArrayBox velFab(velInt, WFab);
      m_molPhysics->consToPrim(WFab, UFab, bx1);
      // For tagging purposes, we don't need highly accurate
      // estimate of vorticity.
      if (SpaceDim == 3)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              FORT_COMPUTEVORT(CHF_FRA1(vecFab, idir),
                               CHF_CONST_FRA(velFab),
                               CHF_BOX(bx),
                               CHF_CONST_REAL(m_dx),
                               CHF_CONST_INT(idir));
            }
        }
      else if (SpaceDim == 2)
        {
          int idir = 0; // ignored by FORT_COMPUTEVORT
          FORT_COMPUTEVORT(CHF_FRA1(vecFab, idir),
                           CHF_CONST_FRA(velFab),
                           CHF_BOX(bx),
                           CHF_CONST_REAL(m_dx),
                           CHF_CONST_INT(idir));
        }
      else
        {
          MayDay::Error("Vorticity requires 2D or 3D");
        }

      FArrayBox& vecMagFab = a_vecMag[dit];
      FORT_MAGNITUDEF(CHF_FRA1(vecMagFab, 0),
                      CHF_CONST_FRA(vecFab),
                      CHF_BOX(bx));
      Real fabMax = vecMagFab.max();
      if (fabMax > vorticityMaxLocal) vorticityMaxLocal = fabMax;
    }
  Real vorticityMax = vorticityMaxLocal;
#ifdef CH_MPI
  MPI_Allreduce(&vorticityMaxLocal, &vorticityMax, 1,
                MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
#endif
  return vorticityMax;
}

//////////////////////////////////////////////////////////////////////////////

// Set up data on this level after regridding
void AMRLevelCons::regrid(const Vector<Box>& a_newGrids)
{
  CH_TIME("AMRLevelCons::regrid");
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::regrid " << m_level << endl;
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

  // Interpolate from coarser level
  if (m_hasCoarser)
    {
      AMRLevelCons* amrConsCoarserPtr = getCoarserLevel();
      m_fineInterp.interpToFine(m_Unew, amrConsCoarserPtr->m_Unew);

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
void AMRLevelCons::initialGrid(const Vector<Box>& a_newGrids)
{
  CH_TIME("AMRLevelCons::initialGrid");
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelCons::initialGrid " << m_level << endl;
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

  // m_Unew.define(m_grids, m_numStates, m_ghostVect);
  // m_Uold.define(m_grids, m_numStates, m_ghostVect);
  defineDataHolders();

  // End application-dependent code - PC.

  // Set up data structures
  levelSetup();
}

//////////////////////////////////////////////////////////////////////////////
void AMRLevelCons::defineDataHolders()
{
  CH_TIME("AMRLevelCons::defineDataHolders");
  m_Unew.define(m_grids, m_numStates, m_ghostVect);
  m_Uold.define(m_grids, m_numStates, m_ghostVect);
}


//////////////////////////////////////////////////////////////////////////////

// Initialize data:  This is a function of AMRLevel
void AMRLevelCons::initialData()
{
  CH_TIME("AMRLevelCons::initialData");
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::initialData " << m_level << endl;
    }
  setData(m_Unew, 0.);
}

//////////////////////////////////////////////////////////////////////////////

// Initialize data
void AMRLevelCons::setData(LevelData<FArrayBox>& a_U,
                           Real a_time) const
{
  CH_TIME("AMRLevelCons::setData");
  // Begin application-dependent code - PC.

  // FIX:  BasicAdvectIBC requires setTime.
  // m_basicIBCPtr->setTime(0.0);
  PhysIBC* physIBCPtr = m_molPhysics->getPhysIBC();
  physIBCPtr->initialize(a_U);

  //  m_Unew.exchange();
  if (m_initialAverage)
    { // call to new function, petermc, 19 Dec 2008
      fourthOrderAverage(a_U, m_problem_domain);
    }
  // m_Unew.exchange();
  m_levelConsOperatorPtr->exchangeGhosts(a_U);
  // End application-dependent code - PC.
}

//////////////////////////////////////////////////////////////////////////////

// Things to do after initialization
void AMRLevelCons::postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::postInitialize " << m_level << endl;
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
void AMRLevelCons::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  CH_TIME("AMRLevelCons::writeCheckpointHeader");
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::writeCheckpointHeader" << endl;
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
void AMRLevelCons::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  CH_TIME("AMRLevelCons::writeCheckpointLevel");
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::writeCheckpointLevel" << endl;
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
  D_TERM6(
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
            header.m_int ["is_periodic_2"] = 0; ,

          if (m_problem_domain.isPeriodic(3))
            header.m_int ["is_periodic_3"] = 1;
          else
            header.m_int ["is_periodic_3"] = 0; ,

          if (m_problem_domain.isPeriodic(4))
            header.m_int ["is_periodic_4"] = 1;
          else
            header.m_int ["is_periodic_4"] = 0; ,

          if (m_problem_domain.isPeriodic(5))
            header.m_int ["is_periodic_5"] = 1;
          else
            header.m_int ["is_periodic_5"] = 0; );

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
void AMRLevelCons::readCheckpointHeader(HDF5Handle& a_handle)
{
  CH_TIME("AMRLevelCons::readCheckpointHeader");
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::readCheckpointHeader" << endl;
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
      MayDay::Error("AMRLevelCons::readCheckpointHeader: checkpoint file does not have num_components");
    }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates*2)
    {
      MayDay::Error("AMRLevelCons::readCheckpointHeader: num_components in checkpoint file does not match solver");
    }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates*2; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      if (header.m_string.find(compStr) == header.m_string.end())
        {
          MayDay::Error("AMRLevelCons::readCheckpointHeader: checkpoint file does not have enough component names");
        }

      stateName = header.m_string[compStr];
      if (stateName != m_stateNames[comp])
        {
          MayDay::Error("AMRLevelCons::readCheckpointHeader: state_name in checkpoint does not match solver");
        }
    }
}


//////////////////////////////////////////////////////////////////////////////

// Read checkpoint data for this level
void AMRLevelCons::readCheckpointLevel(HDF5Handle& a_handle)
{
  CH_TIME("AMRLevelCons::readCheckpointLevel");
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::readCheckpointLevel" << endl;
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
    MayDay::Error("AMRLevelCons::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelCons::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelCons::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("AMRLevelCons::readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
    {
      pout() << "read dt = " << m_dt << endl;
    }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMRLevelCons::readCheckpointLevel: file does not contain time");
    }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
    {
      pout() << "read time = " << m_time << endl;
    }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
    {
      MayDay::Error("AMRLevelCons::readCheckpointLevel: file does not contain prob_domain");
    }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM6(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
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
            isPeriodic[2] = false; ,

          if (!(header.m_int.find("is_periodic_3") == header.m_int.end()))
            isPeriodic[3] =  (header.m_int["is_periodic_3"] == 1);
          else
            isPeriodic[3] = false; ,

          if (!(header.m_int.find("is_periodic_4") == header.m_int.end()))
            isPeriodic[4] =  (header.m_int["is_periodic_4"] == 1);
          else
            isPeriodic[4] = false; ,

          if (!(header.m_int.find("is_periodic_5") == header.m_int.end()))
            isPeriodic[5] =  (header.m_int["is_periodic_5"] == 1);
          else
            isPeriodic[5] = false;);

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle, grids);

  if (gridStatus != 0)
    {
      MayDay::Error("AMRLevelCons::readCheckpointLevel: file does not contain a Vector<Box>");
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
      MayDay::Error("AMRLevelCons::readCheckpointLevel: file does not contain all state data");
    }

  // m_Unew.define(m_grids, m_numStates, m_ghostVect);
  // m_Uold.define(m_grids, m_numStates, m_ghostVect);
  defineDataHolders();

  Interval interval0(0, m_numStates-1);
  Interval interval1(m_numStates, 2*m_numStates-1);
  inData.copyTo(interval0, m_Unew, interval0);
  inData.copyTo(interval1, m_Uold, interval0);

  // Set up data structures
  levelSetup();
}

//////////////////////////////////////////////////////////////////////////////

// Write plotfile header
void AMRLevelCons::writePlotHeader(HDF5Handle& a_handle) const
{
  CH_TIME("AMRLevelCons::writePlotHeader");
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelCons::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  int numStatesOutput = m_numStates;
  if (m_tagVorticity) numStatesOutput++;
  header.m_int["num_components"] = numStatesOutput;

  // Set up the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = m_stateNames[comp];
    }
  if (m_tagVorticity)
    {
      int comp = m_numStates;
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = "HOvorticity";
    }

  // Write the header
  header.writeToFile(a_handle);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  m_molPhysics->expressions(expressions);
  expressions.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

// Write plotfile data for this level
void AMRLevelCons::writePlotLevel(HDF5Handle& a_handle) const
{
  CH_TIME("AMRLevelCons::writePlotLevel");
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelCons::writePlotLevel" << endl;
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

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  write(a_handle, m_Unew.boxLayout());
  if (m_tagVorticity)
    {
      // Vorticity requires 2 ghost layers.
      LevelData<FArrayBox> outData(m_Unew.getBoxes(), m_numStates + 1,
                                   m_ghostVect);

      Interval intervalU(0, m_numStates-1);
      LevelData<FArrayBox> UnewCopy;
      aliasLevelData(UnewCopy, &outData, intervalU);
      // m_Unew.copyTo(UnewCopy);
      DataIterator dit = UnewCopy.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) UnewCopy[dit].copy(m_Unew[dit]);
      // If there is a coarser level interpolate undefined ghost cells
      if (m_hasCoarser)
        {
          const AMRLevelCons* amrLevelCoarserPtr = getCoarserLevel();
          const LevelData<FArrayBox>& UnewCoarser = amrLevelCoarserPtr->m_Unew;
          // need 4 layers of ghost cells for vorticity
          int nghostCoarseFill = 4;
          FourthOrderFillPatch filler(m_Unew.disjointBoxLayout(),
                                      UnewCoarser.disjointBoxLayout(),
                                      m_numStates,
                                      amrLevelCoarserPtr->m_problem_domain,
                                      m_ref_ratio,
                                      nghostCoarseFill,
                                      true); // fixed time
          // No interpolation in time.
          filler.fillInterp(UnewCopy, UnewCoarser, 0, 0, m_numStates);
        }
      UnewCopy.exchange();

      Interval intervalVorticity(m_numStates, m_numStates);
      LevelData<FArrayBox> vorticity;
      aliasLevelData(vorticity, &outData, intervalVorticity);
      computeVorticity(vorticity, UnewCopy);
      write(a_handle, outData, "data");
    }
  else
    {
      write(a_handle, m_Unew, "data");
    }
}

#endif

//////////////////////////////////////////////////////////////////////////////

// Returns the dt computed earlier for this level
Real AMRLevelCons::computeDt()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::computeDt " << m_level << endl;
    }

  return m_dtNew;
}

//////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelCons::computeInitialDt()
{
  CH_TIME("AMRLevelCons::computeInitialDt");
  // Real newDT = computeDt();
  // newDT = min(m_initial_dt_multiplier * m_dx, newDT);
  Real newDT = m_initial_dt_multiplier * m_dx / getMaxWaveSpeed(m_Unew);
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelCons::computeInitialDt on level " << m_level << " = " << newDT << endl;
    }

  return newDT;
}

//////////////////////////////////////////////////////////////////////////////

Real AMRLevelCons::getMaxWaveSpeed(const LevelData<FArrayBox>& a_U)
{
  CH_TIME("AMRLevelCons::getMaxWaveSpeed");
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
      Real speedOverBox = m_molPhysics->getMaxWaveSpeed(Ufab, bx);
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
void AMRLevelCons::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

//////////////////////////////////////////////////////////////////////////////

// Set the spatial order of accuracy
void AMRLevelCons::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}

//////////////////////////////////////////////////////////////////////////////

// sets whether to limit face values in levelOperator
void AMRLevelCons::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

//////////////////////////////////////////////////////////////////////////////

// sets whether to use high-order limiter
void AMRLevelCons::highOrderLimiter(bool a_highOrderLimiter)
{
  m_highOrderLimiter = a_highOrderLimiter;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::initialAverage(bool a_initialAverage)
{
  m_initialAverage = a_initialAverage;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
}

//////////////////////////////////////////////////////////////////////////////
void AMRLevelCons::useArtVisc(bool a_useArtVisc)
{
  m_useArtVisc = a_useArtVisc;
}

//////////////////////////////////////////////////////////////////////////////
void AMRLevelCons::ratioArtVisc(Real a_ratioArtVisc)
{
  m_ratioArtVisc = a_ratioArtVisc;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::forwardEuler(bool a_forwardEuler)
{
  m_forwardEuler = a_forwardEuler;
}

//////////////////////////////////////////////////////////////////////////////

// sets whether to enforce a min value for advected quantity
void AMRLevelCons::enforceMinVal(bool a_enforceMinVal, Real a_minVal)
{
  m_enforceMinVal = a_enforceMinVal;
  m_minVal = a_minVal;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

// Set the physical dimension of the longest side of the domain
void AMRLevelCons::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
}

//////////////////////////////////////////////////////////////////////////////

// Set the refinement threshold
void AMRLevelCons::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}

//////////////////////////////////////////////////////////////////////////////

// Set whether refinement threshold is scaled with dx
void AMRLevelCons::refinementIsScaled(bool a_refinementIsScaled)
{
  m_refinementIsScaled = a_refinementIsScaled;
}

//////////////////////////////////////////////////////////////////////////////

// Set whether to tag on pressure instead of on density
void AMRLevelCons::tagPressure(bool a_tagPressure)
{
  m_tagPressure = a_tagPressure;
}

//////////////////////////////////////////////////////////////////////////////

// Set whether to tag on vorticity instead of on density
void AMRLevelCons::tagVorticity(bool a_tagVorticity)
{
  m_tagVorticity = a_tagVorticity;
}

//////////////////////////////////////////////////////////////////////////////

// Set the tag buffer size
void AMRLevelCons::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::molPhysics(const MOLPhysics* const a_molPhysics)
{
  m_molPhysics = a_molPhysics->new_molPhysics();
}

//////////////////////////////////////////////////////////////////////////////

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelCons::loadBalance(const Vector<Box>& a_grids)
{
  CH_TIME("AMRLevelCons::loadBalance");
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
    pout() << "AMRLevelCons::loadBalance: procesor map: " << endl;
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
void AMRLevelCons::levelSetup()
{
  CH_TIME("AMRLevelCons::levelSetup");
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelCons::levelSetup " << m_level << endl;
  }

  AMRLevelCons* amrConsCoarserPtr = getCoarserLevel();
  AMRLevelCons* amrConsFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrConsCoarserPtr != NULL);
  m_hasFiner   = (amrConsFinerPtr   != NULL);

  m_levelConsOperatorPtr->numGhost(m_numGhost);
  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();

      m_coarseAveragePtr->define(m_grids, m_numStates, nRefCrse);

      m_fineInterp.define(m_grids, m_numStates, nRefCrse, m_problem_domain);

      const DisjointBoxLayout& coarserLevelLayout = amrConsCoarserPtr->m_grids;

      // Maintain levelConsOperator
      m_levelConsOperatorPtr->define(m_grids, coarserLevelLayout,
                                     m_problem_domain,
                                     nRefCrse,
                                     m_dx,
                                     m_molPhysics,
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
      m_levelConsOperatorPtr->define(m_grids, DisjointBoxLayout(),
                                 m_problem_domain,
                                 m_ref_ratio,
                                 m_dx,
                                 m_molPhysics,
                                 m_numStates,
                                 m_hasCoarser, m_hasFiner);
    }

  transferSettingsToLevelOp();
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::transferSettingsToLevelOp()
{
  CH_TIME("AMRLevelCons::transferSettingsToLevelOp");
  // this should happen whether or not there's a coarser level
  m_levelConsOperatorPtr->spaceOrder(m_spaceOrder);
  m_levelConsOperatorPtr->limitFaceValues(m_limitFaceValues);
  m_levelConsOperatorPtr->highOrderLimiter(m_highOrderLimiter);
  m_levelConsOperatorPtr->useFlattening(m_useFlattening);
  m_levelConsOperatorPtr->noPPM(m_noPPM);
  m_levelConsOperatorPtr->doDeconvolution(m_doDeconvolution);
  m_levelConsOperatorPtr->doFaceDeconvolution(m_doFaceDeconvolution);
  m_levelConsOperatorPtr->useArtificialViscosity(m_useArtificialViscosity);
  m_levelConsOperatorPtr->artificialViscosity(m_artificialViscosity);
  m_levelConsOperatorPtr->forwardEuler(m_forwardEuler);
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::transferSettingsFromAMRLevel(AMRLevelCons* a_amrConsPtr)
{
  CH_TIME("AMRLevelCons::transferSettingsFromAMRLevel");
  CH_assert(a_amrConsPtr != NULL);
  // These are AMRLevel members, so we have to use get/set
  a_amrConsPtr->initialDtMultiplier(initialDtMultiplier());
  a_amrConsPtr->verbosity(verbosity());

  // These are AMRLevelCons specific, so we can access member variables
  m_cfl = a_amrConsPtr->m_cfl;
  m_spaceOrder = a_amrConsPtr->m_spaceOrder;
  m_limitFaceValues = a_amrConsPtr->m_limitFaceValues;
  m_highOrderLimiter = a_amrConsPtr->m_highOrderLimiter;
  m_useFlattening = a_amrConsPtr->m_useFlattening;
  m_useArtVisc = a_amrConsPtr->m_useArtVisc;
  m_ratioArtVisc = a_amrConsPtr->m_ratioArtVisc;
  m_enforceMinVal = a_amrConsPtr->m_enforceMinVal;
  m_noPPM = a_amrConsPtr->m_noPPM;
  m_doDeconvolution = a_amrConsPtr->m_doDeconvolution;
  m_doFaceDeconvolution = a_amrConsPtr->m_doFaceDeconvolution;
  m_useArtificialViscosity = a_amrConsPtr->m_useArtificialViscosity;
  m_artificialViscosity = a_amrConsPtr->m_artificialViscosity;
  m_forwardEuler = a_amrConsPtr->m_forwardEuler;
  m_minVal = a_amrConsPtr->m_minVal;
  m_domainLength = a_amrConsPtr->m_domainLength;
  m_refineThresh = a_amrConsPtr->m_refineThresh;
  m_refinementIsScaled = a_amrConsPtr->m_refinementIsScaled;
  m_tagPressure = a_amrConsPtr->m_tagPressure;
  m_tagVorticity = a_amrConsPtr->m_tagVorticity;
  m_tagBufferSize = a_amrConsPtr->m_tagBufferSize;
  m_initialAverage = a_amrConsPtr->m_initialAverage;
  // petermc, 22 May 2012: Remove this because it replaces m_molPhysics
  // on coarser level with that on finer level, which is wrong.
  // m_molPhysics = a_amrConsPtr->m_molPhysics;
}

//////////////////////////////////////////////////////////////////////////////

// Get the next coarser level
AMRLevelCons* AMRLevelCons::getCoarserLevel() const
{
  AMRLevelCons* amrConsCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrConsCoarserPtr = dynamic_cast<AMRLevelCons*>(m_coarser_level_ptr);

    if (amrConsCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelCons::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrConsCoarserPtr;
}

//////////////////////////////////////////////////////////////////////////////

// Get the next finer level
AMRLevelCons* AMRLevelCons::getFinerLevel() const
{
  AMRLevelCons* amrConsFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrConsFinerPtr = dynamic_cast<AMRLevelCons*>(m_finer_level_ptr);

    if (amrConsFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelCons::getFinerLevel: dynamic cast failed");
    }
  }

  return amrConsFinerPtr;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelCons::computeVorticity(LevelData<FArrayBox>& a_vorticity,
                                    const LevelData<FArrayBox>& a_U) const
{
  CH_TIME("AMRLevelCons::computeVorticity");
  // These are actually indices into the PRIMITIVE variables,
  // not the conserved variables.
  // I'll need to add functions to MOLPhysics.
  int densityIndex = m_molPhysics->densityIndex();
  Interval densInt(densityIndex, densityIndex);
  const Interval& velInt = m_molPhysics->velocityInterval();

  const int indx = 0;
  const int indy = 1;

  RealVect dxVect = m_dx * RealVect::Unit;

  // Use this alias so that the function remains const.
  //  LevelData<FArrayBox>& UnewCopy = (LevelData<FArrayBox>&) m_Unew;
  DataIterator dit = a_U.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& boxBase = m_grids[dit];
      Box bx2 = grow(boxBase, m_ghostVect - 2*IntVect::Unit);
      Box bx1 = grow(boxBase, m_ghostVect - IntVect::Unit);

      const FArrayBox& UFab = a_U[dit];
      // can't alias const FArrayBox&
      FArrayBox& UFabAlias = (FArrayBox&) UFab;
      // momentumAvgFab = <rho*(u, v)>
      FArrayBox momentumAvgFab(velInt, UFabAlias);
      // densityAvgFab = <rho>
      FArrayBox densityAvgFab(densInt, UFabAlias);

      // Maybe these FABs should be on a bigger box...

      // velocityAvgFab = cell-averaged (u, v)
      FArrayBox velocityAvgFab(bx1, SpaceDim);
      cellFGToCellF(velocityAvgFab, momentumAvgFab, densityAvgFab, bx1);

      FluxBox velocityAvgFlub(bx2, SpaceDim);
      // This uses the coefficients 7/12, -1/12, -1/12, 7/12.
      fourthOrderCellToFace(velocityAvgFlub, velocityAvgFab);

      FArrayBox dudyAvgFab(bx2, 1);
      FORT_SPACEDERIV(CHF_FRA1(dudyAvgFab, 0),
                      CHF_CONST_FRA1(velocityAvgFlub[indy], indx),
                      CHF_BOX(bx2),
                      CHF_CONST_REALVECT(dxVect),
                      CHF_CONST_INT(indy));

      FArrayBox dvdxAvgFab(bx2, 1);
      FORT_SPACEDERIV(CHF_FRA1(dvdxAvgFab, 0),
                      CHF_CONST_FRA1(velocityAvgFlub[indx], indy),
                      CHF_BOX(bx2),
                      CHF_CONST_REALVECT(dxVect),
                      CHF_CONST_INT(indx));

      // vorticityFab = du/dy - dv/dx at cell centers
      FArrayBox& vorticityFab = a_vorticity[dit];
      vorticityFab.copy(dudyAvgFab);
      vorticityFab -= dvdxAvgFab;
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
}

#include "NamespaceFooter.H"
