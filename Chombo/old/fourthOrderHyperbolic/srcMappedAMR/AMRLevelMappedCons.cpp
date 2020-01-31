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
#include <limits>

#include "parstream.H"

#include "AMRLevelMappedCons.H"
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
#include "NodeAMRIO.H"
#include "TimeInterpolatorRK4.H"
#include "LevelRK4.H"
#include "AMRLevel.H"
#include "computeNorm.H"
#include "FourthOrderUtil.H"
#include "PhysMappedIBC.H"
#include "LevelGridMetrics.H"
#include "MultiBlockFluxRegister.H"

#include "NamespaceHeader.H"

// Define this to write errors in <U> in place of <JU>.  Of course, you must
// be able to compute an exact solution in reportError().  It must be defined
// to your last iteration because currently, plotting the error will destroy
// the solution in <JU>
//#define WRITEERRORPLOT 1

/*==============================================================================
 * Definition of static variables
 *============================================================================*/

int AMRLevelMappedCons::s_finestLevel = 0;
Vector<Real> AMRLevelMappedCons::s_L1NormJUConsvRef;


/*******************************************************************************
 *
 * Class AMRLevelMappedCons: member definitions
 *
 ******************************************************************************/

//////////////////////////////////////////////////////////////////////////////

// Factory constructor
AMRLevelMappedCons::AMRLevelMappedCons(MultiBlockCoordSysFactory *const a_coordSysFact,
                                       const AMRLevelMappedStabilityStrategy* const a_stabilityStrategy,
                                       const AMRLevelMappedTaggingStrategy* const a_taggingStrategy,
                                       const std::string&        a_plotPrefix)
  :
  m_coordSysFact(a_coordSysFact),
  m_coordSys(NULL),
  m_fluxRegister(NULL),
  m_stabilityStrategy(a_stabilityStrategy->clone()),
  m_taggingStrategy(a_taggingStrategy->clone()),
  m_plotPrefix(a_plotPrefix)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons default constructor" << endl;
  }
  m_levelConsOperatorPtr = new LevelMappedConsOperator();
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////

// Destructor
AMRLevelMappedCons::~AMRLevelMappedCons()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons destructor" << endl;
    }
  if (m_levelConsOperatorPtr != NULL)
    {
      delete m_levelConsOperatorPtr;
    }
  if (m_molPhysics != NULL)
    {
      delete m_molPhysics;
      m_molPhysics = NULL;
    }
  if (m_levelGridMetricsPtr != NULL)
    {
      delete m_levelGridMetricsPtr;
    }
  if (m_fluxRegister != NULL)
    {
      delete m_fluxRegister;
    }
  delete m_stabilityStrategy;
  delete m_taggingStrategy;
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
  LevelGridMetrics* coarserLGM = NULL;
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelMappedCons* amrConsPtr =
      dynamic_cast<AMRLevelMappedCons*>(a_coarserLevelPtr);
    if (amrConsPtr != NULL)
      {
        transferSettingsFromAMRLevel(amrConsPtr);
        coarserLGM = amrConsPtr->m_levelGridMetricsPtr;
      }
    else
      {
        MayDay::Error("AMRLevelMappedCons::define: a_coarserLevelPtr is not "
                      "castable to AMRLevelMappedCons*");
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

  if (m_coordSysFact != NULL)
    {
      RealVect dx = m_dx * RealVect::Unit;

      // Assign a multi-block coordinate system.
      m_coordSys = m_coordSysFact->getCoordSys(a_problemDomain, dx);

      m_fluxRegister = new MultiBlockFluxRegister(m_coordSysFact->getCoordSys(a_problemDomain, dx));
      m_levelGridMetricsPtr = new LevelGridMetrics(m_numStates,
                                                   m_spaceOrder);
      const bool haveVectorData =  // Only meaningful for multiblock
        (m_molPhysics->velocityInterval().size() > 0);
      m_levelGridMetricsPtr->define(this,
                                    m_coordSysFact,
                                    coarserLGM,
                                    m_dx*RealVect::Unit,
                                    m_ghostVect,
                                    LevelGridMetrics::TransverseNOpAverage,
                                    haveVectorData);
      m_levelConsOperatorPtr->defineMetricsIndices(
        m_levelGridMetricsPtr->getCoordSys(Box()));
      m_mbUtil.define(m_coordSys);
    }

  // End application-dependent code - dfm
}

//////////////////////////////////////////////////////////////////////////////

void
AMRLevelMappedCons::setDefaultValues()
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
  m_levelStep = 0;
  m_useSourceTerm = false;
  m_sourceTermPtr = NULL;
}

//////////////////////////////////////////////////////////////////////////////

// Advance by one timestep
Real AMRLevelMappedCons::advance()
{
  advanceU();
  // FIXME: Without this statement, large errors occur at domain boundaries! -JNJ
  m_data.getU();
  ++m_levelStep;
  //**FIXME We probably need to think carefully about where this is done -- here
  //**it forces lots of message passing to get <U>.
  Real returnDt = computeNewDt();
  return returnDt;
}


//////////////////////////////////////////////////////////////////////////////

// Advance by one timestep
void AMRLevelMappedCons::advanceU()
{
  if (s_verbosity >= 3)
  {
    const int savedPrec = pout().precision();
    pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
    pout().precision(6);
    pout() << "AMRLevelMappedCons::advance level " << m_level
           << " at time " << m_time << " by " << m_dt << endl;
    pout().setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
    pout().precision(savedPrec);
  }

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
      coarserFR = coarserPtr->m_fluxRegister;

      coarserUold = &coarserPtr->m_data.getJUOld();
      coarserUnew = &coarserPtr->m_data.getJU();

      tCoarserNew = coarserPtr->m_time;
      tCoarserOld = tCoarserNew - coarserPtr->m_dt;
    }

  // A finer level exists
  if (m_hasFiner)
    {
      // Recall that my flux register goes between my level and the next
      // finer level
      finerFR = m_fluxRegister;
      finerFR->setToZero();
    }

//--Fill the invalid ghost of <U> (and 1 layer of <JU>) if there is a coarser
//--level.  This is akin to calling fillGhostsRK4AndComputeU() for the first
//--stage in LevelMappedConsOperator.cpp.  <U> will be required right now if
//--artificial viscosity is added
//**FIXME Perhaps we can also do the computation of <U> for the first stage and
//**avoid it in LevelMappedConsOperator.
//**FIXME If this is the first subcycle, we should be able to just use the
//**coarsened-fine data previously set by the call to postTimeStep() for the
//**next coarser level
//**FIXME If this is not the first subcycle, we should still have the correct
//**coarsened-fine data, set at the end of the previous advance() on this level

  if (m_hasCoarser)
    {
      // Check that current fine-level time "a_time" falls between
      // the old and new coarse times.
      Real dtCoarse = tCoarserNew - tCoarserOld;
      Real alpha = (m_time - tCoarserOld)/dtCoarse;

      // Truncate the fraction to the range [0,1] to remove floating-point
      // subtraction roundoff effects
      const Real eps = 0.04 * dtCoarse/refRatio();

      if (Abs(alpha) < eps)
        {
          alpha = 0.0;
        }
      if (Abs(1.0-alpha) < eps)
        {
          alpha = 1.0;
        }
      const Interval interval(0, m_numStates-1);

      if (!m_forwardEuler)
        {
          // Explicitly mark the CrFn data used to interpolate to this level as
          // invalid (before resetting it)
          m_levelGridMetricsPtr->timeIntermediate(
            m_levelConsOperatorPtr->getTimeInterpolator(),
            alpha,
            0,  // Always stage zero since we are always outside of stages
            m_molPhysics->velocityInterval());
          // m_levelGridMetricsPtr->invalidateCr2ThisInterpolatorCrFnLevData();
          // TimeInterpolatorRK4& timeInterpolator =
          //   m_levelConsOperatorPtr->getTimeInterpolator();
          // timeInterpolator.intermediate(
          //   m_levelGridMetricsPtr->presetCr2ThisInterpolatorCrFnLevJU(),
          //   alpha,
          //   0,  // Always stage zero since we are always outside of stages
          //   interval);
        }
    }

  // Save <JU> in old location.  This means <JU> new is already initialized to
  // the current and this operation can be avoided in RK4LevelAdvance
  m_data.copyJUNewToOld();

//--Add in artificial viscosity

  if (m_useArtificialViscosity)
    {
      LevelData<FArrayBox>& UOld  = m_data.getU(m_numGhost, m_numGhost);
      LevelData<FArrayBox>& JUNew = m_data.getJU();
      // The updates the flux registers and JUNew
      m_levelConsOperatorPtr->addArtificialViscosity(JUNew,
                                                     UOld,
                                                     *finerFR,
                                                     *coarserFR,
                                                     m_time,
                                                     m_dt);
    }

//--Advance conservation law by one time step using 4th-order Runge-Kutta.

  // Real oldTime = m_time;
  if (m_forwardEuler)
    {
      //**FIXME -- not working
      LevelData<FArrayBox> RHSTmp;
      m_levelConsOperatorPtr->defineRHSData(RHSTmp, m_data.getJUOld());
      m_levelConsOperatorPtr->evalCountMax(1); // to make 1 call to evalRHS
      m_levelConsOperatorPtr->resetEvalCount();
      m_levelConsOperatorPtr->evalRHS(RHSTmp, m_data.getJUOld(),
                                      *finerFR, *coarserFR,
                                      *coarserUold, tCoarserOld,
                                      *coarserUnew, tCoarserNew,
                                      m_time, m_dt);
      m_levelConsOperatorPtr->updateODE(m_data.getJU(), RHSTmp, m_dt);
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
            LevelMappedConsOperator>(m_data.getJU(),
                                     m_data.getJUOld(),
                                     timeInterpolator,
                                     *coarserUold, tCoarserOld,
                                     *coarserUnew, tCoarserNew,
                                     *coarserFR, *finerFR,
                                     m_time,
                                     m_dt,
                                     *m_levelConsOperatorPtr,
                                     false);  // Do not initialize JU new
        }
      else
        {
          RK4LevelAdvance<LevelData<FArrayBox>,
            LevelFluxRegister,
            LevelMappedConsOperator>(m_data.getJU(),
                                     m_data.getJUOld(),
                                     *coarserUold, tCoarserOld,
                                     *coarserUnew, tCoarserNew,
                                     *coarserFR, *finerFR,
                                     m_time,
                                     m_dt,
                                     *m_levelConsOperatorPtr,
                                     false);  // Do not initialize JU new
        }
    }

   // Update the time and store the new timestep
   m_time += m_dt;

//--So that we can compute <U> properly in the invalid ghost cells on this
//--level, we have to set the CrFnU to the current time in the interpolator
//**FIXME If this is the last subcycle, then postTimeStep will also update
//**CrFnU.  Can we avoid doing this twice?  ComputeNewDt gets in the way.

   if (m_hasCoarser)
     {
       // Check that current fine-level time "a_time" falls between
       // the old and new coarse times.
       Real dtCoarse = tCoarserNew - tCoarserOld;
       Real alpha = (m_time - tCoarserOld)/dtCoarse;

       // Truncate the fraction to the range [0,1] to remove floating-point
       // subtraction roundoff effects
       const Real eps = 0.04 * dtCoarse/refRatio();

       if (Abs(alpha) < eps)
         {
           alpha = 0.0;
         }
       if (Abs(1.0-alpha) < eps)
         {
           alpha = 1.0;
         }
       const Interval interval(0, m_numStates-1);

       // Explicitly mark the CrFn data used to interpolate to this level as
       // invalid (before resetting it)
       m_levelGridMetricsPtr->timeIntermediate(
         m_levelConsOperatorPtr->getTimeInterpolator(),
         alpha,
         0,  // Always stage zero since we are always outside of stages
         m_molPhysics->velocityInterval());
       // m_levelGridMetricsPtr->invalidateCr2ThisInterpolatorCrFnLevData();
       //  TimeInterpolatorRK4& timeInterpolator =
       //    m_levelConsOperatorPtr->getTimeInterpolator();
       //  timeInterpolator.intermediate(
       //    m_levelGridMetricsPtr->presetCr2ThisInterpolatorCrFnLevJU(),
       //    alpha,
       //    0,  // Always stage zero since we are always outside of stages
       //    interval);
     }

   // All we really know is <JU> on valid cells.
   // The next finer level may perform a regrid as the next operation:
   //   1) If it is on the first sub-cycle, it needs coarse data from this level
   //      from the old time.  This was already stored in a previous
   //      postTimeStep() on this level.  The 'false' argument below indicates
   //      that this data is still valid.
   //   2) If it is on a later sub-cycle, it will have updated the coarsened-
   //      fine information it requires in the block above.
   m_data.invalidate(false);
}

//////////////////////////////////////////////////////////////////////////////

Real AMRLevelMappedCons::computeNewDt()
{
  // See "High-Order, Finite-Volume Methods in Mapped Coordinates", Colella
  // et al., JCP 2011 for derivation of this constraint.
  const Real stabilityConstraint = 1.3925;

  // Compute the maximum of the sum of the velocities in each direction.
  Real maxVelSumLoc = m_stabilityStrategy->maxWaveSpeed(m_problem_domain, m_grids, m_levelGridMetricsPtr,
                                                        *m_molPhysics, m_data, m_dx, m_time);
  Real maxVelSum = maxVelSumLoc;
#ifdef CH_MPI
  MPI_Allreduce(&maxVelSumLoc,
                &maxVelSum,
                1,
                MPI_CH_REAL,
                MPI_MAX,
                Chombo_MPI::comm);
#endif

  m_dtNew = 0.;
  if (maxVelSum > 1.E-15)
    {
      m_dtNew = m_cfl*stabilityConstraint*m_dx/maxVelSum;
    }
  if (s_verbosity >= 3)
    {
      const int savedPrec = pout().precision();
      pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
      pout().precision(std::numeric_limits<Real>::digits10);
      pout() << "AMRLevelMappedCons::computeNewDt on level " << m_level << ": "
             << m_dtNew << std::endl;
      pout().setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      pout().precision(savedPrec);
    }
  // See "High-Order, Finite-Volume Methods in Mapped Coordinates", Colella
  return m_dtNew;
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Things to do after a timestep -- reflux
/** If (!m_hasFiner) then this does nothing except write out some
 *  diagnostics depending on the verbosity level.
 *//*-----------------------------------------------------------------*/

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

  if (m_hasFiner)
    {
      // Reflux
      Real scale = -1.0/m_dx;
      if (m_levelGridMetricsPtr->isMultiblock())
        {
          m_fluxRegister->reflux(m_data.getJU(),
                                 m_molPhysics->velocityInterval(),
                                 scale);
        }
      else
        {
          m_fluxRegister->reflux(m_data.getJU(), scale);
        }

      // Average from finer level data
      AMRLevelMappedCons* amrConsFinerPtr = getFinerLevel();

      // Set m_JUnew at cell c to the mean of amrConsFinerPtr->m_JUnew
      // over all finer cells within c.
      amrConsFinerPtr->m_coarseAverage.averageToCoarse(
        m_data.getJU(),
        amrConsFinerPtr->m_data.getJU());
      // Average metric terms (these only average if required)
      if (m_levelGridMetricsPtr)
        {
          m_levelGridMetricsPtr->postTimeStep();
        }
      m_data.invalidate();
      // Set coarsened-fine data for the interpolator for the next finer level.
      // We rely on proper-nesting to keep this independent from the next
      // coarser level
      //**FIXME check if this can be removed (the code will fail an assertion if
      //**it can't be removed)
      if (m_levelGridMetricsPtr)
        {
          m_levelGridMetricsPtr->presetThis2FnInterpolatorCrFnLevU(
            m_data.getU(0, -1), m_molPhysics->velocityInterval());
        }
    }

  if (s_verbosity >= 2 && m_level == 0)
    {
      int nRefFine = 1;

      pout() << "AMRLevelMappedCons::postTimeStep:" << endl;
      pout() << "  Sums:" << endl;
      for (int comp = 0; comp < m_numStates; comp++)
        {
          Interval curComp(comp,comp);
          Real integral = computeSum(m_data.getJU(), NULL, nRefFine, m_dx,
                                     curComp);

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

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::postTimeStep " << m_level
             << " finished" << endl;
    }
//   if (m_level == 0) //DBG
//     {
//       pout() << "Post time step:\n";
//       reportNorm();
//     }
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

  //**FIXME there may be excessive copying of IVS here (depending on what the
  //**      compiler does.
  IntVectSet localTags =
    m_taggingStrategy->cellsToRefine(m_level, m_problem_domain,
                                     m_grids, m_levelGridMetricsPtr,
                                     *m_molPhysics, m_data, m_dx, m_time);
  m_mbUtil.growIVS(localTags, m_tagBufferSize);
  a_tags = localTags;
}

//////////////////////////////////////////////////////////////////////////////

// Create tags at initialization
void AMRLevelMappedCons::tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

//////////////////////////////////////////////////////////////////////////////

int AMRLevelMappedCons::indexForTagging()
{
  int ind = (m_tagPressure) ?
    m_molPhysics->pressureIndex() :
    m_molPhysics->densityIndex();
  return ind;
}

//////////////////////////////////////////////////////////////////////////////

// Reflux to make the metrics consistent between the old and new grids
void AMRLevelMappedCons::preRegrid(int                         a_base_level,
                                   const Vector<Vector<Box> >& a_new_grids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::preRegrid " << a_base_level << "  "
             << m_level << endl;
    }

  // <U> needs to be up-to-date here (on the old mesh), in the valid cells of
  // the coarse level and in the valid cells + 1 layer of valid ghost cells for
  // this level.  However, if this is the finest level on the new hierarchy,
  // then we don't need any <U>s.  Since this routine traverses fine-to-coarse,
  // we retrieve <U> on the coarser level, if it exists.
  if (m_levelGridMetricsPtr)
    {
      if (a_base_level != m_level)
        {
          m_grids = loadBalance(a_new_grids[m_level]);
        }
      // We do not require <U> and coarser <U> if this is the finest level on
      // the new hierarchy (note that new levels can only be added one at a
      // time).
      // UPtr will be dereferenced but may not be used.  Start by pointing it
      // to rawU.
      LevelData<FArrayBox>* UPtr             = &(m_data.rawU());
      const LevelData<FArrayBox>* coarseUPtr = NULL;
      if (m_level < s_finestLevel)
        {
          if (m_hasCoarser)
            {
              //**FIXME Changed until BlockRegister is used to do snapback
//            coarseUPtr = &(getCoarserLevel()->m_data.getU());
              coarseUPtr = &(getCoarserLevel()->m_data.getU(1, 1));
            }
          UPtr = &(m_data.getU(1, 0));
        }
      else
        {
          s_finestLevel = m_level;
        }
      // Use pointers/references here to clearly avoid side-effects
      LevelData<FArrayBox>& JUNew = m_data.getJU();
      m_levelGridMetricsPtr->preRegrid(a_base_level,
                                       m_grids,
                                       coarseUPtr,
                                       *UPtr,
                                       JUNew);

      // And now everything is invalid because <JU> may have changed
      m_data.invalidate();
    }
}

//////////////////////////////////////////////////////////////////////////////

// Set up data on this level after regridding
void AMRLevelMappedCons::regrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::regrid " << m_level << endl;
    }

//--If we don't receive any boxes, then this level is going away.  This has been
//--separated out to clearly illustrate what should happen

  if (a_newGrids.size() == 0)
    {
      if (s_verbosity >= 3)
        {
          pout() << "AMRLevelMappedCons::regrid: level " << m_level
                 << " is being created with 0 cells. " << endl;
        }
      // Correct tracking of finest level
      if (m_level <= s_finestLevel)
        {
          s_finestLevel = m_level - 1;
        }
      // Sanity check: if using metrics, LGM should have also detected that this
      // level is going away and the metric terms should be undefined
      if (m_levelGridMetricsPtr)
        {
          CH_assert(!m_levelGridMetricsPtr->metricsDefined());
          CH_assert(m_levelGridMetricsPtr->getBoxes().size() == 0);
        }
      m_level_grids = a_newGrids;
      // This must be done even if using metrics (if this level persists,
      // normally it would have been done in preRegrid).
      m_grids = loadBalance(a_newGrids);
      // Everything is built but with a zero-sized box layout
      // levelSetup() will correctly set m_hasFiner and m_hasCoarser for this
      // level and fix up m_hasFiner for the coarser level.
      levelSetup();
      m_data.reshapeJUNew();
      m_data.reshapeJUOld();
      m_data.reshapeU();
      m_data.invalidate();
      return;
    }

  //**FIXME This should be in library code (meshRefine itself).
  // First, intersect the given boxes with each of the blocks in our
  // coordinate system.
  RealVect dx = m_dx * RealVect::Unit;
  const Vector<Box>& blocks = m_coordSys->mappingBlocks();
  Vector<Box> validGrids;
  for (int ig = 0; ig < a_newGrids.size(); ++ig)
  {
    const Box& grid = a_newGrids[ig];
    for (int ib = 0; ib < blocks.size(); ++ib)
    {
      const Box& blockBox = blocks[ib];
      Box intersect = grid & blockBox;
      if (intersect != Box())
        validGrids.push_back(intersect);
    }
  }

  // Save original grids
  m_level_grids = validGrids;

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
  m_data.copyJUNewToOld();

  // Define the new grid if not mapped (otherwise done in preRegrid).
  if (m_levelGridMetricsPtr == 0)
    {
      m_grids = loadBalance(validGrids);
    }

  // Reshape is simple because MappedLevelData already knows about the grid
  // layout from LevelGridMetrics
  // Reshape state with new grids
  m_data.reshapeJUNew();

  // Note that U is now on a different layout from JU.  Any attempts to fill
  // <U> (beyond a simple exchange) will therefore fail.

//--Fill <JU>

  // Interpolate from coarser level
  if (AMRLevel::hasCoarserLevel())  // Use from base class since levelSetup()
    {                               //  has not yet been called
      AMRLevelMappedCons *const crAMRLevel = getCoarserLevel();
      if (m_levelGridMetricsPtr)
        {
          // Use references here to clearly avoid side-effects.  With proper
          // nesting, we only need <U> computed from valid <JU>.  I.e., 1 layer
          // away from coarse-fine interfaces is sufficient.
          LevelData<FArrayBox>& crU     = crAMRLevel->m_data.getU(0, -1);
          LevelData<FArrayBox>& crJUNew = crAMRLevel->m_data.getJU();
          // We don't required any ghosts in JU here
          m_levelGridMetricsPtr->regrid(m_data.getJU(),
                                        crU,
                                        crJUNew,
                                        m_molPhysics->velocityInterval());
        }
    }

  // Copy from old state
  m_data.copyJUOldToNew();

  // Reshape JUold to be consistent with JUnew
  m_data.reshapeJUOld();

  // And <U> needs to be reshaped too
  m_data.reshapeU();

  // And now everything is invalid because <JU> may have changed
  m_data.invalidate();

  // Set up remaining data structures (averagers, interpolators, and the level
  // operator).
  levelSetup();

  // And unit normals if using mapped grids
  if (m_levelGridMetricsPtr)
    {
      m_levelConsOperatorPtr->defineUnitNormals(m_levelGridMetricsPtr->m_N);
    }
}

//////////////////////////////////////////////////////////////////////////////

// Initialize grids
void AMRLevelMappedCons::initialGrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::initialGrid " << m_level << endl;
    }

  // First, intersect the given boxes with each of the blocks in our
  // coordinate system.
  RealVect dx = m_dx * RealVect::Unit;
  const Vector<Box>& blocks = m_coordSys->mappingBlocks();
  Vector<Box> validGrids;
  for (int ig = 0; ig < a_newGrids.size(); ++ig)
  {
    const Box& grid = a_newGrids[ig];
    for (int ib = 0; ib < blocks.size(); ++ib)
    {
      const Box& blockBox = blocks[ib];
      Box intersect = grid & blockBox;
      if (intersect != Box())
        validGrids.push_back(intersect);
    }
  }

  // Save original grids and load balance
  m_level_grids = validGrids;
  m_grids = loadBalance(validGrids);
  m_levelGridMetricsPtr->initialGrid(&m_grids);

  // Assign the coordinate system to our IBC.
  PhysMappedIBC* IBC = dynamic_cast<PhysMappedIBC*>(m_molPhysics->getPhysIBC());
  IBC->setCoordSys(m_coordSys);

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

  // Once LevelGridMetrics knows the layout, we can define the data
  const LevelData<FArrayBox>* coarseUPtr = 0;
  {
    // m_hasCoarser won't be ready until levelSetup is called
    AMRLevelMappedCons* coarserPtr = getCoarserLevel();
    if (coarserPtr)
      {
        coarseUPtr = &(coarserPtr->m_data.rawU());
      }
  }
  // Note the LGP, if used, knows the grid but no metrics are yet defined
  m_data.define(m_levelGridMetricsPtr,
                0,  // 0 means m_data will allocate required memory
                0,  // 0 means m_data will allocate required memory
                0,  // 0 means m_data will allocate required memory
                coarseUPtr,
                m_numStates,
                m_molPhysics->velocityInterval(),
                m_numGhost);

  // Set up data structures
  levelSetup();
}

//////////////////////////////////////////////////////////////////////////////

// Compute the grid metrics for mapped grids
void AMRLevelMappedCons::postInitialGrid(const bool a_restart)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::postInitialGrid " << m_level << endl;
    }

  // Post initial grid is only called for levels with grids so we can use this
  // to get the finest level
  s_finestLevel = std::max(s_finestLevel, m_level);

  if (m_levelGridMetricsPtr)
    {
      // Define the metric terms
      m_levelGridMetricsPtr->postInitialGrid(NULL);
      // Define the unit normals required for the Riemann solver
      m_levelConsOperatorPtr->defineUnitNormals(m_levelGridMetricsPtr->m_N);
    }

//--Things to do if this is a restart

  if (a_restart)
    {
      // (normally performed in initialData())

      // m_basicIBCPtr->setTime(0.0);
      if (m_levelGridMetricsPtr)
        {
          PhysMappedIBC* physIBCPtr =
            dynamic_cast<PhysMappedIBC*>(m_molPhysics->getPhysIBC());
          if (!physIBCPtr)
            {
              MayDay::Error("Coordinate system specified but the PhysIBC class "
                            "does not inherit from PhysMappedIBC");
            }
          physIBCPtr->setTime(m_time);
        }

//--Set coarsened-fine data for the interpolator to the next finer level.

      if (m_levelGridMetricsPtr && m_levelGridMetricsPtr->hasFiner())
        {
          m_levelGridMetricsPtr->presetThis2FnInterpolatorCrFnLevU(
            m_data.getU(0, -1), m_molPhysics->velocityInterval());
        }

      // (normally performed in postInitialize())

//--This is our conservation reference

      if (m_level == 0)
        {
          s_L1NormJUConsvRef.resize(m_numStates);
          computeNorm(0, s_L1NormJUConsvRef);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

// Initialize data
void AMRLevelMappedCons::initialData()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::initialData " << m_level << endl;
    }

  LevelData<FArrayBox>& U = m_data.rawU();

  // FIX:  BasicAdvectIBC requires setTime.
  // m_basicIBCPtr->setTime(0.0);
  if (m_levelGridMetricsPtr)
    {
      PhysMappedIBC* physIBCPtr =
        dynamic_cast<PhysMappedIBC*>(m_molPhysics->getPhysIBC());
      if (!physIBCPtr)
        {
          MayDay::Error("Coordinate system specified but the PhysIBC class "
                        "does not inherit from PhysMappedIBC");
        }
      physIBCPtr->setTime(m_time);
      physIBCPtr->initialize(U);
    }
  else
    {
      PhysIBC* physIBCPtr = m_molPhysics->getPhysIBC();
      physIBCPtr->initialize(U);
    }

//--Approximation to <U> if given point-wise values.  This is needed on valid +1
//--except at physical boundaries so we can calculate <JU> next.

  if (m_initialAverage)
    {
      for (DataIterator dit = m_levelGridMetricsPtr->getDataIterator();
           dit.ok(); ++dit)
        {
          Box box = m_levelGridMetricsPtr->getBoxes()[dit];
          ProblemDomain domain = m_levelGridMetricsPtr->blockDomain(box, 1);
          box.grow(1);
          box &= domain;
          FArrayBox& UFab = U[dit];
          fourthOrderAverageCell(UFab, domain, box);
        }
    }

//--Compute <JU> from <U> and <J>.  One-sided derivates are used for <U> at
//--the physical boundaries

  for (DataIterator dit = m_levelGridMetricsPtr->getDataIterator(); dit.ok();
       ++dit)
    {
      const Box interiorBox = m_levelGridMetricsPtr->getBoxes()[dit];
      FArrayBox& JUFab = m_data.rawJU()[dit];
      const FArrayBox& UFab = U[dit];
      const FArrayBox& JFab = m_levelGridMetricsPtr->m_J[dit];
      fourthOrderCellProd(JUFab, UFab, JFab, interiorBox,
                          m_levelGridMetricsPtr->blockDomain(interiorBox, 1));
    }
  m_data.invalidate();

//--Set coarsened-fine data for the interpolator from the next coarser level.
//--We have to do this from the finer level because the space-interpolator isn't
//--defined until now.  Of course, don't bother if there is actually no
//--grid on this level.

  if (m_hasCoarser && m_levelGridMetricsPtr && m_grids.size() > 0)
    {
      getCoarserLevel()->m_levelGridMetricsPtr
        ->presetThis2FnInterpolatorCrFnLevU(
        getCoarserLevel()->m_data.getU(0, -1),
        m_molPhysics->velocityInterval());
    }
}

//////////////////////////////////////////////////////////////////////////////

// Things to do after initialization
void AMRLevelMappedCons::postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::postInitialize " << m_level << endl;
    }

//--Average from finer level data (as in postTimeStep)

  if (m_hasFiner)
    {
      // Volume weighted average from finer level data
      AMRLevelMappedCons* amrConsFinerPtr = getFinerLevel();
      amrConsFinerPtr->m_coarseAverage.averageToCoarse(
        m_data.getJU(),
        amrConsFinerPtr->m_data.getJU());
      m_data.invalidate();
      // Set coarsened-fine data for the interpolator for the next finer level.
      // We rely on proper-nesting to keep this independent from the next
      // coarser level
      if (m_levelGridMetricsPtr)
        {
          m_levelGridMetricsPtr->presetThis2FnInterpolatorCrFnLevU(
            m_data.getU(0, -1), m_molPhysics->velocityInterval());
        }
    }

//--This is our conservation reference

  if (m_level == 0)
    {
      s_L1NormJUConsvRef.resize(m_numStates);
      computeNorm(0, s_L1NormJUConsvRef);
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

  // We write out all components of m_JUnew and all components of m_JUold.

  // Set up the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Set up the component names
  char compStr[32];
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr, "component_%d", comp);
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
  D_TERM6(if (m_problem_domain.isPeriodic(0))
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
      pout() << header << std::endl;
    }

  // Write the data for this level, current <JU>
  write(a_handle, m_data.rawJU().boxLayout());
  write(a_handle, m_data.rawJU(), "data");
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
      MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: checkpoint file "
                    "does not have num_components");
    }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates)
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: num_components "
                    "in checkpoint file does not match solver");
    }

  // Get the component names
  std::string stateName;
  char compStr[32];
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      if (header.m_string.find(compStr) == header.m_string.end())
        {
          MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: checkpoint "
                        "file does not have enough component names");
        }

      stateName = header.m_string[compStr];
      if (stateName != m_stateNames[comp])
        {
          MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: state_name "
                        "in checkpoint does not match solver");
        }
    }

  if (header.m_int.find("iteration") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain iteration");
    }
  m_levelStep = header.m_int ["iteration"];
}


//////////////////////////////////////////////////////////////////////////////

// Read checkpoint data for this level
void AMRLevelMappedCons::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::readCheckpointLevel" << endl;
    }
  s_finestLevel = std::max(s_finestLevel, m_level);

  // Setup the level string
  char levelStr[32];
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
    MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not "
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
    MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not "
                  "contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not "
                  "contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not "
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
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not "
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
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not "
                    "contain prob_domain");
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
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not "
                    "contain a Vector<Box>");
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

  // Reshape for new grids
  if (m_levelGridMetricsPtr)
    {
      m_levelGridMetricsPtr->initialGrid(&m_grids);
      // LevelGridMetrics now knows the grid and has defined the grids for the
      // space-interpolator
    }
  // Once LevelGridMetrics knows the layout, we can define the data
  const LevelData<FArrayBox>* coarseUPtr = 0;
  {
    // m_hasCoarser won't be ready until levelSetup is called
    AMRLevelMappedCons* coarserPtr = getCoarserLevel();
    if (coarserPtr)
      {
        coarseUPtr = &(coarserPtr->m_data.rawU());
        // A little fixup here for the level step
        //**FIXME Do you want to worry about this for newly appearing levels
        //**during a regrid?
        m_levelStep = coarserPtr->m_levelStep*coarserPtr->refRatio();
      }
  }
  m_data.define(m_levelGridMetricsPtr,
                0,  // 0 means m_data will allocate required memory
                0,  // 0 means m_data will allocate required memory
                0,  // 0 means m_data will allocate required memory
                coarseUPtr,
                m_numStates,
                m_molPhysics->velocityInterval(),
                m_numGhost);

  // Set up data structures
  levelSetup();

  // Read the data
  LevelData<FArrayBox>& JU = m_data.rawJU();
  const int dataStatus = read<FArrayBox>(a_handle,
                                         JU,
                                         "data",
                                         m_grids,
                                         JU.interval(),
                                         false);  // Do not redefine JU
  m_data.invalidate();  // Only know <JU> on valid real cells

  if (dataStatus != 0)
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointLevel: file does not "
                    "contain all state data");
    }
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
  int numStatesOutput = 2*m_numStates + 1;
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
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp + m_numStates);
      char stateNameChars[60];
      sprintf(stateNameChars, "%sJ", m_stateNames[comp].c_str());
      header.m_string[compStr] = stateNameChars;
    }
  { // J
    int comp = 2*m_numStates;
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = "J";
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
void AMRLevelMappedCons::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::writePlotLevel" << endl;
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

// Define write error plot to write errors in <JU> instead of <JU>.  Note that
// this corrupts <JU> and can only be done on conclude.
#ifdef WRITEERRORPLOT
  if (m_level == 0 && m_levelStep == WRITEERRORPLOT)
    {
      reportError();
    }
  else if (m_levelStep == 0)
    {
      m_data.getU();
    }
  const LevelData<FArrayBox>& U = m_data.rawU();
  const LevelData<FArrayBox>& JU = m_data.rawJU();
#else
  const LevelData<FArrayBox>& U = m_data.getU();
  const LevelData<FArrayBox>& JU = m_data.getJU();
#endif
  const DisjointBoxLayout& boxes = U.getBoxes();
  write(a_handle, boxes);
  if (m_tagVorticity)
    {
      MayDay::Error("Tagging vorticity not fully implemented for mapped grids");
      // // Vorticity requires 2 ghost layers.
      // LevelData<FArrayBox> outData(U.getBoxes(), m_numStates + 1,
      //                              IntVect::Unit);

      // Interval intervalU(0, m_numStates-1);
      // LevelData<FArrayBox> UnewCopy;
      // aliasLevelData(UnewCopy, &outData, intervalU);
      // U.copyTo(UnewCopy);
      // DataIterator dit = UnewCopy.dataIterator();
      // for (dit.begin(); dit.ok(); ++dit)
      //   {
      //     UnewCopy[dit].copy(m_data.getJU()[dit]);
      //   }
      // // If there is a coarser level interpolate undefined ghost cells
      // if (m_hasCoarser)
      //   {
      //     const AMRLevelMappedCons* amrLevelCoarserPtr = getCoarserLevel();
      //     const LevelData<FArrayBox>& UnewCoarser = amrLevelCoarserPtr->m_data.getJU();
      //     // need 4 layers of ghost cells for vorticity
      //     int nghostCoarseFill = 4;
      //     FourthOrderFillPatch filler(m_JUnew.disjointBoxLayout(),
      //                                 UnewCoarser.disjointBoxLayout(),
      //                                 m_numStates,
      //                                 amrLevelCoarserPtr->m_problem_domain,
      //                                 m_ref_ratio,
      //                                 nghostCoarseFill,
      //                                 true); // fixed time
      //     // No interpolation in time.
      //     filler.fillInterp(UnewCopy, UnewCoarser, 0, 0, m_numStates);
      //   }
      // UnewCopy.exchange();

      // Interval intervalVorticity(m_numStates, m_numStates);
      // LevelData<FArrayBox> vorticity;
      // aliasLevelData(vorticity, &outData, intervalVorticity);
      // computeVorticity(vorticity, UnewCopy);
      // write(a_handle, outData, "data");
    }
  else
    {
      LevelData<FArrayBox> outData(boxes, 2*m_numStates + 1);
      U.copyTo(U.interval(),
               outData,
               Interval(0, m_numStates - 1));
      JU.copyTo(JU.interval(),
                outData,
                Interval(m_numStates, 2*m_numStates - 1));
      m_levelGridMetricsPtr->m_J.copyTo(m_levelGridMetricsPtr->m_J.interval(),
                                        outData,
                                        Interval(2*m_numStates, 2*m_numStates));
      write(a_handle, outData, "data");
    }
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::writeMappedPlotFile() const
{
  // only do this on level 0
  if (m_level == 0)
    {
      int realDim = SpaceDim;
      if (s_verbosity >= 3)
        {
          pout() << "AMRLevelMappedCons::writeMappedPlotFile" << endl;
        }

      // gather AMR Levels and create node-centered dataset of
      // node locations of mapped grids

      Vector<AMRLevel*> vectAMRLevels;
      {
        // cast away const for this to call this function
        AMRLevelMappedCons* nonConstThis =
          const_cast<AMRLevelMappedCons*>(this);

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

          const MultiBlockCoordSys& levelCS = levelPtr->m_levelGridMetricsPtr->getCoordSys();
          realDim = levelCS.realDim();
          // Real levelDx = levelPtr->m_dx;
          vectNodeLoc[lev] = new LevelData<NodeFArrayBox>(levelGrids,
                                                          realDim);
          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const Box& bxBase = levelGrids[dit];

              const NewCoordSys* levelCSblock =
                levelCS.getCoordSys(bxBase);

              const RealVect& levelDx = levelCSblock->dx();

              NodeFArrayBox& thisNodeFAB = (*vectNodeLoc[lev])[dit];
              // node-centered FAB
              FArrayBox& thisFAB = thisNodeFAB.getFab();
              const Box& thisBox = thisFAB.box();
              FArrayBox xiFab(thisBox, SpaceDim);
              BoxIterator bit(thisBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  // location in index space
                  RealVect nodeIndexLoc(iv);
                  nodeIndexLoc *= levelDx;
                  // now convert to real space location
                  // RealVect nodeRealLoc = levelCSblock->realCoord(nodeIndexLoc);
                  D_EXPR6(xiFab(iv,0) = nodeIndexLoc[0],
                          xiFab(iv,1) = nodeIndexLoc[1],
                          xiFab(iv,2) = nodeIndexLoc[2],
                          xiFab(iv,3) = nodeIndexLoc[3],
                          xiFab(iv,4) = nodeIndexLoc[4],
                          xiFab(iv,5) = nodeIndexLoc[5]);
                }
              levelCSblock->realCoord(thisFAB, xiFab, thisBox);
            } // end loop over grids

          // advance to next level
          levelPtr = (AMRLevelMappedCons*) levelPtr->getFinerLevel();
        } // end loop over levels

      // create names
      Vector<string> locationNames(realDim);
      if (realDim > 0) locationNames[0] = "x";
      if (realDim > 1) locationNames[1] = "y";
      if (realDim > 2) locationNames[2] = "z";
      if (realDim > 3) locationNames[3] = "u";
      if (realDim > 4) locationNames[4] = "v";
      if (realDim > 5) locationNames[5] = "w";

      // create filename
      char suffix[64];
      sprintf(suffix, "%06d.%dd.map.hdf5", m_levelStep, SpaceDim);
      string fileName(m_plotPrefix);
      fileName += suffix;

      // Choose a dx
      Real dx = m_levelGridMetricsPtr->dxVect()[
        m_levelGridMetricsPtr->dxVect().maxDir(false)];
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

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Conclude by testing conservation if necessary
/*--------------------------------------------------------------------*/

void
AMRLevelMappedCons::conclude(int a_step) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::conclude " << m_level  << endl;
    }

  if (m_level == 0)
    {
      reportNorm();
#ifndef WRITEERRORPLOT
      if (m_levelGridMetricsPtr)
        {
          PhysMappedIBC* physIBCPtr = dynamic_cast<PhysMappedIBC*>(
            m_molPhysics->getPhysIBC());
          if (!physIBCPtr)
            {
              MayDay::Error("Coordinate system specified but the PhysIBC class "
                            "does not inherit from PhysMappedIBC");
            }
          if (physIBCPtr->haveExactSoln())
            {
              reportError();
            }
        }
#endif
    }
}

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
  Real newDT = m_initial_dt_multiplier*computeNewDt()/m_cfl;
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::computeInitialDt on level " << m_level
             << " = " << newDT << endl;
    }
  return newDT;
}

//////////////////////////////////////////////////////////////////////////////

Real AMRLevelMappedCons::getMaxWaveSpeed(const LevelData<FArrayBox>& a_U)
{
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
      Real speedOverBox = m_molPhysics->getMaxWaveSpeed(Ufab, bx);
      // Compute a running maximum
      speed = Max(speed, speedOverBox);
    }
  // Return the maximum wave speed local to this processor
  return speed;
}


//////////////////////////////////////////////////////////////////////////////

// Set the CFL number
void AMRLevelMappedCons::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Set the spatial order of accuracy
/** Can be 2 or 4 (default) FIXME: Doubtful if 2 works anymore...
 *//*-----------------------------------------------------------------*/

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

// sets whether to use high-order limiter
void AMRLevelMappedCons::highOrderLimiter(bool a_highOrderLimiter)
{
  m_highOrderLimiter = a_highOrderLimiter;
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

// Set whether refinement threshold is scaled with dx
void AMRLevelMappedCons::refinementIsScaled(bool a_refinementIsScaled)
{
  m_refinementIsScaled = a_refinementIsScaled;
}

//////////////////////////////////////////////////////////////////////////////

// Set whether to tag on pressure instead of on density
void AMRLevelMappedCons::tagPressure(bool a_tagPressure)
{
  m_tagPressure = a_tagPressure;
}

//////////////////////////////////////////////////////////////////////////////

// Set whether to tag on vorticity instead of on density
void AMRLevelMappedCons::tagVorticity(bool a_tagVorticity)
{
  m_tagVorticity = a_tagVorticity;
}

//////////////////////////////////////////////////////////////////////////////

// Set the tag buffer size
void AMRLevelMappedCons::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::useSourceTerm(bool a_useSourceTerm)
{
  m_useSourceTerm = a_useSourceTerm;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::sourceTerm(const LevelSourceTerm* const a_sourceTerm)
{
  if (a_sourceTerm != NULL)
    {
      m_sourceTermPtr = a_sourceTerm->new_sourceTerm();
    }
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::molPhysics(const MOLPhysics* const a_molPhysics)
{
  m_molPhysics = a_molPhysics->new_molPhysics();
}

//////////////////////////////////////////////////////////////////////////////

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelMappedCons::loadBalance(const Vector<Box>& a_grids)
{
  // Load balance and create boxlayout
  Vector<int> procMap;

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
    {
      pout() << "AMRLevelMappedCons::loadBalance: processor map: " << endl;
      for (int igrid = 0; igrid < a_grids.size(); ++igrid)
        {
          pout() << igrid << ": " << procMap[igrid] << "  "
                 << a_grids[igrid].volume() << endl;
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

  // In general, this can only look at coarser levels since during
  // initialization, the levels are traversed in order of coarse to fine.
  AMRLevelMappedCons* amrConsCoarserPtr = getCoarserLevel();

  // These describe whether or not a level is actually used (i.e. has a non-zero
  // number of cells) rather than just if those levels were defined.
  m_hasCoarser = AMRLevel::hasCoarserLevel();
  // During initialization, traversing coarse to fine, this will most likely be
  // set to false even if a finer level will be set up next.  So we fix up
  // what a coarser level thinks from the finer level (see below)
  m_hasFiner   = AMRLevel::hasFinerLevel();
  if (m_hasCoarser)
    {
      amrConsCoarserPtr->m_hasFiner = (m_grids.size() > 0);
    }

  // Our flux register is undefined unless a finer level is set up.  This only
  // changes a flag.
  m_fluxRegister->undefine();

  // If there is no grid on this level, don't set up anything else
  if (m_grids.size() == 0)
    {
      return;
    }

  int nRefCrse = 0;
  int numInterpolatorCrFnGhost = 0;
  const DisjointBoxLayout* coarserGridsPtr = 0;
  if (m_hasCoarser)
    {
      nRefCrse = amrConsCoarserPtr->refRatio();
      coarserGridsPtr = &amrConsCoarserPtr->m_grids;
      m_coarseAverage.define(m_grids,
                             *coarserGridsPtr,
                             m_numStates,
                             nRefCrse);

      // This may look twisted but you have to do this this way because the
      // coarser levels get set up before the finer levels so, since a flux
      // register lives between this level and the next FINER level, the finer
      // level has to do the setup because it is the only one with the
      // information at the time of construction.

      // Maintain flux registers
      amrConsCoarserPtr->m_fluxRegister->define(m_grids,
                                                *coarserGridsPtr,
                                                m_problem_domain,
                                                nRefCrse,
                                                m_numStates);
      amrConsCoarserPtr->m_fluxRegister->setToZero();

      // Find the number of ghost cells on the coarsened-fine mesh required to
      // fill the invalid ghost cells in the fine mesh.
      const IntVect& interpolatorCrFnGhostVect =
        amrConsCoarserPtr->m_levelGridMetricsPtr
        ->interpolatorCrFnNumGhost(true);
      // Assumes ghosts in all directions are equal
      numInterpolatorCrFnGhost = interpolatorCrFnGhostVect[0];
    }

  // Maintain levelMappedConsOperator
  m_levelConsOperatorPtr->define(m_levelGridMetricsPtr,
                                 &m_data.rawU(),
                                 &m_data.getUExchangeCopier(),
                                 &m_data.getJUExchangeCopier(),
                                 coarserGridsPtr,
                                 m_problem_domain,
                                 m_numGhost,
                                 numInterpolatorCrFnGhost,
                                 nRefCrse,
                                 m_dx,
                                 m_molPhysics,
                                 m_numStates,
                                 m_hasCoarser);
  if (m_sourceTermPtr != NULL)
    m_sourceTermPtr->define(m_coordSys, m_molPhysics, m_grids);

  transferSettingsToLevelOp();
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::transferSettingsToLevelOp()
{
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
  m_levelConsOperatorPtr->useSourceTerm(m_useSourceTerm);

  // If using source term, transfer m_sourceTermPtr.
  m_levelConsOperatorPtr->setSourceTerm(m_sourceTermPtr);
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::transferSettingsFromAMRLevel(AMRLevelMappedCons* a_amrConsPtr)
{
  CH_assert(a_amrConsPtr != NULL);
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
  m_minVal = a_amrConsPtr->m_minVal;
  m_domainLength = a_amrConsPtr->m_domainLength;
  m_refineThresh = a_amrConsPtr->m_refineThresh;
  m_refinementIsScaled = a_amrConsPtr->m_refinementIsScaled;
  m_tagPressure = a_amrConsPtr->m_tagPressure;
  m_tagVorticity = a_amrConsPtr->m_tagVorticity;
  m_tagBufferSize = a_amrConsPtr->m_tagBufferSize;
  m_useSourceTerm = a_amrConsPtr->m_useSourceTerm;
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

//////////////////////////////////////////////////////////////////////////////

#if 0
void AMRLevelMappedCons::computeVorticity(LevelData<FArrayBox>& a_vorticity,
                                    const LevelData<FArrayBox>& a_U) const
{
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
  //  LevelData<FArrayBox>& UnewCopy = (LevelData<FArrayBox>&) m_JUnew;
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
#endif

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
// Compute a norm
/*--------------------------------------------------------------------*/

void
AMRLevelMappedCons::computeNorm(const int a_type, Vector<Real>& a_norm) const
{
  CH_assert(m_level == 0);

  // Gather AMR Levels
  Vector<AMRLevel*> vAMRLevels;
  {
    // cast away const for this to call this function
    AMRLevelMappedCons* nonConstThis = const_cast<AMRLevelMappedCons*>(this);
    vAMRLevels = nonConstThis->getAMRLevelHierarchy();
  }
  int numLevels = vAMRLevels.size();
  // Gather ref ratios
  Vector<int> vRefRatios(numLevels - 1);
  for (int iLev = 0; iLev < numLevels - 1; ++iLev)
    {
      vRefRatios[iLev] =
        static_cast<AMRLevelMappedCons*>(vAMRLevels[iLev])->refRatio();
    }

  switch (a_type)
    {

//--Sum of all f_i for testing conservation

    case 0:
    {
      Vector<LevelData<FArrayBox>*> vPhi(numLevels);
      for (int iLev = 0; iLev != numLevels; ++iLev)
        {
          vPhi[iLev] = &(dynamic_cast<AMRLevelMappedCons*>(vAMRLevels[iLev])
                         ->m_data.getJU());
        }
      for (int iComp = 0; iComp != m_numStates; ++iComp)
        {
          const Interval compIntv(iComp, iComp);
          a_norm[iComp] = ::computeSum(vPhi,
                                        vRefRatios,
                                        1.,
                                        compIntv,
//                                        1,
                                        0);
        }
      return;
    }

    default:

//--Nothing else yet supported

    {
      MayDay::Error("computeNorm type not understood!");
      break;
    }
    }
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Compute and write out norm of conserved quantities
/** Mostly useful for debugging.  Can only be called for level 0.  The
 *  reference norm should have already been computed (i.e. only call
 *  this after post-initialize.
 *//*-----------------------------------------------------------------*/

void
AMRLevelMappedCons::reportNorm() const
{
  CH_assert(m_level == 0);
  Vector<Real> L1NormJUConsvTest(m_numStates);
  computeNorm(0, L1NormJUConsvTest);

  const int savedPrec = pout().precision();
  pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
  pout().precision(std::numeric_limits<Real>::digits10);
  pout() << "Reference and current L1 norm of conserved quantities:\n";
  for (int iComp = 0; iComp != m_numStates; ++iComp)
    {
      pout() << "  Component " << iComp << ": " << s_L1NormJUConsvRef[iComp]
             << " and " << L1NormJUConsvTest[iComp] << std::endl;
    }
  pout().setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
  pout().precision(savedPrec);
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
// Compute the error
/*--------------------------------------------------------------------*/

void
AMRLevelMappedCons::reportError() const
{
  CH_assert(m_level == 0);

  // Gather AMR Levels
  Vector<AMRLevel*> vAMRLevels;
  {
    // cast away const for this to call this function
    AMRLevelMappedCons* nonConstThis = const_cast<AMRLevelMappedCons*>(this);
    vAMRLevels = nonConstThis->getAMRLevelHierarchy();
  }
  int numLevels = vAMRLevels.size();
  // Gather ref ratios
  Vector<int> vRefRatios(numLevels - 1);
  for (int iLev = 0; iLev < numLevels - 1; ++iLev)
    {
      vRefRatios[iLev] =
        static_cast<AMRLevelMappedCons*>(vAMRLevels[iLev])->refRatio();
    }

  Vector<LevelData<FArrayBox>*> vPhi(numLevels);
  for (int iLev = 0; iLev != numLevels; ++iLev)
    {
      AMRLevelMappedCons* amrConsPtr =
        dynamic_cast<AMRLevelMappedCons*>(vAMRLevels[iLev]);
#ifdef WRITEERRORPLOT
      amrConsPtr->m_data.getU();  // Ensure <U> is known on all levels for
                                  // writing plots (ifndef WRITEERRORPLOT, this
                                  // is performed in writePlotLevel())
#endif
      vPhi[iLev] = &(amrConsPtr->m_data.getJU());
      LevelData<FArrayBox> Ux(vPhi[iLev]->disjointBoxLayout(),
                              amrConsPtr->m_numStates,
                              2*IntVect::Unit);  // U on +2, <U> on +1, <JU>
      if (amrConsPtr->m_levelGridMetricsPtr)
        {
          PhysMappedIBC* physIBCPtr = dynamic_cast<PhysMappedIBC*>(
            amrConsPtr->m_molPhysics->getPhysIBC());
          if (!physIBCPtr)
            {
              MayDay::Error("Coordinate system specified but the PhysIBC class "
                            "does not inherit from PhysMappedIBC");
            }
          physIBCPtr->setTime(m_time);
          physIBCPtr->initialize(Ux);
          for (DataIterator dit = vPhi[iLev]->dataIterator(); dit.ok(); ++dit)
            {
              const Box interiorBox = vPhi[iLev]->disjointBoxLayout()[dit];
              // Ux -> <Ux>
              Box avgBox = grow(interiorBox, 1);
              avgBox &= amrConsPtr->m_problem_domain;
              FArrayBox& UxFab = Ux[dit];
              fourthOrderAverageCell(UxFab,
                                     amrConsPtr->m_problem_domain,
                                     avgBox);
              // <JUx>
              FArrayBox JUxFab(interiorBox, amrConsPtr->m_numStates);
              const FArrayBox& JFab =
                amrConsPtr->m_levelGridMetricsPtr->m_J[dit];
              fourthOrderCellProd(JUxFab,
                                  UxFab,
                                  JFab,
                                  interiorBox,
                                  amrConsPtr->m_problem_domain);
              (*vPhi[iLev])[dit] -= JUxFab;
            }
        }
      else
        {
          MayDay::Error("Error norms only available for mapped grids");
        }
    }

  const int savedPrec = pout().precision();
  pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
  pout().precision(6);
  pout() << "Error norms at time " << m_time << std::endl;
  Real errorNorm[3];
  for (int iComp = 0; iComp != m_numStates; ++iComp)
  {
    const Interval compIntv(iComp, iComp);
    errorNorm[0] = ::computeNorm(vPhi,
                                 vRefRatios,
                                 m_dx,
                                 compIntv,
                                 0);
    errorNorm[1] = ::computeNorm(vPhi,
                                 vRefRatios,
                                 m_dx,
                                 compIntv,
                                 1);
    errorNorm[2] = ::computeNorm(vPhi,
                                 vRefRatios,
                                 m_dx,
                                 compIntv,
                                 2);
    pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
    pout().precision(std::numeric_limits<Real>::digits10);
    pout() << "  Component " << iComp << std::endl;
    pout() << "    L1 (err): " << errorNorm[1] << std::endl;
    pout() << "    L2 (err): " << errorNorm[2] << std::endl;
    pout() << "    max(err): " << errorNorm[0] << std::endl;
    pout().setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
    pout().precision(savedPrec);
  }
}

#include "NamespaceFooter.H"
