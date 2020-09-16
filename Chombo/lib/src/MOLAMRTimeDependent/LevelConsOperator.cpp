#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "LevelConsOperator.H"
#include "FourthOrderUtil.H" // functions with names beginning "fourthOrder"
// #include "AMRIO.H"
#include "CellToEdge.H"
#include "DivergenceF_F.H"
#include "AdvectOpF_F.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
LevelConsOperator::LevelConsOperator()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchConsOperatorPtr = new PatchConsOperator();
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
LevelConsOperator::~LevelConsOperator()
{
  delete m_patchConsOperatorPtr;
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void LevelConsOperator::define(const DisjointBoxLayout&  a_thisLayout,
                               const DisjointBoxLayout&  a_coarserLayout,
                               const ProblemDomain&      a_domain,
                               const int&                a_refineCoarse,
                               const Real&               a_dx,
                               const MOLPhysics* const a_molPhysics,
                               const int&                a_numFields,
                               const bool&               a_hasCoarser,
                               const bool&               a_hasFiner)
{
  CH_TIME("LevelConsOperator::define");
  // Sanity checks
  CH_assert(a_refineCoarse > 0);
  CH_assert(a_dx > 0.0);

  // Make a copy of the current grids
  m_grids  = a_thisLayout;

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;

  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;
  m_numFields = a_numFields;

  // Set up the interpolator if there is a coarser level
  // WAS QuadCFInterp m_patcher;
  if (m_hasCoarser)
    {
      // m_patcher for interpolation from next coarser level to this level
      //      m_patcher.define(m_grids,
      //                       &a_coarserLayout,
      //                       m_dx,
      //                       a_refineCoarse,
      //                       m_numFields,
      //                       a_domain);
      ProblemDomain domainCoarse = coarsen(m_domain, a_refineCoarse);
      // FourthOrderFillPatch m_patcher;
      m_patcher.define(m_grids,
                       a_coarserLayout,
                       m_numFields,
                       domainCoarse,
                       m_refineCoarse,
                       m_numGhost);
    }

  defineFlattening();

  m_doDeconvolution = true;
  m_noPPM = false;
  m_useArtificialViscosity = false;

  m_numCons   = a_molPhysics->numConserved();
  m_numFluxes = a_molPhysics->numFluxes();

  // allocate storage for fluxes
  // LevelData<FluxBox> m_fluxes;
  m_fluxes.define(m_grids, m_numFluxes);

  definePatch(a_molPhysics);

  defineExchanger();

  // Everything is defined now.
  m_defined = true;
}


//////////////////////////////////////////////////////////////////////////////
// Evaluate the operator (-div(F) ) at time a_time,
// on conserved variables a_U:  at coarser level a_UcoarseOld, a_UcoarseNew.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// The flux registers are incremented with the normal derivatives of a_U
// at the boundary, times a fraction of the time step corresponding to
// which stage of an explicit time-stepping scheme (e.g. Runge-Kutta)
// is being invoked.  PATCH BY PATCH.
void
LevelConsOperator::evalRHS(
                           LevelData<FArrayBox>&       a_LofU,
                           LevelData<FArrayBox>&       a_U,
                           LevelFluxRegister&          a_finerFluxRegister,
                           LevelFluxRegister&          a_coarserFluxRegister,
                           const LevelData<FArrayBox>& a_UcoarseOld,
                           const Real&                 a_timeCoarseOld,
                           const LevelData<FArrayBox>& a_UcoarseNew,
                           const Real&                 a_timeCoarseNew,
                           Real                        a_time,
                           Real                        a_weight)
{
  CH_TIME("LevelConsOperator::evalRHS");
  // Make sure everything is defined
  CH_assert(m_defined);
  m_evalCount++;

  // Added by petermc, 4 Feb 2010.
  // What an ugly kludge.  But it preserves the function interface.
  // RK4LevelAdvance calls this function, and we want to get
  // timeFineOrig and dtOrig out of a_time and a_weight.
  Real timeFineOrig; // beginning of fine timestep
  Real dtOrig; // full fine timestep
  switch (m_evalCount)
    {
    case 1:
      // a_time == timeFineOrig, a_weight == dtOrig/6
      dtOrig = 6 * a_weight;
      timeFineOrig = a_time;
      break;
    case 2:
    case 3:
      // a_time == timeFineOrig + dtOrig/2, a_weight == dtOrig/3
      dtOrig = 3 * a_weight;
      timeFineOrig = a_time - dtOrig / 2.;
      break;
    case 4:
      // a_time == timeFineOrig + dtOrig, a_weight == dtOrig/6
      dtOrig = 6 * a_weight;
      timeFineOrig = a_time - dtOrig;
      break;
    default:
      MayDay::Error("LevelConsOperator::evalRHS miscounting step");
    }

  if (m_forwardEuler)
    {
      fillGhosts(a_U, a_time, a_timeCoarseOld, a_timeCoarseNew);
    }
  else
    {
      fillGhostsRK4(a_U, timeFineOrig, m_evalCount-1,
                    a_timeCoarseOld, a_timeCoarseNew);
    }

  m_patchConsOperatorPtr->setCurrentTime(a_time);

  evalRHSpatches(a_LofU, a_U,
                 a_finerFluxRegister, a_coarserFluxRegister, a_weight);
}


//////////////////////////////////////////////////////////////////////////////
void
LevelConsOperator::evalRHSpatches(
                                  LevelData<FArrayBox>&       a_LofU,
                                  const LevelData<FArrayBox>& a_U,
                                  LevelFluxRegister&          a_finerFluxRegister,
                                  LevelFluxRegister&          a_coarserFluxRegister,
                                  Real                        a_weight)
{
  CH_TIME("LevelConsOperator::evalRHSpatches");
  bool setFlattening = (m_useFlattening && (m_evalCount == 1));
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());

      FluxBox FfaceAvg;
      m_patchConsOperatorPtr->getNormalFlux(FfaceAvg,
                                            a_U[dit],
                                            a_weight,
                                            setFlattening,
                                            m_flattening[dit]);

      /*
        Update total fluxes m_fluxes from FfaceAvg:
        m_fluxes[dit] += a_weight * FfaceAvg.

        Also, if there is a finer level,
        update a_finerFluxRegister using FfaceAvg.

        Also, if there is a coarser level,
        update a_coarserFluxRegister using FfaceAvg.
      */
      updateFluxTotalsAndRegisters(FfaceAvg,
                                   a_finerFluxRegister, a_coarserFluxRegister,
                                   dit(), a_weight);
      /*
        Set LofUFab = - div(FfaceAvg).
       */
      FArrayBox& LofUFab = a_LofU[dit];
      // Set LofUFab = div(FfaceAvg).
      m_patchConsOperatorPtr->getFluxDivergence(LofUFab, FfaceAvg);
      // Actually want -div.
      LofUFab.negate();
    }
  // added 22 Oct 2008:  these change nothing
  // a_LofU.exchange();
  // a_U.exchange();
}


//////////////////////////////////////////////////////////////////////////////
void
LevelConsOperator::updateODE(LevelData<FArrayBox>& a_soln,
                             const LevelData<FArrayBox>& a_rhs,
                             Real a_dt)
{
  CH_TIME("LevelConsOperator::updateODE");
  DataIterator dit = a_soln.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& solnFab = a_soln[dit];
      const FArrayBox& rhsFab = a_rhs[dit];
      m_patchConsOperatorPtr->updateODE(solnFab, rhsFab, a_dt);
    }
}


//////////////////////////////////////////////////////////////////////////////
/// reset fluxes contained in this object to zero
void
LevelConsOperator::resetFluxes()
{
  CH_TIME("LevelConsOperator::resetFluxes");
  DataIterator dit = m_fluxes.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_fluxes[dit].setVal(0.0);
    }
}


//////////////////////////////////////////////////////////////////////////////
/// returns reference to fluxes computed and accumulated by this operator.
LevelData<FluxBox>&
LevelConsOperator::getFluxes()
{
  return m_fluxes;
}


//////////////////////////////////////////////////////////////////////////////
void
LevelConsOperator::defineSolnData(LevelData<FArrayBox>& a_newSoln,
                                  const LevelData<FArrayBox>& a_existingSoln)
{
  CH_TIME("LevelConsOperator::defineSolnData");
  int nComp = a_existingSoln.nComp();
  const DisjointBoxLayout& grids = a_existingSoln.getBoxes();
  const IntVect& ghostVect = a_existingSoln.ghostVect();

  a_newSoln.define(grids, nComp, ghostVect);
}


//////////////////////////////////////////////////////////////////////////////
void
LevelConsOperator::defineRHSData(LevelData<FArrayBox>& a_newRHS,
                                 const LevelData<FArrayBox>& a_existingSoln)
{
  CH_TIME("LevelConsOperator::defineRHSData");
  // same as defineSolnData, but w/o ghost cells
  int nComp = a_existingSoln.nComp();
  const DisjointBoxLayout& grids = a_existingSoln.getBoxes();
  const IntVect ghostVect = IntVect::Zero;

  a_newRHS.define(grids, nComp, ghostVect);
}


//////////////////////////////////////////////////////////////////////////////
/// copy data src->dest
void
LevelConsOperator::copySolnData(LevelData<FArrayBox>& a_dest,
                                const LevelData<FArrayBox>& a_src)
{
  CH_TIME("LevelConsOperator::copySolnData");
  if (a_dest.disjointBoxLayout() == a_src.disjointBoxLayout())
    { // Do this to copy all the ghost cells too.
      for (DataIterator dit = a_src.dataIterator(); dit.ok(); ++dit)
        {
          const FArrayBox& srcFab = a_src[dit];
          FArrayBox& destFab = a_dest[dit];
          destFab.copy(srcFab);
        }
    }
  else
    {
      a_src.copyTo(a_dest);
    }
}


//////////////////////////////////////////////////////////////////////////////
void
LevelConsOperator::exchangeGhosts(LevelData<FArrayBox>&   a_U)
{
  CH_TIME("LevelConsOperator::exchangeGhosts");
  a_U.exchange(m_exchangeCopier);
}

//////////////////////////////////////////////////////////////////////////////
//  set spatial order of accuracy
/** Can be 2 or 4 (default) FIXME: Doubtful if 2 works anymore...*/
void
LevelConsOperator::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
  if (m_defined) m_patchConsOperatorPtr->spaceOrder(m_spaceOrder);
}

//////////////////////////////////////////////////////////////////////////////
void
LevelConsOperator::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
  if (m_defined) m_patchConsOperatorPtr->limitFaceValues(m_limitFaceValues);
}

//////////////////////////////////////////////////////////////////////////////
void
LevelConsOperator::highOrderLimiter(bool a_highOrderLimiter)
{
  m_highOrderLimiter = a_highOrderLimiter;
  if (m_defined) m_patchConsOperatorPtr->highOrderLimiter(m_highOrderLimiter);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
  if (m_defined) m_patchConsOperatorPtr->useFlattening(m_useFlattening);
  defineFlattening();
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
  if (m_defined) m_patchConsOperatorPtr->noPPM(m_noPPM);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
  if (m_defined) m_patchConsOperatorPtr->doDeconvolution(m_doDeconvolution);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
  if (m_defined) m_patchConsOperatorPtr->doFaceDeconvolution(m_doFaceDeconvolution);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
  if (m_defined) m_patchConsOperatorPtr->useArtificialViscosity(m_useArtificialViscosity);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
  if (m_defined) m_patchConsOperatorPtr->artificialViscosity(m_artificialViscosity);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::forwardEuler(bool a_forwardEuler)
{
  m_forwardEuler = a_forwardEuler;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::numGhost(int a_numGhost)
{
  m_numGhost = a_numGhost;
  if (m_defined) m_patchConsOperatorPtr->numGhost(m_numGhost);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::defineFlattening()
{
  CH_TIME("LevelConsOperator::defineFlattening");
  // Need to define m_flattening because we'll take m_flattening[dit].
  //  if (m_useFlattening)
  // The boxes of flatteningLayout are
  // the boxes of m_grids grown by 1 and intersected with m_domain.
  BoxLayout flatteningLayout;
  flatteningLayout.deepCopy(m_grids);
  flatteningLayout.grow(1);
  flatteningLayout &= m_domain;
  flatteningLayout.close();
  m_flattening.define(flatteningLayout, 1);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::definePatch(const MOLPhysics* const   a_molPhysics)
{
  CH_TIME("LevelConsOperator::definePatch");
  // PatchConsOperator m_patchConsOperator;
  m_patchConsOperatorPtr->define(m_domain, m_dx, a_molPhysics, m_numFields);
  m_patchConsOperatorPtr->spaceOrder(m_spaceOrder);
  m_patchConsOperatorPtr->limitFaceValues(m_limitFaceValues);
  m_patchConsOperatorPtr->highOrderLimiter(m_highOrderLimiter);
  m_patchConsOperatorPtr->useFlattening(m_useFlattening);
  m_patchConsOperatorPtr->noPPM(m_noPPM);
  m_patchConsOperatorPtr->doDeconvolution(m_doDeconvolution);
  m_patchConsOperatorPtr->doFaceDeconvolution(m_doFaceDeconvolution);
  m_patchConsOperatorPtr->useArtificialViscosity(m_useArtificialViscosity);
  m_patchConsOperatorPtr->numGhost(m_numGhost);
}


//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::defineExchanger()
{
  CH_TIME("LevelConsOperator::defineExchanger");
  m_exchangeCopier.exchangeDefine(m_grids, m_numGhost*IntVect::Unit);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::resetEvalCount()
{
  m_evalCount = 0;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::evalCountMax(int a_evalCountMax)
{
  m_evalCountMax = a_evalCountMax;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::addArtificialViscosity(LevelData<FArrayBox>&         a_Unew,
                                          const LevelData<FArrayBox>&   a_Uold,
                                          LevelFluxRegister&            a_finerFluxRegister,
                                          LevelFluxRegister&            a_coarserFluxRegister,
                                          Real  a_oldTime,
                                          Real  a_weight)
{
  CH_TIME("LevelConsOperator::addArtificialViscosity");
  CH_assert(m_useArtificialViscosity);
  m_patchConsOperatorPtr->setCurrentTime(a_oldTime);

  DataIterator dit = a_Uold.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& UnewFab = a_Unew[dit];
      const FArrayBox& UoldFab = a_Uold[dit];

      setPatchIndex(dit());
      FluxBox currentFlux;
      m_patchConsOperatorPtr->addArtificialViscosity(UnewFab,
                                                     UoldFab,
                                                     currentFlux,
                                                     a_weight);

      /*
        Update total fluxes m_fluxes from FfaceAvg:
        m_fluxes[dit] += a_weight * currentFlux.

        Also, if there is a finer level,
        update a_finerFluxRegister using currentFlux.

        Also, if there is a coarser level,
        update a_coarserFluxRegister using currentFlux.
      */
      updateFluxTotalsAndRegisters(currentFlux,
                                   a_finerFluxRegister, a_coarserFluxRegister,
                                   dit(), a_weight);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::updateFluxTotalsAndRegisters(// used to update flux registers:  it is shifted and shifted back, but data remain unchanged
                                                FluxBox&  a_Fface,
                                                LevelFluxRegister&   a_finerFluxRegister,
                                                LevelFluxRegister&   a_coarserFluxRegister,
                                                const DataIndex&     a_dataIndex,
                                                Real                 a_weight)
{
  CH_TIME("LevelConsOperator::updateFluxTotalsAndRegisters");
  // thisTotalFlux[idir] += a_weight * a_Fface[idir]
  FluxBox& thisTotalFlux = m_fluxes[a_dataIndex];

  Interval intvlF(0, m_numFluxes-1);
  // Do flux register updates
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // This gets shifted and shifted back, with no data change.
      FArrayBox& currentFluxDir = a_Fface[idir];

      // First update m_fluxes:  thisTotalFluxDir += a_weight * currentFluxDir
      FArrayBox& thisTotalFluxDir = thisTotalFlux[idir];
      thisTotalFluxDir.plus(currentFluxDir, a_weight);

      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level
      if (m_hasFiner && a_finerFluxRegister.isDefined())
        {
          a_finerFluxRegister.incrementCoarse(currentFluxDir, a_weight,
                                              a_dataIndex,
                                              intvlF, intvlF, idir);
        }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
      if (m_hasCoarser && a_coarserFluxRegister.isDefined())
        {
          a_coarserFluxRegister.incrementFine(currentFluxDir, a_weight,
                                              a_dataIndex,
                                              intvlF, intvlF, idir);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

TimeInterpolatorRK4&
LevelConsOperator::getTimeInterpolator()
{
  CH_assert(m_defined);
  return m_patcher.getTimeInterpolator();
}

//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::fillGhosts(LevelData<FArrayBox>&   a_U,
                              const Real&             a_time,
                              const Real&             a_timeCoarseOld,
                              const Real&             a_timeCoarseNew)
{
  CH_TIME("LevelConsOperator::fillGhosts");
  if (m_hasCoarser)
    {
      // Check that current fine-level time "a_time" falls between
      // the old and new coarse times.
      Real dtCoarse = a_timeCoarseNew - a_timeCoarseOld;
      Real alpha = (a_time - a_timeCoarseOld) / dtCoarse;

      // Truncate the fraction to the range [0,1] to remove floating-point
      // subtraction roundoff effects
      Real eps = 0.04 * dtCoarse / m_refineCoarse;

      if (Abs(alpha) < eps)
        {
          alpha = 0.0;
        }
      if (Abs(1.0-alpha) < eps)
        {
          alpha = 1.0;
        }

      // Current time before old coarse time
      if (alpha < 0.0)
        {
          MayDay::Error( "LevelConsOperator::evalRHS: alpha < 0.0");
        }
      // Current time after new coarse time
      if (alpha > 1.0)
        {
          MayDay::Error( "LevelConsOperator::evalRHS: alpha > 1.0");
        }

      // OK, now we have alpha in the interval [0:1].

      // Interpolate ghost cells from next coarser level using both space
      // and time interpolation.

      // We need 4 ghost cell layers of fine grids.
      // To get those, we need 3 ghost cell layers of coarsened fine grids.

      // first, interpolate coarse data to the current time
      //        LevelData<FArrayBox> Ucoarse;
      //        Ucoarse.define(a_UcoarseOld);
      // timeInterp() is a function in BoxTools
      //        timeInterp(Ucoarse       ,a_time,
      //                   a_UcoarseOld  ,a_timeCoarseOld,
      //                   a_UcoarseNew  ,a_timeCoarseNew,
      //                   intvlU);
      // use current-time coarse data to fill fine ghost cells
      // m_patcher.coarseFineInterp(a_U, Ucoarse);
      m_patcher.fillInterp(a_U, alpha, 0, 0, m_numCons);
    }

  // do domain BC's
  // FIX:  In periodic case, this is a no-op.
  // PhysIBC* physIBCPtr = m_molPhysics->getPhysIBC();
  // physIBCPtr->ghostCellBC(a_U, a_time);

  // Exchange all the ghost cell data between grids at this level
  exchangeGhosts(a_U);
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}


//////////////////////////////////////////////////////////////////////////////

void
LevelConsOperator::fillGhostsRK4(LevelData<FArrayBox>&   a_U,
                                 const Real&             a_time,
                                 int                     a_stage,
                                 const Real&             a_timeCoarseOld,
                                 const Real&             a_timeCoarseNew)
{
  CH_TIME("LevelConsOperator::fillGhostsRK4");
  if (m_hasCoarser)
    {
      // Check that current fine-level time "a_time" falls between
      // the old and new coarse times.
      Real dtCoarse = a_timeCoarseNew - a_timeCoarseOld;
      Real alpha = (a_time - a_timeCoarseOld) / dtCoarse;

      // Truncate the fraction to the range [0,1] to remove floating-point
      // subtraction roundoff effects
      Real eps = 0.04 * dtCoarse / m_refineCoarse;

      if (Abs(alpha) < eps)
        {
          alpha = 0.0;
        }
      if (Abs(1.0-alpha) < eps)
        {
          alpha = 1.0;
        }
      // Current time before old coarse time
      if (alpha < 0.0)
        {
          MayDay::Error( "LevelConsOperator::evalRHS: alpha < 0.0");
        }
      // Current time after new coarse time
      if (alpha > 1.0)
        {
          MayDay::Error( "LevelConsOperator::evalRHS: alpha > 1.0");
        }
      m_patcher.fillRK4Intermediate(a_U, alpha, a_stage, 0, 0, m_numCons);
    }
  // Exchange all the ghost cell data between grids at this level
  exchangeGhosts(a_U);
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////
void LevelConsOperator::setPatchIndex(const DataIndex&  a_ind) const
{
  CH_TIME("LevelConsOperator::setPatchIndex");
  const Box& curBox = m_grids[a_ind];
  m_patchConsOperatorPtr->setCurrentBox(curBox);
}

#include "NamespaceFooter.H"
