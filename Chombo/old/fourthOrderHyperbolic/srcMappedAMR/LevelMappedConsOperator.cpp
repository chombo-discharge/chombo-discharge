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

#include "LevelMappedConsOperator.H"
#include "FourthOrderUtil.H" // functions with names beginning "fourthOrder"
#include "CellToEdge.H"
#include "DivergenceF_F.H"
#include "AdvectOpF_F.H"
#include "UnitNormalsF_F.H"
#include "LevelGridMetrics.H"
#include "BlockRegister.H"
#include "GrowInBlock.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Default constructor
/** Set up some defaults.  Object requires define() to be called
 *  before all other functions.
 *//*-----------------------------------------------------------------*/

LevelMappedConsOperator::LevelMappedConsOperator()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchConsOperatorPtr = new PatchConsOperator(SpaceDim);
  m_useSourceTerm = false;
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Patch argument constructor
/** Set up some defaults.  Object requires define() to be called
 *  before all other functions.
 *//*-----------------------------------------------------------------*/

LevelMappedConsOperator::LevelMappedConsOperator(PatchConsOperator* a_operator)
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchConsOperatorPtr = a_operator;
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Destructor
/** Destroys all objects created by define(). Passed in data
 *  references of define() are left alone.
 *//*-----------------------------------------------------------------*/

LevelMappedConsOperator::~LevelMappedConsOperator()
{
  delete m_patchConsOperatorPtr;
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Define the object so that time stepping can begin (actual
//  constructor)
/** a_refine is the refinement ratio between this level and the next
 *  coarser level.  For the coarsest level, an empty DisjointBoxLayout
 *  is passed in for coarserDisjointBoxLayout.
 *  \param[in]  a_numInterpolatorCrFnGhost
 *                      The number of ghost cells from the next
 *                      coarser level required by the space
 *                      interpolator to fill the invalid ghost cells
 *                      on this level.  This is in terms of /<JU/>.
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::define(
  LevelGridMetrics *const     a_levelGridMetricsPtr,
  LevelData<FArrayBox> *const a_UPtr,
  const Copier *const         a_UExchangeCopierPtr,
  const Copier *const         a_JUExchangeCopierPtr,
  const DisjointBoxLayout*    a_coarserGrids,
  const ProblemDomain&        a_domain,
  const int                   a_numGhost,
  const int                   a_numInterpolatorCrFnGhost,
  const int                   a_refineCoarse,
  const Real                  a_dx,
  const MOLPhysics* const     a_molPhysics,
  const int                   a_numFields,
  const bool                  a_hasCoarser)
{
  // Sanity checks
  CH_assert(!a_hasCoarser || a_refineCoarse > 0);
  CH_assert(a_dx > 0.0);

  m_levelGridMetricsPtr = a_levelGridMetricsPtr;
  m_UPtr = a_UPtr;
  m_UExchangeCopierPtr = a_UExchangeCopierPtr;
  m_JUExchangeCopierPtr = a_JUExchangeCopierPtr;

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;
  m_numGhost = a_numGhost;

  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_numFields = a_numFields;  // What is the difference with numCons?

  m_doDeconvolution = true;
  m_noPPM = false;
  m_useArtificialViscosity = false;

  m_numCons      = a_molPhysics->numConserved();
  m_numFluxes    = a_molPhysics->numFluxes();
  m_velocityIntv = a_molPhysics->velocityInterval();

  // Define the time interpolator
  if (a_hasCoarser)
    {
      // Hack for multiblock -- the fine grid is a refinement of the coarse
      // grid
      if (m_levelGridMetricsPtr->isMultiblock())
        {
          // See normal case below.  Here, the "fine" grid completely covers
          // the coarse grid.  We do this so we can compute <U> everywhere
          // and then do a copyTo between blocks.  In other words, we are
          // tricking m_timeInterpolator into not really doing a copyTo so
          // that we can do it ourselves later...
          DisjointBoxLayout dummyFineGrids;
          refine(dummyFineGrids, *a_coarserGrids, m_refineCoarse);
          m_timeInterpolator.define(dummyFineGrids,
                                    *a_coarserGrids,
                                    m_domain,
                                    m_refineCoarse,
                                    m_numCons,
                                    0);  // Proper nesting means we shouldn't
                                         // need any ghosts
        }
      else
        {
          // For a single block, the copyTo performed internally in
          // m_timeInterpolator is sufficient (no need for special multiblock
          // copyTo)
          m_timeInterpolator.define(m_levelGridMetricsPtr->getBoxes(),
                                    *a_coarserGrids,
                                    m_domain,
                                    m_refineCoarse,
                                    m_numCons,
                                    a_numInterpolatorCrFnGhost);
        }
    }

  // If multiblock, define additional required members
  if (m_levelGridMetricsPtr->isMultiblock())
    {
      defineMultiblockMbrs();
    }

  GrowInBlock growTransform(m_levelGridMetricsPtr, 1);
  m_grow1inDomainLayout.deepCopy(m_levelGridMetricsPtr->getBoxes());
  m_grow1inDomainLayout.transform(growTransform);
  m_grow1inDomainLayout.closeNoSort();

  defineFlattening(); // m_flattening defined on m_grow1inDomainLayout

  //
  definePatch(a_molPhysics);

  // Everything is defined now.
  m_defined = true;
}

//////////////////////////////////////////////////////////////////////////////

// Define indices in to the metrics matrix
void
LevelMappedConsOperator::defineMetricsIndices(
  const NewFourthOrderCoordSys *const a_coordSysPtr)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      IntVect& metricTermComponentsDir = m_metricTermComponents[idir];
      for (int comp = 0; comp < SpaceDim; comp++)
        {
          metricTermComponentsDir[comp] =
            a_coordSysPtr->getNcomponent(comp, idir);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

// Find unit normals for applying the Riemann problem on mapped grids
void
LevelMappedConsOperator::defineUnitNormals(LevelData<FluxBox>& a_NLev)
{
  m_unitNormalLay.define(m_levelGridMetricsPtr->getBoxes());
  for (DataIterator dit = m_levelGridMetricsPtr->getDataIterator(); dit.ok();
       ++dit)
    {
      // We need +1 since the Riemann problem is solved on box+1
      const Box box = grow(m_levelGridMetricsPtr->getBoxes()[dit], 1);
      FluxBox& unitNormalFxb = m_unitNormalLay[dit];
      unitNormalFxb.define(box, SpaceDim*SpaceDim);
      const FluxBox &NFxb = a_NLev[dit];
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          const IntVect& metricTermComponentsDir = m_metricTermComponents[dir];
          FArrayBox& unitNormalFab = unitNormalFxb[dir];
          const FArrayBox& NFab = NFxb[dir];
          FORT_GETUNITNORMALS(CHF_FRA(unitNormalFab),
                              CHF_CONST_FRA(NFab),
                              CHF_CONST_INTVECT(metricTermComponentsDir),
                              CHF_CONST_INT(dir),
                              CHF_BOX(unitNormalFab.box()));
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Evaluate the operator \f$(-div(F))\f$ at a given time.
/** For the coarsest level, JUcoarseOld and JUcoarseNew are empty
 *  LevelData<FArrayBox> objects.  Also, either JUcoarseOld or
 *  JUcoarseNew might be empty to indicate that t(nf) = t(nc) the
 *  one grid is at the current time and interpolation is not required
 *  for boundary condition generation.  JU must be defined on the same
 *  DisjointBoxLayouts as were used in define().  Coarse flux register
 *  is flux register with the next coarser level.  Fine flux register
 *  is the flux register with the next finer level.  To the finer
 *  level FR, this level is the coarse level.  To the coarser level
 *  FR, this level is the fine level.
 *
 *  The flux registers are incremented with the normal derivatives of
 *  a_U at the boundary, times a fraction of the time step
 *  corresponding to which stage of an explicit time-stepping scheme
 *  (e.g. Runge-Kutta) is being invoked.  PATCH BY PATCH.
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::evalRHS(
  LevelData<FArrayBox>&       a_LofJU,
  LevelData<FArrayBox>&       a_JU,
  LevelFluxRegister&          a_finerFluxRegister,
  LevelFluxRegister&          a_coarserFluxRegister,
  const LevelData<FArrayBox>& a_JUcoarseOld,  // Useless
  const Real&                 a_timeCoarseOld,
  const LevelData<FArrayBox>& a_JUcoarseNew,  // Useless
  const Real&                 a_timeCoarseNew,
  Real                        a_time,
  Real                        a_weight)
{
  // Make sure everything is defined
  CH_assert(m_defined);
  m_evalCount++;

  // Added by petermc, 4 Feb 2010.
  // What an ugly kludge.  But it preserves the function interface.
  Real timeFineOrig; // beginning of fine timestep
  Real dtOrig; // full fine timestep
  switch (m_evalCount)
    {
    case 1:
      dtOrig = 6. * a_weight;
      timeFineOrig = a_time;
      break;
    case 2:
    case 3:
      dtOrig = 3. * a_weight;
      timeFineOrig = a_time - dtOrig / 2.;
      break;
    case 4:
      dtOrig = 6. * a_weight;
      timeFineOrig = a_time - dtOrig;
      break;
    default:
      MayDay::Error("LevelMappedConsOperator::evalRHS miscounting step");
    }

  // Find <U> in all cells and ghosts, store in LevelData<FArrayBox>* m_UPtr.
  fillGhostsRK4AndComputeU(a_JU, timeFineOrig, m_evalCount-1,
                           a_timeCoarseOld, a_timeCoarseNew);

  bool setFlattening = (m_useFlattening && (m_evalCount == 1));
  m_patchConsOperatorPtr->setCurrentTime(a_time);

  // Allocate a BlockRegister to resolve fluxes on block boundaries
  // so that they are single-valued.
  const DisjointBoxLayout& grids = m_levelGridMetricsPtr->getBoxes();

  BlockRegister blockRegister;
  const bool isMultiblock = m_levelGridMetricsPtr->isMultiblock();
  if (isMultiblock)
    {
      const MultiBlockCoordSys& coordSys = m_levelGridMetricsPtr->getCoordSys();
      blockRegister.define(&coordSys, grids, 0);
    }

  // Loop over the boxes ("patches") and evaluate the fluxes.
  LevelData<FluxBox> NtFAvgAll(grids,
                               m_numFluxes,
                               IntVect::Unit);
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
  {
    NtFAvgAll[dit].setVal(0.0);

    setPatchIndex(dit());

    //--Compute <F> using <U>

    const Box FBox = m_levelGridMetricsPtr->getBoxes()[dit];
    // const Box FBox1 = grow(FBox, 1);
    // const Box FBox1inDom = FBox1 & m_domain;
    // This is the full fourth-order <F>, stored as (dir, var), i.e., with
    // var contiguous
    //**Here's a question:  Why is FfaceAvg build with FBox1inDom instead of
    //**currentBox?  It seems like only current box is filled.
    const Box& FBox1inDomain = m_grow1inDomainLayout[dit];
    FluxBox FfaceAvg(FBox1inDomain, SpaceDim*m_numFluxes);
    // This is a second-order <F> computed from WfaceAvg, and used for
    // gradients
    FluxBox FfromWfaceAvg(FBox1inDomain, SpaceDim*m_numFluxes);
    m_patchConsOperatorPtr->evalFlux((*m_UPtr)[dit],
        FfaceAvg,
        FfromWfaceAvg,
        a_weight,
        setFlattening,
        m_flattening[dit]);
//cout << "U = " << endl;
//dumpFAB(&(*m_UPtr)[dit]);
#if 0
for (int d = 0; d < SpaceDim; ++d)
{
cout << "FfaceAvg[" << d << "] = " << endl;
dumpFAB(&FfaceAvg[d]);
cout << "FfromW[" << d << "] = " << endl;
dumpFAB(&FfromWfaceAvg[d]);
}
#endif
    //--Convert <F> to <NtF>

    FluxBox NtFAvg(FBox, m_numFluxes);
    const Interval NtFInterval(0, m_numFluxes-1);
    // FInterval pretends there are no components -- basically, it gives the
    // start location in FfaceAvg
    const Interval FInterval(0, m_numFluxes-1);
    const FluxBox& N = m_levelGridMetricsPtr->m_N[dit];
    const NewFourthOrderCoordSys* thisCoordSys =
      m_levelGridMetricsPtr->getCoordSys(FBox);
    thisCoordSys->computeMetricTermProductAverage(
                                                  NtFAvg,
                                                  FfaceAvg,
                                                  N,
                                                  SpaceDim,
                                                  FfromWfaceAvg,
                                                  FBox,
                                                  true,
                                                  NtFInterval,
                                                  FInterval,
                                                  0,
                                                  &m_domain);
    NtFAvgAll[dit].copy(NtFAvg, FBox);

    if (isMultiblock)
      {
        // Update the block register.
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            for (SideIterator sit; sit.ok(); ++sit)
              {
                Side::LoHiSide side = sit();
                if (blockRegister.hasInterface(dit(), idir, side))
                  {
                    blockRegister.storeFlux(NtFAvg[idir], dit(), idir, side);
                  }
              }
          }
      }
  }

  if (isMultiblock)
    {
      // Set the single-valued flux on block boundaries.
      blockRegister.close();
      setCommonFlux(NtFAvgAll, blockRegister);
    }

  // Update the flux registers.
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
  {
    FluxBox& NtFAvg = NtFAvgAll[dit];
    setPatchIndex(dit());

    // Include metrics in flux values
    /*
       Update total fluxes m_fluxes from FfaceAvg:
       m_fluxes[dit] += a_weight * NtFAvg.

       Update a_finerFluxRegister and a_coarserFluxRegister from FfaceAvg.
     */
#if 0
for (int d = 0; d < SpaceDim; ++d)
{
cout << "NtFAvg[" << d << "] = " << endl;
dumpFAB(&NtFAvg[d]);
}
#endif
    updateFluxTotalsAndRegisters(NtFAvg,
        a_finerFluxRegister, a_coarserFluxRegister,
        dit(), a_weight);

    // Find the RHS (= -div F)
    FArrayBox& LofUFab = a_LofJU[dit];
    m_patchConsOperatorPtr->evalRHS(LofUFab, NtFAvg);
  }

  if (m_useSourceTerm)
    {
      // a_LofJU += <J * S(*m_UPtr)>
      m_sourceTermPtr->addSourceTerm(a_LofJU, *m_UPtr);
    }
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  update solution -- soln += dt*rhs
/** (required by LevelRK4)
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::updateODE(LevelData<FArrayBox>& a_soln,
                                   const LevelData<FArrayBox>& a_rhs,
                                   Real a_dt)
{
  DataIterator dit = a_soln.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& solnFab = a_soln[dit];
      const FArrayBox& rhsFab = a_rhs[dit];
      m_patchConsOperatorPtr->updateODE(solnFab, rhsFab, a_dt);
    }
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  define newSoln to match existingSoln, including ghost cells
/** (required by LevelRK4)
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::defineSolnData(LevelData<FArrayBox>& a_newSoln,
                                  const LevelData<FArrayBox>& a_existingSoln)
{
  int nComp = a_existingSoln.nComp();
  const DisjointBoxLayout& grids = a_existingSoln.getBoxes();
  const IntVect& ghostVect = a_existingSoln.ghostVect();

  a_newSoln.define(grids, nComp, ghostVect);
}


//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  define RHS data based on existingSoln (in this case, w/o ghost
//  cells)
/** (required by LevelRK4)
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::defineRHSData(LevelData<FArrayBox>& a_newRHS,
                                 const LevelData<FArrayBox>& a_existingSoln)
{
  // same as defineSolnData, but w/o ghost cells
  int nComp = a_existingSoln.nComp();
  const DisjointBoxLayout& grids = a_existingSoln.getBoxes();
  const IntVect ghostVect = IntVect::Zero;

  a_newRHS.define(grids, nComp, ghostVect);
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  copy data from src->dest
/** (required by LevelRK4)
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::copySolnData(LevelData<FArrayBox>& a_dest,
                                      const LevelData<FArrayBox>& a_src)
{
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
/*--------------------------------------------------------------------*/
//  set spatial order of accuracy
/** Can be 2 or 4 (default) FIXME: Doubtful if 2 works anymore...
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
  if (m_defined)
    {
      m_patchConsOperatorPtr->spaceOrder(m_spaceOrder);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
  if (m_defined)
    {
      m_patchConsOperatorPtr->limitFaceValues(m_limitFaceValues);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::highOrderLimiter(bool a_highOrderLimiter)
{
  m_highOrderLimiter = a_highOrderLimiter;
  if (m_defined)
    {
      m_patchConsOperatorPtr->highOrderLimiter(m_highOrderLimiter);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
  if (m_defined)
    {
      m_patchConsOperatorPtr->useFlattening(m_useFlattening);
    }
  defineFlattening();
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
  if (m_defined)
    {
      m_patchConsOperatorPtr->noPPM(m_noPPM);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
  if (m_defined)
    {
      m_patchConsOperatorPtr->doDeconvolution(m_doDeconvolution);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
  if (m_defined)
    {
      m_patchConsOperatorPtr->doFaceDeconvolution(m_doFaceDeconvolution);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
  if (m_defined)
    {
      m_patchConsOperatorPtr->useArtificialViscosity(m_useArtificialViscosity);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
  if (m_defined)
    {
      m_patchConsOperatorPtr->artificialViscosity(m_artificialViscosity);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::useSourceTerm(bool a_useSourceTerm)
{
  m_useSourceTerm = a_useSourceTerm;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::defineFlattening()
{
  // Need to define m_flattening because we'll take m_flattening[dit].
  // if (m_useFlattening)
  // The boxes of flatteningLayout are
  // the boxes of the layout grown by 1 and intersected with m_domain.
  m_flattening.define(m_grow1inDomainLayout, 1);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::definePatch(const MOLPhysics* const   a_molPhysics)
{
  // PatchConsOperator m_patchConsOperator;
  m_patchConsOperatorPtr->define(m_domain, m_dx, m_levelGridMetricsPtr, a_molPhysics, m_numFields);
  m_patchConsOperatorPtr->spaceOrder(m_spaceOrder);
  m_patchConsOperatorPtr->limitFaceValues(m_limitFaceValues);
  m_patchConsOperatorPtr->highOrderLimiter(m_highOrderLimiter);
  m_patchConsOperatorPtr->useFlattening(m_useFlattening);
  m_patchConsOperatorPtr->noPPM(m_noPPM);
  m_patchConsOperatorPtr->doDeconvolution(m_doDeconvolution);
  m_patchConsOperatorPtr->doFaceDeconvolution(m_doFaceDeconvolution);
  m_patchConsOperatorPtr->useArtificialViscosity(m_useArtificialViscosity);
  m_patchConsOperatorPtr->numGhost(m_numGhost);
  // Need give access to the unit normals
  if (m_levelGridMetricsPtr)
    {
      m_patchConsOperatorPtr->unitNormals(&m_unitNormalLay);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::resetEvalCount()
{
  m_evalCount = 0;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::evalCountMax(int a_evalCountMax)
{
  m_evalCountMax = a_evalCountMax;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::addArtificialViscosity(
  LevelData<FArrayBox>&       a_JUnew,
  const LevelData<FArrayBox>& a_Uold,
  LevelFluxRegister&          a_finerFluxRegister,
  LevelFluxRegister&          a_coarserFluxRegister,
  Real                        a_oldTime,
  Real                        a_weight)
{
  CH_assert(m_useArtificialViscosity);
  m_patchConsOperatorPtr->setCurrentTime(a_oldTime);

  DataIterator dit = a_Uold.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());  // For mapped grids, we pass in the box but this
                             // is still needed to get the unit normals
      const Box& box = m_levelGridMetricsPtr->getBoxes()[dit];

      FArrayBox& JUnewFab      = a_JUnew[dit];
      const FArrayBox& UoldFab = a_Uold[dit];
      FluxBox& Nflb            = m_levelGridMetricsPtr->m_N[dit];
      FArrayBox& Jfab          = m_levelGridMetricsPtr->m_J[dit];

      FluxBox artViscNtF(box, m_numFluxes);

      // Compute the flux and update JUnew
      m_patchConsOperatorPtr->addMappedArtificialViscosity(JUnewFab,
                                                           UoldFab,
                                                           artViscNtF,
                                                           Nflb,
                                                           Jfab,
                                                           box,
                                                           a_weight);

      // Update the flux registers
      updateFluxTotalsAndRegisters(artViscNtF,
                                   a_finerFluxRegister, a_coarserFluxRegister,
                                   dit(),
                                   a_weight);
    }
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  update a_finerFluxRegister and a_coarserFluxRegister
/** \param[in]  a_Fface used to update flux registers:  it is shifted
 *                      and shifted back, but data remain unchanged
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::updateFluxTotalsAndRegisters(
  FluxBox&           a_Fface,
  LevelFluxRegister& a_finerFluxRegister,
  LevelFluxRegister& a_coarserFluxRegister,
  const DataIndex&   a_dataIndex,
  Real               a_weight)
{
  Interval intvlF(0, m_numFluxes-1);
  // Do flux register updates
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // This gets shifted and shifted back, with no data change.
      FArrayBox& currentFluxDir = a_Fface[idir];

      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level.  The finer flux register must not
      // be defined if there is no finer level.
      if (a_finerFluxRegister.isDefined())
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
LevelMappedConsOperator::getTimeInterpolator()
{
  CH_assert(m_defined);
  return m_timeInterpolator;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::fillGhostsRK4AndComputeU(
  LevelData<FArrayBox>& a_JU,
  const Real&           a_time,
  int                   a_stage,
  const Real&           a_timeCoarseOld,
  const Real&           a_timeCoarseNew)
{

//--Fill invalid ghosts of <U> (and 1 layer of <JU>) if there is a coarser
//--level

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
          MayDay::Error( "LevelMappedConsOperator::evalRHS: alpha < 0.0");
        }
      // Current time after new coarse time
      if (alpha > 1.0)
        {
          MayDay::Error( "LevelMappedConsOperator::evalRHS: alpha > 1.0");
        }
      const Interval interval(0, m_numCons-1);

      // The time interpolator directly modified a coarsened-fine <JU>.  Use the
      // time-interpolator to directly set the coarsened-fine data in the
      // space-interpolator
      //**FIXME If this is stage 1, this was probably already done.
      m_levelGridMetricsPtr->timeIntermediate(m_timeInterpolator,
                                              alpha,
                                              a_stage,
                                              m_velocityIntv);
      // m_levelGridMetricsPtr->invalidateCr2ThisInterpolatorCrFnLevData();
      // m_timeInterpolator.intermediate(
      //   m_levelGridMetricsPtr->presetCr2ThisInterpolatorCrFnLevJU(),
      //   alpha,
      //   a_stage,
      //   interval);

      // Amazingly, we can now use one routine which will find <U> and fill
      // the ghosts.  Note that this also fills 1 layer of ghosts in <JU>.
      m_levelGridMetricsPtr->fillFineGhostCells(*m_UPtr, a_JU);
    }

//**DON'T CHANGE THIS!  We just filled the invalid ghosts and now have to
//**compute <U> in valid cells, valid ghosts, and extra-block ghosts

//--Compute <U> in all valid cells.  For this, we need to exchange 1 layer of
//--ghosts in <JU>

  a_JU.exchange(*m_JUExchangeCopierPtr);
  // Note: block boundaries are treated the same as domain boundaries
  m_levelGridMetricsPtr->computeValidU(*m_UPtr, a_JU);

//--Fill all the valid ghosts of <U> by exchange

  m_UPtr->exchange(*m_UExchangeCopierPtr);

//--Fill the extra-block ghosts of <U>

  m_levelGridMetricsPtr->multiblockExchangeU(*m_UPtr, m_velocityIntv);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::setPatchIndex(const DataIndex&  a_ind) const
{
  const Box& curBox = m_levelGridMetricsPtr->getBoxes()[a_ind];
  m_patchConsOperatorPtr->setCurrentBox(a_ind, curBox);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::setSourceTerm(LevelSourceTerm* a_sourceTermPtr)
{
  m_sourceTermPtr = a_sourceTermPtr;
}

//////////////////////////////////////////////////////////////////////////////

#include "NamespaceFooter.H"
