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
  m_useSourceTerm = false;
  m_sourceTermPtr = NULL;
  m_useArtificialDissipation = false;
  m_artificialDissipation = 0.;
  m_patchMappedConsOperatorPtr = NULL;
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Patch argument constructor
/** Set up some defaults.  Object requires define() to be called
 *  before all other functions.
 *//*-----------------------------------------------------------------*/

LevelMappedConsOperator::LevelMappedConsOperator(PatchMappedConsOperator* a_operator)
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchConsOperatorPtr = (PatchConsOperator*) m_patchMappedConsOperatorPtr;
}

//////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------*/
//  Destructor
/** Destroys all objects created by define(). Passed in data
 *  references of define() are left alone.
 *//*-----------------------------------------------------------------*/

LevelMappedConsOperator::~LevelMappedConsOperator()
{
  delete m_patchMappedConsOperatorPtr;
  // Since we copied this into the base class cons op ptr, 
  // just set to NULL so the base class destructor completes
  m_patchConsOperatorPtr = NULL;
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
  CH_TIME("LevelMappedConsOperator::define");
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
  m_vectorFluxIntv = a_molPhysics->vectorFluxInterval();

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
  CH_TIME("LevelMappedConsOperator::defineMetricsIndices");
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
  CH_TIME("LevelMappedConsOperator::defineUnitNormals");
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
  CH_TIME("LevelMappedConsOperator::evalRHS");
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
  m_patchMappedConsOperatorPtr->setCurrentTime(a_time);

  // Allocate a BlockRegister to resolve fluxes on block boundaries
  // so that they are single-valued.
  const DisjointBoxLayout& grids = m_levelGridMetricsPtr->getBoxes();

  BlockRegister blockRegister;
  const bool isMultiblock = m_levelGridMetricsPtr->isMultiblock();
  if (isMultiblock)
    {
      const MultiBlockCoordSys& coordSys = m_levelGridMetricsPtr->getCoordSys();
      // Need a ghost cell when and only when there are vector fluxes.
      int blockRegisterGhost = (m_vectorFluxIntv.size() > 0) ? 1 : 0;
      blockRegister.define(&coordSys, grids, blockRegisterGhost);
    }

  if (m_useArtificialDissipation)
    {
      Real maxSpeed = getMaxSpeed();
      m_patchMappedConsOperatorPtr->setMaxSpeed(maxSpeed);
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
    m_patchMappedConsOperatorPtr->evalFlux((*m_UPtr)[dit],
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
    ProblemDomain probDomain(m_domain);
    if (isMultiblock)
      {
        probDomain = ProblemDomain(FBox1inDomain);
      }
    thisCoordSys->
      computeMetricTermProductAverage(
                                      NtFAvg,
                                      FfaceAvg,
                                      N,
                                      SpaceDim,
                                      FfromWfaceAvg,
                                      FBox,
                                      true, // 4th-order
                                      NtFInterval,
                                      FInterval,
                                      0,
                                      &probDomain);
    NtFAvgAll[dit].copy(NtFAvg, FBox);

    FluxBox NtFforGrad;
    if (isMultiblock)
      {
        if (m_vectorFluxIntv.size() > 0)
          {
            // Seems very wasteful to calculate these on whole FluxBoxes,
            // when we need them only on block boundaries,
            // but that's what you do if you want to use functions that
            // take FluxBoxes.

            // NtFforGrad will have 1 tangential ghost.
            NtFforGrad.define(FBox1inDomain, m_numFluxes);
            thisCoordSys->
              computeMetricTermProductAverage(
                                              NtFforGrad,
                                              FfromWfaceAvg,
                                              N,
                                              SpaceDim, 
                                              FfromWfaceAvg, // unused
                                              FBox1inDomain,
                                              false, // 2nd-order
                                              NtFInterval, // restrict to m_vectorFluxIntv only?
                                              FInterval, // restrict to m_vectorFluxIntv only?
                                              0,
                                              &probDomain);

            // I would like to alias FluxBoxes so that we don't have to copy,
            // but FluxBox aliasing is not in Chombo yet.

            int numVectorFluxes = m_vectorFluxIntv.size();
            Interval vectorFluxIntv0(0, numVectorFluxes-1);

            FluxBox NtFVectorForGrad(FBox1inDomain, numVectorFluxes);
            NtFVectorForGrad.copy(NtFforGrad,
                                  m_vectorFluxIntv.begin(), 0,
                                  numVectorFluxes);

            // NtFVectorCen will have 0 ghosts.
            // NtFVectorCen := NtFAvg[m_vectorFluxIntv] -
            //    m_dx^2/24 * Laplacian(NtFVectorforGrad)
            FluxBox NtFVectorCen(FBox, numVectorFluxes);
            NtFVectorCen.copy(NtFAvg,
                              m_vectorFluxIntv.begin(), 0,
                              numVectorFluxes);
            MOLUtilities molUtil;
            molUtil.define(probDomain, m_dx);
            molUtil.deconvolveFace(NtFVectorCen, NtFVectorForGrad, FBox, -1);

            // NtFVectorForGrad has 1 tangential ghost.
            // Transform its vector components to orthonormal frame.
            // NtFVectorforGrad is second-order accurate, so is close enough to
            // be either face-centered or face-averaged.
            thisCoordSys->orthonormalize(NtFVectorForGrad, vectorFluxIntv0);

            // NtFVectorOrthoCen will have 0 ghosts.  It is face-CENTERED.
            // Copy from NtFVectorCen, and then
            // transform its vector components to orthonormal frame.
            FluxBox NtFVectorOrthoCen(FBox, numVectorFluxes);
            NtFVectorOrthoCen.copy(NtFVectorCen);
            thisCoordSys->orthonormalize(NtFVectorOrthoCen, vectorFluxIntv0);

            // NtFVectorOrthoAvg has 0 ghosts.
            // Copy from NtFVectorOrthoCen, and then
            // convert from face-CENTERED to face-AVERAGED.
            FluxBox NtFVectorOrthoAvg(FBox, numVectorFluxes);
            NtFVectorOrthoAvg.copy(NtFVectorOrthoCen);
            molUtil.deconvolveFace(NtFVectorOrthoAvg, NtFVectorForGrad, FBox, 1);

            // We will store the vector components of NtFAvg and NtFforGrad.
            // Save the face-averaged versions in the orthonormal frame.
            NtFAvg.copy(NtFVectorOrthoAvg,
                        0, m_vectorFluxIntv.begin(),
                        numVectorFluxes);

            NtFforGrad.copy(NtFVectorForGrad,
                            0, m_vectorFluxIntv.begin(),
                            numVectorFluxes);
          }

        // Update the block register.
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            for (SideIterator sit; sit.ok(); ++sit)
              {
                Side::LoHiSide side = sit();
                if (blockRegister.hasInterface(dit(), idir, side))
                  {
                    if (m_vectorFluxIntv.size() > 0)
                      {
                        blockRegister.storeAux(NtFforGrad[idir], dit(), idir, side);
                      }
                    // Very important:
                    // The vector flux components (if any) stored here
                    // are face-AVERAGED, and in the ORTHONORMAL FRAME.
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
// FIXME: Will this be correct on block boundaries?
    updateFluxTotalsAndRegisters(NtFAvg,
        a_finerFluxRegister, a_coarserFluxRegister,
        dit(), a_weight);

    // Find the RHS (= -div F)
    FArrayBox& LofUFab = a_LofJU[dit];
    m_patchMappedConsOperatorPtr->evalRHS(LofUFab, NtFAvg);
  }

  if (m_useSourceTerm)
    {
      // a_LofJU += <J * S(*m_UPtr)>
      m_sourceTermPtr->addSourceTerm(a_LofJU, *m_UPtr, m_dx*RealVect::Unit);
    }
}



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
  CH_TIME("LevelMappedConsOperator::definePatch");
  if (m_patchMappedConsOperatorPtr == NULL)
  {
    m_patchMappedConsOperatorPtr = new PatchMappedConsOperator(SpaceDim);
    // Copy the pointer into the base class PatchConsOperator, too
    m_patchConsOperatorPtr = (PatchConsOperator*) m_patchMappedConsOperatorPtr;
  }

// PatchMappedConsOperator m_patchMappedConsOperator;
  m_patchMappedConsOperatorPtr->define(m_domain, m_dx, m_levelGridMetricsPtr, a_molPhysics, m_numFields);
  m_patchMappedConsOperatorPtr->spaceOrder(m_spaceOrder);
  m_patchMappedConsOperatorPtr->limitFaceValues(m_limitFaceValues);
  m_patchMappedConsOperatorPtr->highOrderLimiter(m_highOrderLimiter);
  m_patchMappedConsOperatorPtr->useFlattening(m_useFlattening);
  m_patchMappedConsOperatorPtr->noPPM(m_noPPM);
  m_patchMappedConsOperatorPtr->doDeconvolution(m_doDeconvolution);
  m_patchMappedConsOperatorPtr->doFaceDeconvolution(m_doFaceDeconvolution);
  m_patchMappedConsOperatorPtr->useArtificialViscosity(m_useArtificialViscosity);
  m_patchMappedConsOperatorPtr->useArtificialDissipation(m_useArtificialDissipation);
  m_patchMappedConsOperatorPtr->artificialDissipation(m_artificialDissipation);
  m_patchMappedConsOperatorPtr->numGhost(m_numGhost);
  // Need give access to the unit normals
  if (m_levelGridMetricsPtr)
    {
      m_patchMappedConsOperatorPtr->unitNormals(&m_unitNormalLay);
    }
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
  CH_TIME("LevelMappedConsOperator::addArtificialViscosity");
  CH_assert(m_useArtificialViscosity);
  m_patchMappedConsOperatorPtr->setCurrentTime(a_oldTime);

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
      m_patchMappedConsOperatorPtr->addMappedArtificialViscosity(JUnewFab,
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
  CH_TIME("LevelMappedConsOperator::updateFluxTotalsAndRegisters");
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
  CH_TIME("LevelMappedConsOperator::fillGhostsRK4AndComputeU");

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
  m_patchMappedConsOperatorPtr->setCurrentBox(a_ind, curBox);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::setSourceTerm(LevelSourceTerm* a_sourceTermPtr)
{
  m_sourceTermPtr = a_sourceTermPtr;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::useArtificialDissipation(bool a_useArtificialDissipation)
{
  m_useArtificialDissipation = a_useArtificialDissipation;
  if (m_patchMappedConsOperatorPtr != NULL)
    {
      m_patchMappedConsOperatorPtr->useArtificialDissipation(m_useArtificialDissipation);
    }
}


//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::artificialDissipation(Real a_artificialDissipation)
{
  m_artificialDissipation = a_artificialDissipation;
  if (m_patchMappedConsOperatorPtr != NULL)
    {
      m_patchMappedConsOperatorPtr->artificialDissipation(m_artificialDissipation);
    }
}


/******************************************************************************/
/**
 * Member functions that are only used for multiblock calculations are defined
 * here.
 *//*+*************************************************************************/

/*--------------------------------------------------------------------*/
//  Define members required for multiblock grids
/** Such as m_facesToFill which is required for use of BlockRegister
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::defineMultiblockMbrs()
{
  CH_TIME("LevelMappedConsOperator::defineMultiblockMbrs");
//--Define m_facesToFill

  const MultiBlockCoordSys& coordSys = m_levelGridMetricsPtr->getCoordSys();
  const DisjointBoxLayout& grids = m_levelGridMetricsPtr->getBoxes();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSys.boundaries();
  // Tuple< LayoutData<IntVectSet>, 2*SpaceDim> m_facesToFill;
  int faceID = 0;
  for (SideIterator sit; sit.ok(); ++sit)
    {
      Side::LoHiSide side = sit();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_facesToFill[faceID].define(grids);
          for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
            {
              const Box& baseBox = grids[dit];
              int baseBlockNum = coordSys.whichBlock(baseBox);
              const Box& blockBox = coordSys.mappingBlocks()[baseBlockNum];
              Box blockBdryBox = adjCellBox( blockBox, idir, side, 1 );
              Box clipBox = adjCellBox( baseBox, idir, side, 1 );
              clipBox &= blockBdryBox;
              const BlockBoundary& bb = boundaries[baseBlockNum][faceID];
              if ( clipBox.isEmpty() || bb.isDomainBoundary() )
                { // empty or domain boundary
                  m_facesToFill[faceID][dit] = IntVectSet();
                }
              else
                {
                  const IndicesTransformation& it = bb.getTransformation();
                  Box clipOtherBox = it.transform(clipBox);
                  // ivsMissing will hold all IntVects of clipBox
                  // that are NOT present in the other block.
                  IntVectSet ivsMissing(clipBox);
                  // Loop through all boxes, and remove their back transforms
                  // from ivsMissing.
                  LayoutIterator lit = grids.layoutIterator();
                  for (lit.begin(); lit.ok(); ++lit)
                    {
                      const Box& baseOtherBox = grids[lit];
                      if (baseOtherBox.intersectsNotEmpty(clipOtherBox))
                        {
                          Box baseThisBox = it.transformBack(baseOtherBox);
                          ivsMissing -= baseThisBox;
                        }
                    }
                  // Got ivsMissing.
                  // ivsPresent will hold all IntVects of clipBox
                  // that ARE present in the other block.
                  IntVectSet ivsPresent(clipBox);
                  ivsPresent -= ivsMissing;
                  m_facesToFill[faceID][dit] = IntVectSet(ivsPresent);
                }
            }
          faceID++;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Sets single-valued fluxes at block boundaries
/** Using the contents of a BlockRegister, set the single-valued
 *  fluxes at block boundaries.
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::setCommonFlux(
  LevelData<FluxBox>&   a_flux,
  const BlockRegister&  a_blockRegister) const
{
  CH_TIME("LevelMappedConsOperator::setCommonFlux");
  const MultiBlockCoordSys& coordSys = m_levelGridMetricsPtr->getCoordSys();
  const DisjointBoxLayout& grids = m_levelGridMetricsPtr->getBoxes();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSys.boundaries();
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
  {
    FluxBox& FfaceAvg = a_flux[dit];
    const Box& baseBox = grids[dit];
    int baseBlockNum = coordSys.whichBlock(baseBox);
    const NewFourthOrderCoordSys* thisCoordSys =
      m_levelGridMetricsPtr->getCoordSys(baseBox);
    // const Box& FBox1inDomain = m_grow1inDomainLayout[dit];
    int faceID = 0;
    for (SideIterator sit; sit.ok(); ++sit)
    {
      Side::LoHiSide side = sit();
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (a_blockRegister.hasInterface(dit(), idir, side))
        {
          // maybe better if this is done inside BlockRegister
          const BlockBoundary& bb = boundaries[baseBlockNum][faceID];
          int reorientFace = bb.reorientFace(idir);
          Box faceBox = adjCellBox(baseBox, idir, side, 1);
          // faceBox.grow(faceGrowVect);
          // if Lo, then shift +1; if Hi, then shift -1
          faceBox.shiftHalf(idir, -sign(side));
          Side::LoHiSide sideOther = flip(side);
          // Need to define these FABs.
          FArrayBox fluxThisFab(faceBox, m_numFluxes); // from this block
          FArrayBox fluxOtherFab(faceBox, m_numFluxes); // from other block
          a_blockRegister.getFlux(fluxThisFab, dit(),
                                  idir, side, sideOther);
          a_blockRegister.getFlux(fluxOtherFab, dit(),
                                  idir, side, side);
          fluxOtherFab.mult(reorientFace * 0.5);
          fluxThisFab.mult(0.5);
          fluxThisFab += fluxOtherFab;
          // Now fluxThisFab contains mean of 
          // this face's flux and other face's flux;
          // these fluxes are all face-AVERAGED.
          if (m_vectorFluxIntv.size() > 0)
            {
              Box face1Box(faceBox);
              // Grow face1Box in tangential directions.
              // FIXME: Should really cut off face1Box at physical boundaries.
              face1Box.grow(IntVect::Unit - BASISV(idir));

              // auxiliary data contain ghosted fluxes in orthonormal frame,
              // good enough for taking laplacians.
              FArrayBox auxThisFab(face1Box, m_numFluxes); // from this block
              FArrayBox auxOtherFab(face1Box, m_numFluxes); // from other block
              a_blockRegister.getAux(auxThisFab, dit(),
                                     idir, side, sideOther);
              a_blockRegister.getAux(auxOtherFab, dit(),
                                     idir, side, side);

              auxOtherFab.mult(reorientFace * 0.5);
              auxThisFab.mult(0.5);
              auxThisFab += auxOtherFab;
              // Now auxThisFab contains mean of 
              // this face's aux and other face's aux;
              // these auxes are all face-AVERAGED.

              // Note that components m_vectorFluxIntv of fluxThisFab
              // are face-averaged fluxes in the ORTHONORMAL frame.

              // On the components m_vectorFluxIntv of fluxThisFab,
              // we must perform the following operations:
              // 1. convert from face-averaged to face-centered
              //    in the orthonormal frame;
              // 2. deorthonormalize (face-centered);
              // 3. convert from face-centered to face-averaged
              //    in the frame of this block.

              // 1. convert from face-averaged to face-centered
              //    in the orthonormal frame;
              // fluxThisFab -= Laplacian(auxThisFab) * m_dx^2 / 24
              // I wish I could use MOLUtilities::deconvolveFace(),
              // but that's for FluxBoxes.
              // I can't use MOLUtilities::deconvolve(), because that's for
              // cell-centered FABs and takes derivatives in all directions.
              // And FourthOrderUtil functions do not take the additional
              // "for grad" argument.
              // FIXME: This won't work when there are physical boundaries.

              FArrayBox fluxVectorThisFab, auxVectorThisFab; // aliases
              fluxVectorThisFab.define(m_vectorFluxIntv, fluxThisFab);
              auxVectorThisFab.define(m_vectorFluxIntv, auxThisFab);
              FArrayBox d2auxVectorThisFab(faceBox, SpaceDim);
              Real d2scaling = 1./24.;

              d2auxVectorThisFab.setVal(0.);
              for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
                {
                  if (tanDir != idir)
                    {
                      // With e = unit vector in direction tanDir,
                      // d2auxVectorThisFab[i] :=
                      // auxVectorThisFab[i - e]
                      // - 2*auxVectorThisFab[i]
                      // + auxVectorThisFab[i + e].
                      auxVectorThisFab.shift(tanDir, -1);
                      d2auxVectorThisFab.plus(auxVectorThisFab, 1.);

                      auxVectorThisFab.shift(tanDir, +1);
                      d2auxVectorThisFab.plus(auxVectorThisFab, -2.);

                      auxVectorThisFab.shift(tanDir, +1);
                      d2auxVectorThisFab.plus(auxVectorThisFab, 1.);

                      // shift auxVectorThisFab back to where it started
                      auxVectorThisFab.shift(tanDir, -1);
                    }
                }
              fluxVectorThisFab.plus(d2auxVectorThisFab, -d2scaling);

              // 2. deorthonormalize (face-centered);
              thisCoordSys->deorthonormalize(auxThisFab, m_vectorFluxIntv);
              thisCoordSys->deorthonormalize(fluxThisFab, m_vectorFluxIntv);

              // 3. convert from face-centered to face-averaged
              //    in the frame of this block.
              // fluxThisFab += Laplacian(auxThisFab) * m_dx^2 / 24

              d2auxVectorThisFab.setVal(0.);
              for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
                {
                  if (tanDir != idir)
                    {
                      // With e = unit vector in direction tanDir,
                      // d2auxVectorThisFab[i] :=
                      // auxVectorThisFab[i - e]
                      // - 2*auxVectorThisFab[i]
                      // + auxVectorThisFab[i + e].
                      auxVectorThisFab.shift(tanDir, -1);
                      d2auxVectorThisFab.plus(auxVectorThisFab, 1.);

                      auxVectorThisFab.shift(tanDir, +1);
                      d2auxVectorThisFab.plus(auxVectorThisFab, -2.);

                      auxVectorThisFab.shift(tanDir, +1);
                      d2auxVectorThisFab.plus(auxVectorThisFab, 1.);

                      // shift auxVectorThisFab back to where it started
                      auxVectorThisFab.shift(tanDir, -1);
                    }
                }
              fluxVectorThisFab.plus(d2auxVectorThisFab, d2scaling);
            }

          // FfaceAvg[idir].copy(fluxThisFab);
          FArrayBox& faceFillFab = FfaceAvg[idir];
          int ncomp = faceFillFab.nComp();
          const IntVectSet& ivs = m_facesToFill[faceID][dit];
          for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
            {
              // ivClip is index of the ghost cell just outside (side, idir)
              // of the block.
              // (Store cells, not faces, because IntVectSet is for cells.)
              IntVect ivClip = ivsit();
              // ivFace is index of the cell face on the boundary,
              // or ivClip shifted by "half" in (-side, idir) direction.
              // If side is high, then ivFace is ivClip shifted down from
              // cell to face in direction idir.  So no change in coordinates.
              // If side is low, then ivFace is ivClip shifted up from
              // cell to face in direction idir.  So +1 in idir coordinate.
              IntVect ivFace = ivClip;
              if (side == Side::Lo) ivFace += BASISV(idir);
              for (int icomp = 0; icomp < ncomp; icomp++)
                {
                  faceFillFab(ivFace, icomp) = fluxThisFab(ivFace, icomp);
                }
            }
        }
        faceID++;
      } // iterate over dimensions
    } // iterate over sides
  } // iterate over patches
}

//////////////////////////////////////////////////////////////////////////////

Real
LevelMappedConsOperator::getMaxSpeed()
{
  CH_TIME("LevelMappedConsOperator::getMaxSpeed");
  const DisjointBoxLayout& grids = m_levelGridMetricsPtr->getBoxes();
  Real speedMaxLocal = 0.;
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());
      Box box = grids[dit];
      FArrayBox& fabU = (*m_UPtr)[dit];
      // kludge because this class does not save MOLMappedPhysics
      Real speed = m_patchMappedConsOperatorPtr->getMaxWaveSpeed(fabU, box);
      speedMaxLocal = Max(speed, speedMaxLocal);
    }
  Real speedMax = speedMaxLocal;
#ifdef CH_MPI
  MPI_Allreduce(&speedMaxLocal, &speedMax, 1,
                MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
#endif
  return speedMax;
}

//////////////////////////////////////////////////////////////////////////////

#include "NamespaceFooter.H"
