#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


/******************************************************************************/
/**
 * \file LevelGridMetrics.cpp
 *
 * \brief Non-inline definitions for classes in LevelGridMetrics.H
 *
 *//*+*************************************************************************/

#include <limits>

#include "LevelGridMetrics.H"
#include "BoxIterator.H"
#include "CoarseAverage.H"
#include "PiecewiseLinearFillPatch.H"
#include "FourthOrderUtil.H"
#include "TimeInterpolatorRK4.H"
#include "BlockRegister.H"

#include "NamespaceHeader.H"

// //--Debug support

// #include "DebugOut.H"
// int LGMDebug = 0;


/*******************************************************************************
 *
 * Class IntersectionAverageFace: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/** The constructor performs the average and stores the intermediate
 *  averaged result in the averager.
 *  \param[in]  a_A     The fine data to be averaged
 *  \param[in]  a_nRef  The refinement ratio between the fine and
 *                      coarse levels
 *//*-----------------------------------------------------------------*/

IntersectionAverageFace::IntersectionAverageFace(
  const LevelData<FluxBox>& a_A,
  const int                 a_nRef)
  :
  m_nComp(a_A.nComp()),
  m_nRef(a_nRef)
{
  const int nComp = a_A.nComp();
  const DisjointBoxLayout& dblA = a_A.getBoxes();
  coarsen(m_crDblA, dblA, a_nRef);
  m_averager.define(dblA, nComp, a_nRef);
  m_averager.average(a_A);
}

/*--------------------------------------------------------------------*/
//  Copy the average using a mask
/** The average is copied to the layout of the coarse mesh (E) and
 *  then copied to (C) only at the intersection of Cr(A^(l+1)),
 *  Cr(B^(l+1)), and C^(l)
 *  \param[out] a_C     The fine data to be averaged
 *  \param[in]  a_dblB  The disjoint layout of B
 *//*-----------------------------------------------------------------*/

void
IntersectionAverageFace::copyTo(LevelData<FluxBox>&      a_C,
                                const DisjointBoxLayout& a_dblB)
{
  CH_TIME("IntersectionAverageFace::copyTo");
  DisjointBoxLayout crDblB;
  coarsen(crDblB, a_dblB, m_nRef);
  const DisjointBoxLayout& dblC = a_C.getBoxes();

  // Needs at least 1 ghost cell in E for the copy to work
  IntVect reqGhostVect = a_C.ghostVect();
  if (reqGhostVect.product() == 0)
    {
      reqGhostVect = IntVect::Unit;
    }
  LevelData<FluxBox> E(dblC, m_nComp, reqGhostVect);
  m_averager.copyTo(E);

  // Now compute the mask and copy from E to a_C based on the mask
  LayoutIterator litA = m_crDblA.layoutIterator();
  LayoutIterator litB = crDblB.layoutIterator();
  DataIterator ditC = dblC.dataIterator();
  for (ditC.begin(); ditC.ok(); ++ditC)
    {
      Box box = dblC[ditC];
      // Because we will eventually switch to faces (i.e.,
      /// surroundingNodes), grow the box on the high-side by 1
      IntVectSet maskA(Box(box.smallEnd(), box.bigEnd() + 1));
      maskA.makeEmptyBits();
      IntVectSet maskB(maskA);
      const int boxBigEnd0 = box.bigEnd(0);
      for (litA.begin(); litA.ok(); ++litA)
        {
          Box otherBox = m_crDblA[litA];
          if (otherBox.smallEnd(0) > boxBigEnd0)
            {
              // Can break out of loop, since we know that the smallEnd of
              // all the remaining boxes are lexigraphically beyond this
              // ghost box.
              break;
            }
          maskA |= (box & otherBox);
        }
      if (maskA.isEmpty())
        {
          // If maskA is empty then no intersections so continue with next
          // box in C
          continue;
        }
      for (litB.begin(); litB.ok(); ++litB)
        {
          Box otherBox = crDblB[litB];
          if (otherBox.smallEnd(0) > boxBigEnd0)
            {
              // Can break out of loop, since we know that the smallEnd of
              // all the remaining boxes are lexigraphically beyond this
              // ghosted box.
              break;
            }
          maskB |= (box & otherBox);
        }
      if (maskB.isEmpty())
        {
          // If maskB is empty then no intersections so continue with next
          // box in C
          continue;
        }
      maskA &= maskB;

      // Now for each direction (maybe use Fortran?)
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          FArrayBox& fabC = a_C[ditC][idir];
          FArrayBox& fabE = E[ditC][idir];
          // Switch the IVS to face-based in direction 'idir' using 'maskB'.
          // This is equivalent to adding IntVects to the high-side of each
          // existing IntVect.
          maskB = maskA;
          maskB.growHi(idir);
          IVSIterator ivs(maskB);
          for (ivs.begin(); ivs.ok(); ++ivs)
            {
              IntVect iv = ivs();
              fabC(iv) = fabE(iv);
            }
        }
    }
}


/*******************************************************************************
 *
 * Class LevelGridMetrics: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_nComp Number of components in \f$<U>\f$ and
 *                      \f$<JU>\f$
 *  \param[in]  a_spaceOrder
 *                      Pretends that you can set the spatial order.
 *                      This is where you learn it must be 4. The
 *                      multiblock exchanger asks for this, hence the
 *                      charade.
 *//*-----------------------------------------------------------------*/

LevelGridMetrics::LevelGridMetrics(const int a_nComp,
                                   const int a_spaceOrder)
  :
  m_scrNAverager(NULL),
  m_NAverager(NULL),
  m_NtFDiffAverager(NULL),
  m_coordSys(NULL),
  m_finerLevelGridMetrics(NULL),
  m_coarserLevelGridMetrics(NULL),
  m_spaceOrder(4),
  m_interpolator(NULL),
  m_cachedBlockIdx(-2),
  m_cachedBlockNG(-1),
  m_transverseNOpType(TransverseNOpNone),
  m_metricsChanged(false),
  m_metricsDefined(false),
  m_useScriptN(true),
  m_isMultiblock(false),
  m_haveMultiBlockVectorData(false)
{
  CH_TIME("LevelGridMetrics::LevelGridMetrics");
  m_interpolator = new FourthOrderMappedFineInterp(3, a_nComp);
  if (a_spaceOrder != 4)
    {
      MayDay::Error("LevelGridMetrics only supports fourth-order in space");
    }
  for (int iAge = 0; iAge != AgeNum; ++iAge)
    {
      m_NtXAverager[iAge] = NULL;
    }
  // We access some static members of the CodimBox<FArrayBox>  class before it may be
  // automatically intiallized so initialize it now.
  CodimBox<FArrayBox> ::initialize();
}

/*--------------------------------------------------------------------*/
//  Destructor
/*--------------------------------------------------------------------*/

LevelGridMetrics::~LevelGridMetrics()
{
  delete m_interpolator;
  clear();
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Free memory
/*--------------------------------------------------------------------*/

void LevelGridMetrics::clear()
{
  if (m_coordSys)
    {
      delete m_coordSys;
    }
  m_coordSys = NULL;
}

/*--------------------------------------------------------------------*/
//  Does a finer level exist?
/** Only if there is a pointer and the metrics have been defined
 *  \note
 *  <ul>
 *    <li> hasFiner() is not ready for use until after the metrics
 *         have been defined (usually in postInitialGrid())
 *  </ul>
 *//*-----------------------------------------------------------------*/

bool
LevelGridMetrics::hasFiner() const
{
  return ((m_finerLevelGridMetrics != NULL) && 
          m_finerLevelGridMetrics->m_metricsDefined);
}

/*--------------------------------------------------------------------*/
//  Does a coarser level exist?
/** Sufficient that there is a pointer
 *  \note
 *  <ul>
 *    <li> hasCoarser() ready for use after define() is called.
 *  </ul>
 *//*-----------------------------------------------------------------*/

bool
LevelGridMetrics::hasCoarser() const
{
  return (m_coarserLevelGridMetrics);
}

/*--------------------------------------------------------------------*/
//  Define the class.  No metrics terms are yet available.
/** \param[in]  a_AMRLevel
 *                      Parent AMR level
 *  \param[in]  a_coordSysFact
 *                      Factory for constructing the multi-block
 *                      coordinate system.
 *  \param[out] a_coarserLevelGridMetrics
 *                      Modified to point to us as a finer level of
 *                      grid metrics
 *  \param[in]  a_dxVect
 *                      Computational mesh spacing in each direction
 *  \param[in]  a_ghostVect
 *                      Number of ghost cells
 *  \param[in]  a_transverseNOpType
 *                      Method for determining transverse components
 *                      of N
 *  \param[in]  a_numBlocks
 *                      Blocks have different coordinate systems
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::define(
  AMRLevel *const                  a_AMRLevel,
  MultiBlockCoordSysFactory *const a_coordSysFact,
  LevelGridMetrics *const          a_coarserLevelGridMetrics,
  const RealVect&                  a_dxVect,
  const IntVect&                   a_ghostVect,
  const TransverseNOpType          a_transverseNOpType,
  const bool                       a_haveMultiBlockVectorData,
  const int                        a_numBlocks)
{
  // Allocate the multi-block coordinate system.
  MultiBlockCoordSys *const coordSys =
    a_coordSysFact->getCoordSys(a_AMRLevel->problemDomain(), a_dxVect);
  define(a_AMRLevel,
         coordSys,
         a_coarserLevelGridMetrics,
         a_dxVect,
         a_ghostVect,
         a_transverseNOpType,
         a_haveMultiBlockVectorData,
         a_numBlocks);
}

/*--------------------------------------------------------------------*/
//  Define the class.  No metrics terms are yet available.
/** \param[in]  a_AMRLevel
 *                      Parent AMR level
 *  \param[in]  a_coordSysFact
 *                      Factory for constructing the multi-block
 *                      coordinate system.
 *  \param[out] a_coarserLevelGridMetrics
 *                      Modified to point to us as a finer level of
 *                      grid metrics
 *  \param[in]  a_dxVect
 *                      Computational mesh spacing in each direction
 *  \param[in]  a_ghostVect
 *                      Number of ghost cells
 *  \param[in]  a_transverseNOpType
 *                      Method for determining transverse components
 *                      of N
 *  \param[in]  a_numBlocks
 *                      Blocks have different coordinate systems
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::define(
  AMRLevel *const           a_AMRLevel,
  MultiBlockCoordSys *const a_coordSys,
  LevelGridMetrics *const   a_coarserLevelGridMetrics,
  const RealVect&           a_dxVect,
  const IntVect&            a_ghostVect,
  const TransverseNOpType   a_transverseNOpType,
  const bool                a_haveMultiBlockVectorData,
  const int                 a_numBlocks)
{
  CH_TIME("LevelGridMetrics::define");
  clear();  // Free memory from previous define
  m_parentAMRLevel = a_AMRLevel;
  m_ghostVect = a_ghostVect;
  m_transverseNOpType = a_transverseNOpType;
  m_haveMultiBlockVectorData = a_haveMultiBlockVectorData;
  m_dxVect = a_dxVect;

  if (verbosity() >= 3)
    {
      pout() << "LevelGridMetrics::define " << level() << std::endl;
    }

  // Save the multi-block coordinate system.
  m_coordSys = a_coordSys;

  m_useScriptN = m_coordSys->useScriptN();
  m_isMultiblock = (m_coordSys->numBlocks() > 1);
  m_cachedBlockIdx = -2;
  m_cachedBlockNG  = -1;
  if (!isMultiblock())
    {
      m_cachedBlockDomain = problemDomain();
    }

  // Set the finer pointer in the coarser level to point to us
  if (a_coarserLevelGridMetrics)
    {
      m_coarserLevelGridMetrics = a_coarserLevelGridMetrics;
      a_coarserLevelGridMetrics->m_finerLevelGridMetrics = this;
    }

  // Re-definition of the class forces all metrics to be undefined
  m_metricsChanged = false;
  m_metricsDefined = false;
}

/*--------------------------------------------------------------------*/
//  After a time step, average down \f$N\f$ and \f$J\f$
/** This is only performed if the mesh changed on a finer level.
 *  Since the boundary of the finer mesh did not change so all we have
 *  to do is averaging.  If the metrics changed on the finer mesh,
 *  average down \f$\mathbf(N)\f$ and J.
 *  \note
 *  <ul>
 *    <li> Transverse terms in N are also averaged.
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::postTimeStep()
{
  CH_TIME("LevelGridMetrics::postTimeStep");
  if (verbosity() >= 3)
    {
      pout() << "LevelGridMetrics::postTimeStep " << level();
    }

  if (hasFiner())
    {
      if (didFinerChange())
        {
          if (verbosity() >= 3)
            {
              pout() << ":  averaging metrics";
            }

          // Average N
          const int numNTerms = SpaceDim*SpaceDim;
          CoarseAverageFace NAverager(m_finerLevelGridMetrics->getBoxes(),
                                      numNTerms,
                                      refRatio());
          NAverager.averageToCoarse(m_N,
                                    m_finerLevelGridMetrics->m_N);
          // Average J
          CoarseAverage JAverager(m_finerLevelGridMetrics->getBoxes(),
                                  1,
                                  refRatio());
          JAverager.averageToCoarse(m_J,
                                    m_finerLevelGridMetrics->m_J);
          // Propagate the change flags
          m_metricsChanged = true;
          m_finerLevelGridMetrics->m_metricsChanged = false;
        }
    }

  if (verbosity() >= 3)
    {
      pout() << std::endl;
    }
}

/*--------------------------------------------------------------------*/
//  Set up initial grid
/** The option exists to set the new grid now or in postInitialGrid.
 *  \param[in]  a_newGrid
 *                      The layout of the new boxes.  Use NULL to
 *                      avoid defining the mesh now.
 *  \note
 *  <ul>
 *    <li> Note that this routine traverses coarse to fine
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::initialGrid(const DisjointBoxLayout *const a_newGrid)
{
  CH_TIME("LevelGridMetrics::initialGrid");
  if (verbosity() >= 3)
    {
      pout() << "LevelGridMetrics::initialGrid " << level() << std::endl;
    }

  // Define the mesh
  if (a_newGrid != NULL)
    {
      defineGrids(*a_newGrid);
      defineGridsNP(*a_newGrid);

      // The space interpolator lives on the coarser level.  Only define if
      // the grid on this level has more than zero boxes (otherwise it is
      // unused).
      if (hasCoarser() && a_newGrid->size() > 0)
        {
          const int cr2ThisRefRatio = m_coarserLevelGridMetrics->refRatio();
          IntVect refRatioVec(D_DECL6(cr2ThisRefRatio,
                                      cr2ThisRefRatio,
                                      cr2ThisRefRatio,
                                      cr2ThisRefRatio,
                                      cr2ThisRefRatio,
                                      cr2ThisRefRatio));
          m_coarserLevelGridMetrics->m_interpolator->defineStencils(
            refRatioVec,
            dxVect());
          const MultiBlockLevelGeom *const mbgeoPtr =
            (isMultiblock()) ? &(m_coarserLevelGridMetrics->m_mbgeo) : 0;
          m_coarserLevelGridMetrics->m_interpolator->defineGrids(
            getBoxes(),
            m_coarserLevelGridMetrics->getBoxes(),
            m_coarserLevelGridMetrics->problemDomain(),
            m_ghostVect,
            true,
            mbgeoPtr);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Set up the initial metric terms.
/** postInitialGrid has a fine to coarse traversal.
 *  \param[in]  a_newGrid
 *                      The layout of the new boxes.  Use NULL to
 *                      avoid defining the mesh now but it must have
 *                      been set somewhere else.
 *  \note
 *  <ul>
 *    <li> Note that this routine traverses fine to coarse
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::postInitialGrid(const DisjointBoxLayout *const a_newGrid)
{
  CH_TIME("LevelGridMetrics::postInitialGrid");
  if (verbosity() >= 3)
    {
      pout() << "LevelGridMetrics::postInitialGrid " << level() << std::endl;
    }

  // Define the mesh
  if (a_newGrid != NULL)
    {
      // We expect a grid with > 0 boxes in postInitialGrid
      CH_assert(a_newGrid->size() > 0);
      defineGrids(*a_newGrid);
      defineGridsNP(*a_newGrid);
    }

  // Define the metrics for this level
  defineMetrics();

  // Define the interpolator if a finer level exists
  if (hasFiner())
    {
      if (a_newGrid != NULL)
        {
          // The space interpolator lives on the coarser level
          const int this2FnRefRatio = refRatio();
          IntVect refRatioVec(D_DECL6(this2FnRefRatio,
                                      this2FnRefRatio,
                                      this2FnRefRatio,
                                      this2FnRefRatio,
                                      this2FnRefRatio,
                                      this2FnRefRatio));
          m_interpolator->defineStencils(
            refRatioVec,
            dxVect());
          const MultiBlockLevelGeom *const mbgeoPtr =
            (isMultiblock()) ? &m_mbgeo : 0;
          m_interpolator->defineGrids(
            m_finerLevelGridMetrics->getBoxes(),
            getBoxes(),
            problemDomain(),
            m_finerLevelGridMetrics->m_ghostVect,
            true,
            mbgeoPtr);
        }
      m_interpolator->defineCrFnJ(m_J);
      m_interpolator->defineCFInterface(this);
    }
}

/*--------------------------------------------------------------------*/
//  Compute new metrics and correct the solution on the coarser meshes
/** Pre-regrid has a fine to coarse traversal so we can immediately
 *  average down the metrics.  This performs "snap-back" corrections
 *  in preparation for interpolation.  The level being corrected is
 *  unannotated, i.e., a_UOld; a_CrUOldPtr is the next coarser level
 *  and there must also be a finer level for there to be any effect.
 *
 *  \param[in]  a_baseLevel
 *                      The base level for the regrid
 *  \param[in]  a_newGrid
 *                      The layout of the new boxes
 *  \param[in]  a_CrUOldPtr
 *                      The physical solution on the old mesh at the next
 *                      coarser level.
 *  \param[in]  a_UOld  The physical solution on the old mesh.  This must have
 *                      valid interior ghosts filled by exchange.
 *  \param[out] a_UOld  One layer of invalid interior ghosts are filled by
 *                      piecewise constant interpolation and one layer of domain
 *                      boundary ghosts are filled by zeroth order
 *                      extrapolation.
 *  \param[in]  a_JUOld \<JU\> on the old mesh
 *  \param[out] a_JUOld \<JU\> still on the old mesh but corrected to
 *                      the new definition of the metrics
 *                      handled within.
 *  \note
 *  <ul>
 *    <li> Upon exit, all data members are updated to the new mesh.
 *    <li> Beware of the modification to U.  The next use of U is expected
 *         to be for interpolation from the coarse to fine mesh (as in
 *         ::regrid).  Because of proper nesting and boundary stencils
 *         the modifications made here will not affect the interpolation.
 *  </ul>
 */
/*
//  *   \param[out] a_JOld if (a_claimJOld = T), this is \<J\> on the old
//  *                      mesh but corrected to the new definition of
//  *                      the metrics.  Upon exit, 'm_J' in this class
//  *                      will contain \<J\> for the new mesh on this
//  *                      level.  The calling routine must take
//  *                      responsibility for deleting a_JOld.
//  *                      if (a_claimJOld = F), not modified.
//  *   \param[in]  a_claimJOld
//  *                      T - the calling routine will take control of
//  *                          the memory allocated to m_J.
//  *                      F - a_JOld is unmodified and memory for m_J is
//  *
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::preRegrid(const int                         a_baseLevel,
                            const DisjointBoxLayout&          a_newGrid,
                            const LevelData<FArrayBox> *const a_CrUOldPtr,
                            LevelData<FArrayBox>&             a_UOld,
                            LevelData<FArrayBox>&             a_JUOld)
{
  CH_TIME("LevelGridMetrics::preRegrid");
  if (verbosity() >= 3)
    {
      pout() << "LevelGridMetrics::preRegrid " << level() << std::endl;
    }

  const bool isBase = (level() == a_baseLevel);

  // If there is a finer level but the average volume flux from that level has
  // not been averaged, then for some reason the volume flux was not computed on
  // an old finer level.  This is probably because after re-grid this older
  // finer level will cease to exist.  Traverse any finer levels to get a
  // cascade of average volume flux from the old mesh.
  bool undefFinerLevels = false;
  if (hasFiner() && !m_NtXAverager[AgeNew])
    {
      m_finerLevelGridMetrics->volFluxCascade();
      // At the end of this routine, mark these finer levels as undefined.
      undefFinerLevels = true;
    }

  LevelData<FluxBox> NtXOldLev;
  LevelData<FluxBox> NtXNewLev;

  // Recompute the old volume flux
  if (m_metricsDefined)
    {
      NtXOldLev.define(m_gridsNP, 1, getJReqGhostVect());
      getVolFlux(NtXOldLev, m_N); // NtXOldLev is contiguous
      averageVolFlux(isBase, NtXOldLev, AgeOld); // needs contiguous
    }
  else
    {
      // This must be a new level.  Compute all metric terms.  This will
      // allocate averageVolFlux(NtXNewLev, AgeNew) if a coarser level exists.
      defineGrids(a_newGrid);
      defineGridsNP(a_newGrid);
      defineMetrics();
      // Bail before getting into any of the snapback stuff.
      return;
    }

  // Convert the stored non-periodic grids to the new mesh (the old NP grids
  // should be unused from here on).
  if (!isBase)  // No change if base
    {
      defineGridsNP(a_newGrid);
    }

  if (m_useScriptN)
    {

//--N is calculated using script N

      // Update script N to the new mesh
      LevelData<CodimBox<FArrayBox> > scrN(m_gridsNP,
                               SpaceDim,
                               getNReqGhostVect(),
                               CodimBoxFactory<FArrayBox> (2));
      getScriptN(scrN);
      averageScriptN(isBase, scrN);

      // Update N to the new mesh
      if (!isBase)  // Grid didn't change for the base
        {
          const int numNTerms = SpaceDim*SpaceDim;
          m_N.define(a_newGrid, numNTerms, getNReqGhostVect());
        }
      faceNormalN(m_N, scrN);
      // Option for transverse terms
      // If averaging, transverse components are calculated from normal
      // components and you should exchange first (note that if averaging,
      // transverse components are only computed on the faces of valid cells).
      // Otherwise exchange after.
      if (m_transverseNOpType != TransverseNOpAverage)
        {
          faceTransverseN(m_N);
        }
      //**FIXME This should only exchange across periodic boundaries (not
      //**everywhere)
      if (problemDomain().isPeriodic())
        {
          m_N.exchange();
        }
      if (m_transverseNOpType == TransverseNOpAverage)
        {
          faceTransverseN(m_N);
        }
    }
  else
    {

//--N is calculated directly (not freestream preserving)

      // Update N to the new mesh
      if (!isBase)  // Grid didn't change for the base
        {
          const int numNTerms = SpaceDim*SpaceDim;
          m_N.define(a_newGrid, numNTerms, getNReqGhostVect());
        }
      getN(m_N);
      averageN(isBase, m_N);
    }

  // Get the new volume flux (m_gridsNP now has new grids)
  NtXNewLev.define(m_gridsNP, 1, getJReqGhostVect());
  getVolFlux(NtXNewLev, m_N); // NtXNewLev contiguous
  averageVolFlux(isBase, NtXNewLev, AgeNew); // needs contiguous

  // added by petermc, 21 Nov 2011
  // Make NtXOldLev and NtXNewLev NON-contiguous in cubed-sphere.
  m_coordSys->separateVolFlux(NtXOldLev);
  m_coordSys->separateVolFlux(NtXNewLev);

  // Snapback calculations are not required on the finest mesh and are only
  // performed if metrics are defined using scriptN.  Current examples (such
  // as the mapped climate code) use an exact analytic definition of N and
  // hence do not require corrections.
  if (hasFiner() && m_useScriptN)
    {
      // Begin computation of the difference in flux, or "snapback" flux.
      // Start with the old volume flux overwritten by the new volume flux
      // on the old mesh.
      LevelData<FluxBox> NtXDiffLev; // should be NON-contiguous in cubed-sphere
      NtXDiffLev.define(NtXOldLev);
      NtXNewLev.copyTo(NtXDiffLev);

      const int nUComp = a_UOld.nComp();
      LevelData<FluxBox> NtFDiffLev;
      NtFDiffLev.define(getBoxes(), nUComp, IntVect::Zero);

      // Set up invalid ghost cells in a_UOld
      if (hasCoarser() && a_CrUOldPtr)
        {
          // Piecewise-constant interpolation for interior ghosts
          // FIXME: problemDomain won't be right for multiblock.
          PiecewiseLinearFillPatch plfp(
            a_UOld.getBoxes(),
            a_CrUOldPtr->getBoxes(),
            nUComp,
            m_coarserLevelGridMetrics->problemDomain(),
            m_coarserLevelGridMetrics->refRatio(),
            1,
            true);
          plfp.fillInterpPWConstSpace(a_UOld,
                                      *a_CrUOldPtr,
                                      0,
                                      0,
                                      nUComp);
          // Exterior ghosts (updated in loop below)
        }

      DataIterator ditOld = getDataIterator();
      for (ditOld.begin(); ditOld.ok(); ++ditOld)
        {
          // Only working on physical cells.
          const Box box = getBoxes()[ditOld];

          // Subtract the old volume flux to complete the snapback volume flux
          FluxBox& NtXDiff = NtXDiffLev[ditOld];
          NtXDiff.minus(NtXOldLev[ditOld], box, 0, 0, 1);

          FArrayBox& UOld = a_UOld[ditOld];
          // Zeroth-order extrapolation for exterior ghosts.
          // This means that on each cell just outside the block domain,
          // UOld is set to the value on the cell just inside the block domain.
          zerothOrderCellExtrapAtDomainBdry(UOld,
                                            box,
                                            blockDomain(box)); // 0 ghosts

          // BEGIN DEBUG
          // int blockNum = m_coordSys->whichBlock(box);
          // Box blockBox = m_coordSys->mappingBlocks()[blockNum];
          // END DEBUG

          // Now compute the snapback solution flux
          FluxBox& NtFDiff = NtFDiffLev[ditOld];
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              FArrayBox& NtXDiffDir = NtXDiff[idir];
              FArrayBox& NtFDiffDir = NtFDiff[idir];
              const Box fbox = NtFDiffDir.box(); // FACE-centered
              BoxIterator bit(fbox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect ivFace(bit());
                  Real NtXDiffVal = NtXDiffDir(ivFace);

                  // Set ivCell to a cell adjacent to ivFace:
                  // if NtXDiffVal > 0. then above; if < 0. then below.
                  IntVect ivCell(ivFace);
                  if (NtXDiffVal < 0.)
                    {
                      ivCell.shift(idir, -1);
                    }

                  for (int comp = 0; comp < nUComp; ++comp)
                    {
                      NtFDiffDir(ivFace, comp) = -NtXDiffVal*UOld(ivCell, comp);
                    }
                }  // Loop over fbox
            }  // Loop over dimensions
        }  // Loop over patches

      if (isMultiblock())
        {
          // Copy the snapback metrics into a block boundary register
          BlockRegister blockRegNtFDiff(m_coordSys, getBoxes(), 0); // 0 ghosts
          for (ditOld.begin(); ditOld.ok(); ++ditOld)
            {
              FluxBox& NtFDiff = NtFDiffLev[ditOld];
              for (int idir = 0; idir < SpaceDim; ++idir)
                {
                  FArrayBox& NtFDiffDir = NtFDiff[idir];
                  // Update the block register.
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      Side::LoHiSide side = sit();
                      if (blockRegNtFDiff.hasInterface(ditOld(), idir, side))
                        {
                          blockRegNtFDiff.storeFlux(
                            NtFDiffDir, ditOld(), idir, side);
                        }
                    }
                }
            }

          blockRegNtFDiff.close();

          // Make NtFDiffLev single-valued on block boundaries.
          for (ditOld.begin(); ditOld.ok(); ++ditOld)
            {
              FluxBox& NtFDiff = NtFDiffLev[ditOld];
              const FluxBox& NtXDiff = NtXDiffLev[ditOld];
              Box baseBox = getBoxes()[ditOld];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  const FArrayBox& NtXDiffDir = NtXDiff[idir];
                  FArrayBox& NtFDiffDir = NtFDiff[idir];

                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      Side::LoHiSide side = sit();
                      if (blockRegNtFDiff.hasInterface(ditOld(), idir, side))
                        {
                          // Low face:  Whenever NtXDiffDir < 0., reach across
                          // and set NtFDiffDir from the other block.
                          // High face:  Whenever NtXDiffDir > 0., reach across
                          // and set NtFDiffDir from the other block.
                          Box faceBox = adjCellBox(baseBox, idir, side, 1);
                          faceBox.shiftHalf(idir, -sign(side));
                          FArrayBox fluxFab(faceBox, nUComp);  // From other block
                          blockRegNtFDiff.getFlux(
                            fluxFab, ditOld(), idir, side, side);
                          BoxIterator bit(faceBox);
                          for (bit.begin(); bit.ok(); ++bit)
                            {
                              IntVect ivFace(bit());
                              if (NtXDiffDir(ivFace, 0)*sign(side) > 0.)
                                for (int comp = 0; comp < nUComp; ++comp)
                                  {
                                    NtFDiffDir(ivFace, comp) =
                                      fluxFab(ivFace, comp);
                                  }
                            }
                        }
                    }
                }  // Loop over dimensions
            }  // Loop over patches
        }  // If multiblock

      // Update the snapback solution flux with averages from finer levels.
      averageSnapbackSolFlux(isBase, NtFDiffLev);

      /*
        Update
        m_J += div(NtXDiffLev)
        a_JUOld -= div(NtFDiffLev)
        on each cell, where div is computed by face differences.
      */
      for (ditOld.begin(); ditOld.ok(); ++ditOld)
        {
          const FluxBox& NtXDiff = NtXDiffLev[ditOld];
          const FluxBox& NtFDiff = NtFDiffLev[ditOld];
          // Only working on physical cells.
          const Box box = getBoxes()[ditOld];

          //**FIXME only update J if the user needs the old mesh
          // Update J and JU with the snapback fluxes
          FArrayBox& JUOld = a_JUOld[ditOld];
          FArrayBox& JOld  = m_J[ditOld];
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              const FArrayBox& NtXDiffDir = NtXDiff[idir];
              const FArrayBox& NtFDiffDir = NtFDiff[idir];
              Real dxinv = 1./m_dxVect[idir];
              BoxIterator bit(box);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit(); // cell
                  IntVect ivLo = iv; // idir-face below cell iv
                  IntVect ivHi = iv + BASISV(idir); // idir-face above cell iv
                  JOld(iv) += dxinv *
                    ( NtXDiffDir(ivHi) - NtXDiffDir(ivLo) );
                  for (int comp = 0; comp < nUComp; ++comp)
                    {
                      JUOld(iv, comp) += dxinv *
                        ( NtFDiffDir(ivLo, comp) - NtFDiffDir(ivHi, comp) );
                    }
                }
            }
        }

    }  // Finest level bail out

  // Update J to the new mesh.
  //**FIXME Does the user need the corrected old mesh?
  if (!isBase)
    {
      defineGrids(a_newGrid);
      // m_gridsNP already defined
      m_J.define(getBoxes(), 1, getJReqGhostVect());
    }

  // Update m_J
  setAvgJfromNtX(NtXNewLev);

  //**FIXME This should only exchange across periodic boundaries (not
  //**everywhere)
  if (problemDomain().isPeriodic())
    {
      m_J.exchange();
    }

  // In the circumstance where finer levels are being completely removed, they
  // are now marked as undefined.
  if (undefFinerLevels)
    {
      m_finerLevelGridMetrics->undefineMetrics();
    }

  // Only keep the metrics changed in the base level because we've accounted for
  // all the finer levels
  if (!isBase)
    {
      m_metricsChanged = false;
    }
}

/*--------------------------------------------------------------------*/
//  Regrid operations -- performs a coarse to fine interpolation
/** This routine interpolates from a coarse grid to a fine grid
 *  using constrained least-squares with considerations for the
 *  metrics.  The interpolator, 'FourthOrderMappedFineInterp' can
 *  alternatively be used externally with 'm_J'.
 *  \param[out] a_JU    Interpolated \<JU\> for this level
 *  \param[in]  a_CrU   \<U\> (solution) from the next coarser level
 *  \param[in]  a_CrJU  \<JU\> from the next coarser level
 *  \param[in]  a_vectorIntv
 *                      An interval in components of U consisting of
 *                      vector data
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::regrid(LevelData<FArrayBox>&       a_JU,
                         const LevelData<FArrayBox>& a_CrU,
                         const LevelData<FArrayBox>& a_CrJU,
                         const Interval&             a_vectorIntv)
{
  CH_TIME("LevelGridMetrics::regrid");
  if (verbosity() >= 3)
    {
      pout() << "LevelGridMetrics::regrid " << level() << std::endl;
    }
  if (hasCoarser())
    {
      IntVect refRatioVec(D_DECL6(
                            m_coarserLevelGridMetrics->refRatio(),
                            m_coarserLevelGridMetrics->refRatio(),
                            m_coarserLevelGridMetrics->refRatio(),
                            m_coarserLevelGridMetrics->refRatio(),
                            m_coarserLevelGridMetrics->refRatio(),
                            m_coarserLevelGridMetrics->refRatio()));
      // The interpolator lives on the coarser level
      FourthOrderMappedFineInterp* interpolator =
        m_coarserLevelGridMetrics->m_interpolator;
      interpolator->defineStencils(
        refRatioVec,
        m_coarserLevelGridMetrics->dxVect());
      const MultiBlockLevelGeom *const mbgeoPtr =
        (isMultiblock()) ? &(m_coarserLevelGridMetrics->m_mbgeo) : 0;
      interpolator->defineGrids(
        getBoxes(),
        m_coarserLevelGridMetrics->getBoxes(),
        m_coarserLevelGridMetrics->problemDomain(),
        m_ghostVect,
        true,
        mbgeoPtr);
      interpolator->defineCrFnJ(
        m_coarserLevelGridMetrics->m_J);
      interpolator->defineCFInterface(m_coarserLevelGridMetrics);
      interpolator->interpToFine(
        a_JU,
        m_J,
        a_CrU,
        a_CrJU,
        m_coarserLevelGridMetrics,
        a_vectorIntv);
    }
}

/*--------------------------------------------------------------------*/
//  Post regrid operations -- currently unused.
/*--------------------------------------------------------------------*/

// void
// LevelGridMetrics::postRegrid(int a_base_level)
// {
// }

/*--------------------------------------------------------------------*/
//  Fill invalid fine ghost cells (with preset CrFnLevU)
/** This routine interpolates from a coarse grid to a fine grid
 *  using constrained least-squares with considerations for the
 *  metrics.  The coarsened-fine data \<U\> for the interpolator must
 *  have already been set, presumably using getCrFnLevU().  The
 *  interpolator, 'FourthOrderMappedFineInterp' can alternatively be
 *  used externally with 'm_J'.
 *  \param[out] a_U     Interpolated \<U\> on invalid ghost cells for
 *                      this level
 *  \param[out] a_JU    Interpolated \<JU\> on invalid ghost cells for
 *                      this level
 *  \param[in] a_time   Time at which we would like the coarse data
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::fillFineGhostCells(LevelData<FArrayBox>& a_U,
                                     LevelData<FArrayBox>& a_JU,
                                     Real a_time)
{
  CH_TIME("LevelGridMetrics::fillFineGhostCells");
  if (verbosity() >= 4)
    {
      pout() << "LevelGridMetrics::fillFineGhostCells " << level() << std::endl;
    }
  CH_assert(hasCoarser());
  // The interpolator lives on the coarser level
  FourthOrderMappedFineInterp* interpolator =
    m_coarserLevelGridMetrics->m_interpolator;
  interpolator->fillGhosts(
    a_U,
    a_JU,
    m_J,
    m_coarserLevelGridMetrics);
}

/*--------------------------------------------------------------------*/
//  Computes <U> in valid cells
/** Compute the new \<U\>.  Use the product formula and let
 *  Grad(U) = Grad(\<JU\>/\<J\>).
 *  \param[out] a_U     \<U\> on all valid cells
 *  \param[in]  a_JU    \<JU\> on all valid cells plus 1 layer of
 *                      valid ghost cells except at boundaries of
 *                      either the problem domain or block.  At those
 *                      locations, one-sided derivatives are used.
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::computeValidU(LevelData<FArrayBox>&       a_U,
                                const LevelData<FArrayBox>& a_JU)
{
  CH_TIME("LevelGridMetrics::computeValidU");
  for (DataIterator dit = getDataIterator(); dit.ok(); ++dit)
    {
      const Box interiorBox = getBoxes()[dit];
      cellFGToCellF(a_U[dit],
                    a_JU[dit],
                    m_J[dit],
                    interiorBox,
                    blockDomain(interiorBox));
    }
}

/*--------------------------------------------------------------------*/
//  Compute the minimum grid buffer size for a fourth-order
//  interpolation
/** \param[in] a_refRatios
 *                      Refinement ratios
 *  \param[in] a_maxLevel
 *                      Index of top level (base level == 0)
 *  \param[in] a_numGhost
 *                      Number of ghost cells (fine ghost cells that
 *                      must be filled at coarse-fine boundaries
 *  \return             Minimum grid buffer size
 *//*-----------------------------------------------------------------*/

int
LevelGridMetrics::bufferSize4thO(const std::vector<int>& a_refRatio,
                                 const int               a_maxLevel,
                                 const int               a_numGhost)
{
  if (a_maxLevel == 0)
    {
      return 1;
    }
  CH_assert(a_maxLevel <= a_refRatio.size());
  // Find the minimum ref ratio
  int minNRef = std::numeric_limits<int>::max();
  for (int i = 0; i != a_maxLevel; ++i)
    {
      minNRef = std::min(minNRef, a_refRatio[i]);
    }
  CH_assert(minNRef > 0);
  // Compute the buffer size
  const int buffer =
    // Coarse cells underlying the fine cells
    (a_numGhost + minNRef - 1)/minNRef
    // Addition for stencils
    + 2
    // Addition for computing <U> from <JU>
    + 1;
  return buffer;
}

/*--------------------------------------------------------------------*/
//  Exchange \<U\> across multiblock boundaries
/** \param[in]  a_U     \<U\> known in valid cells
 *  \param[out] a_U     \<U\> updated in extra-block ghost cells
 *  \param[in]  a_vectorIntv
 *                      An interval in components of U consisting of
 *                      vector data
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::multiblockExchangeU(LevelData<FArrayBox>& a_U,
                                      const Interval&       a_vectorIntv) const
{
  if (isMultiblock())
    {
      m_mbex.interpGhostsAllWithVector(a_U, a_vectorIntv);
    }
}

/*--------------------------------------------------------------------*/
//  Supports a call to TimeInterpolatorRK4::intermediate for multi-
//  block grids
/** This is here to support a hack so that we can convert \<JU\> to
 *  \<U\> and do a multiblock copyTo.  For single-block grids, the
 *  time interpolator does a copyTo of \<JU\>  which needs no
 *  alteration
 *  \param[in]  a_timeIterpolator
 *                      RK4 time interpolator
 *  \param[in]  a_timeInterpCoeff
 *                      Time interpolation coefficient in range [0:1]
 *  \param[in]  a_stage Which RK4 stage:  0, 1, 2, 3
 *  \param[in]  a_vectorIntv
 *                      An interval in components of U consisting of
 *                      vector data
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::timeIntermediate(
  const TimeInterpolatorRK4& a_timeInterpolator,
  const Real&                a_timeInterpCoeff,
  const int&                 a_stage,
  const Interval&            a_vectorIntv)
{
  CH_TIME("LevelGridMetrics::timeIntermediate");
  CH_assert(hasCoarser());
  invalidateCr2ThisInterpolatorCrFnLevData();
  if (isMultiblock())
    {
      const DisjointBoxLayout crGrid = m_coarserLevelGridMetrics->getBoxes();
      // Expensive allocation --- but I hope this whole idea is temporary
      // We only need <U> on some properly nested subset of the coarser-level.
      // We obtain <JU> on all coarse cells, exchange to fill ghosts, and
      // compute <U> everywhere but there will be garbage in the layer of cells
      // next to the interferce with the even coarser mesh.
      // CrFnLevJU has been sized to the coarse grid by the interpolator
      LevelData<FArrayBox>& crLevJU = presetCr2ThisInterpolatorCrFnLevJU();
      LevelData<FArrayBox> crLevU(crGrid,
                                  m_interpolator->nComp(),
                                  IntVect::Zero);
      Interval intv(0, m_interpolator->nComp()-1);
      a_timeInterpolator.intermediate(
        crLevJU,
        a_timeInterpCoeff,
        a_stage,
        intv);
      crLevJU.exchange();  // To set the invalid ghosts

      const LevelData<FArrayBox>& crLevJ = m_coarserLevelGridMetrics->m_J;
      for (DataIterator dit = crGrid.dataIterator(); dit.ok(); ++dit)
        {
          const Box interiorBox = crGrid[dit];
          // Compute Cr <U> in the domain of the block.  This calculation will
          // give erroneous results in cells adjacent to a coarser mesh (but
          // proper nesting means we don't need <U> from there).
          cellFGToCellF(crLevU[dit],
                        crLevJU[dit],
                        crLevJ[dit],
                        interiorBox,
                        m_coarserLevelGridMetrics->blockDomain(interiorBox));
          //**CHECKME   true); - FourthOrderMappedFineInterp::fillGhosts has
          //                     this set to true but why??  Metrics should be
          //                     available everywhere to allow for centered
          //                     gradients of <J>.
        }

//--Here we finally get <U> over to a coarsened representation of the fine
//--grid

      LevelData<FArrayBox>& crFnLevU = presetCr2ThisInterpolatorCrFnLevU();
      // Normal copy
      crLevU.copyTo(crFnLevU,
                    m_coarserLevelGridMetrics->m_interpolator->m_copier);
      // Across multiblock boundaries
      MultiBlockLevelCopierAverage& mbcp =
        m_coarserLevelGridMetrics->m_interpolator->m_mbcp;
      //**FIXME Grumble... vector data in between scalars.  Should we copy it or
      //**      change the indices in LGintergrator.H?
      const Interval Uintv = crFnLevU.interval();
      if (a_vectorIntv.size() == 0)
        { // In this case, all components are scalar, none vector.
          mbcp.copyTo(crLevU, crFnLevU);
        }
      else
        {
          // Ugly ... ow ow ow!
          if (Uintv.begin() < a_vectorIntv.begin())
            { // scalar component indices before a_vectorIntv
              Interval scalarIntv(Uintv.begin(), a_vectorIntv.begin()-1);
              mbcp.copyTo(crLevU, crFnLevU, scalarIntv);
            }
          mbcp.copyToVector(crLevU, crFnLevU, a_vectorIntv);
          if (a_vectorIntv.end() < Uintv.end())
            { // scalar component indices after a_vectorIntv
              Interval scalarIntv(a_vectorIntv.end()+1, Uintv.end());
              mbcp.copyTo(crLevU, crFnLevU, scalarIntv);
            }
        }
    }
  else  // Sanity
    {
      Interval intv(0, m_interpolator->nComp()-1);
      a_timeInterpolator.intermediate(
        presetCr2ThisInterpolatorCrFnLevJU(),
        a_timeInterpCoeff,
        a_stage,
        intv);
    }
}

/*--------------------------------------------------------------------*/
//  Return the I/O verbosity
/*--------------------------------------------------------------------*/

int
LevelGridMetrics::verbosity()
{
  return AMRLevel::verbosity();
}

/*--------------------------------------------------------------------*/
//  Completely define the metrics for this level
/*--------------------------------------------------------------------*/

void
LevelGridMetrics::defineMetrics()
{
  CH_TIME("LevelGridMetrics::defineMetrics");
  if (verbosity() >= 3)
    {
      pout() << "LevelGridMetrics::defineMetrics " << level() << std::endl;
    }

  if (m_useScriptN)
    {

//--N is calculated using script N

      // Script N
      LevelData<CodimBox<FArrayBox> > scrN(m_gridsNP,
                               SpaceDim,
                               getNReqGhostVect(),
                               CodimBoxFactory<FArrayBox> (2));
      getScriptN(scrN);
      averageScriptN(false, scrN);

      // N
      const int numNTerms = SpaceDim*SpaceDim;
      m_N.define(getBoxes(), numNTerms, getNReqGhostVect());
      faceNormalN(m_N, scrN);
      // Option for transverse terms
      // If averaging, transverse components are calculated from normal
      // components and you should exchange first (note that if averaging,
      // transverse components are only computed on the faces of valid cells).
      // Otherwise exchange after.
      if (m_transverseNOpType != TransverseNOpAverage)
        {
          faceTransverseN(m_N);
        }
      //**FIXME This should only exchange across periodic boundaries (not
      //**everywhere)
      if (problemDomain().isPeriodic())
        {
          m_N.exchange();
        }
      if (m_transverseNOpType == TransverseNOpAverage)
        {
          faceTransverseN(m_N);
        }
    }
  else
    {

//--N is calculated directly (not freestream preserving)

      // N
      const int numNTerms = SpaceDim*SpaceDim;
      m_N.define(getBoxes(), numNTerms, getNReqGhostVect());
      getN(m_N);
      averageN(false, m_N);
    }

  // NtX
  LevelData<FluxBox> NtXLev(m_gridsNP, 1, getJReqGhostVect());
  getVolFlux(NtXLev, m_N); // NtXLev is contiguous in cubed-sphere.
  averageVolFlux(false, NtXLev, AgeNew);
  // petermc, 21 Nov 2011: Make NtXLev NON-contiguous in cubed-sphere.
  m_coordSys->separateVolFlux(NtXLev);

  // J
  m_J.define(getBoxes(), 1, getJReqGhostVect());
  // Update m_J
  setAvgJfromNtX(NtXLev);

  //**FIXME This should only exchange across periodic boundaries (not
  //**everywhere)
  if (problemDomain().isPeriodic())
    {
      m_J.exchange();
    }

  m_metricsDefined = true;
  if (hasFiner())
    {
      m_finerLevelGridMetrics->m_metricsChanged = false;
    }
}

/*--------------------------------------------------------------------*/
//  Label this and all finer levels as having undefined metrics
/** This only happens when fine levels are completely removed from an
 *  AMR simulation
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::undefineMetrics()
{
  CH_TIME("LevelGridMetrics::undefineMetrics");
  if (hasFiner())
    {
      m_finerLevelGridMetrics->undefineMetrics();
    }

  if (verbosity() >= 3)
    {
      pout() << "LevelGridMetrics::undefineMetrics " << level() << std::endl;
    }

  // Assign no boxes to the grids being undefined
  Vector<Box> noBox(0);
  Vector<int> noProc(0);
  DisjointBoxLayout zeroBoxes(noBox, noProc, problemDomain());
  zeroBoxes.close();
  // Grids and multiblock structures are defined with zero boxes
  defineGrids(zeroBoxes);
  m_metricsDefined = false;
  // Invalidate all CrFn data in the interpolator
  invalidateThis2FnInterpolatorCrFnLevData();
}

/*--------------------------------------------------------------------*/
//  Compute \f$N\f$ on the hyperfaces from \f$\mathcal(N)^s\f$
/** These are the normal components of \f$\mathbf(N)\f$, i.e.,
 *  a specification of a column \f$d\f$ of \f$\mathbf(N)^s_d\f$
 *  for face \f$d\f$.  The are obtained from \f$\mathcal(N)^s\f$ on
 *  the hyperedges.
 *  \param[out] a_N     \f$N\f$
 *  \param[in]  a_scrN  \f$\mathcal(N)\f$ for building \f$N\f$
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::faceNormalN(LevelData<FluxBox>&        a_N,
                              const LevelData<CodimBox<FArrayBox> >& a_scrN)
{
  CH_TIME("LevelGridMetrics::faceNormalN");
  CH_assert(m_useScriptN);  // Should only be called if using script N
  DataIterator dit = a_N.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box box = a_N.disjointBoxLayout()[dit()];
      const NewFourthOrderCoordSys *const coordSys = getCoordSys(box);
      const CodimBox<FArrayBox> & scrN = a_scrN[dit];
      FluxBox& N = a_N[dit];
      N.setVal(0.);
      // Loop over all hyperedge orientations
      const int nOr = scrN.getNumOrient();
      for (int iOr = 0; iOr != nOr; ++iOr)
        {
          // A codimension 2 object is defined by 2 directions.  Get them.
          int dir[2];
          scrN.getDirections(iOr, dir);
          const FArrayBox& hyperEdge = scrN.getSequential(iOr);
          // Loop over each hyperface adjacent to the hyperedge
          for (int iFace = 0; iFace < 2; ++iFace)
            {
              int faceDir = dir[iFace];   // d
              int edgeDir = dir[!iFace];  // d'
              // Begin of column of N (N^s) for this face (normal components)
              const int iNBeg = coordSys->getNcomponent(0, faceDir);
              // Script N is stored for d=dir[0], d'=dir[1].  If the indices are
              // swapped, change the sign.
              const int sign = 1 - 2*iFace;
              const Real mult =  // Includes sign and 1/h^(D-1)
                sign/std::pow(m_dxVect[faceDir], SpaceDim - 1);
              FArrayBox& hyperFace = N[faceDir];
              Box faceBox = hyperFace.box();       // Box orientated to the face
              // Loop over each hyperface in faceDir in the box
              BoxIterator bit(faceBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect ivFace = bit();
                  IntVect ivLo = ivFace;
                  IntVect ivHi = ivFace + BASISV(edgeDir);
                  // Loop over rows 's' and subtract/add low/high scrN
                  for (int s = 0; s < SpaceDim; ++s)
                    {
                      hyperFace(ivFace, iNBeg + s) += mult *
                        ( hyperEdge(ivHi, s) - hyperEdge(ivLo, s) );
                    }
                }
            }
        }
    }
  m_metricsChanged = true;
}

/*--------------------------------------------------------------------*/
//  Compute \f$N\f$ on the hyperfaces from \f$\mathcal{N}^s\f$
/** All components of \f$\mathbf{N}\f$ are computed directly.  This is
 *  not freestream-preserving.
 *  \param[out] a_N     \f$N\f$
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::getN(LevelData<FluxBox>& a_N)
{
  CH_TIME("LevelGridMetrics::getN");
  CH_assert(!m_useScriptN);  // Should only be called if not using script N
  DataIterator dit = a_N.dataIterator();
  const DisjointBoxLayout& layout = a_N.disjointBoxLayout();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box baseBox = layout[dit];
      const NewFourthOrderCoordSys *const coordSys = getCoordSys(baseBox);
      FluxBox& N = a_N[dit];
      const Box bx = N.box();
      coordSys->getN(N, bx);
      for (int idir=0; idir<SpaceDim; idir++)
        N[idir] /= m_dxVect[idir];
    }
  m_metricsChanged = true;
}

/*--------------------------------------------------------------------*/
//  Add transverse components to \f$N\f$
/** These are the transverse components of \f$\mathbf{N}\f$, i.e.,
 *  a specification of the columns \f$d'\ne d\f$ of
 *  \f$\mathbf(N)^s_(d'}\f$ for face \f$d\f$.
 *  \param[out] a_N     \f$N\f$
 *  \note
 *  <ul>
 *    <li> Transverse terms can be computed analytically or by
 *         averaging normal components from faces in other directions
 *    <li> If determined analytically, they should be available on the
 *         faces of all cells
 *    <li> If determined by averages, the transverse terms are only
 *         determined on the faces of grow(N.box(), -1)
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::faceTransverseN(LevelData<FluxBox>& a_N)
{
  CH_TIME("LevelGridMetrics::faceTransverseN");
  switch(m_transverseNOpType)
    {
    case TransverseNOpNone:
      break;
    case TransverseNOpAnalytic:
      for (DataIterator dit = a_N.dataIterator(); dit.ok(); ++dit)
        {
          Box box = a_N.disjointBoxLayout()[dit()];
          const NewFourthOrderCoordSys *const coordSys = getCoordSys(box);
          FluxBox& N = a_N[dit];
          // Loop over faces (d)
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              // factor = 1/faceArea
              const Real factor = m_dxVect[faceDir]/m_dxVect.product();
              FArrayBox& hyperFace = N[faceDir];
              // Loop over component directions d'
              for (int dp = 0; dp != SpaceDim; ++dp)
                {
                  if (dp != faceDir)
                    {
                      coordSys->computeTransverseFaceMetric(hyperFace,
                                                            hyperFace.box(),
                                                            faceDir,
                                                            dp);
                      // This result is multiplied by the face area -- get the
                      // actual <N> by dividing out the face area
                      const int startTransverseComp =
                        coordSys->getNcomponent(0, dp);
                      hyperFace.mult(factor, startTransverseComp, SpaceDim);
                    }
                }
            }
        }
      break;
    case TransverseNOpAverage:
      for (DataIterator dit = a_N.dataIterator(); dit.ok(); ++dit)
        {
          FluxBox& N = a_N[dit];
          const Box box = grow(N.box(), -1);
          // Loop over faces (d)
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              FArrayBox& hyperFace = N[faceDir];
              // Loop over component directions d'
              for (int dp = 0; dp != SpaceDim; ++dp)
                {
                  if (dp != faceDir)
                    {
                      FArrayBox& transverseHyperFace = N[dp];
                      averageTransverseN(hyperFace,
                                         transverseHyperFace,
                                         box,
                                         faceDir,
                                         dp);
                    }
                }
            }
        }
      break;
    }
}

/*--------------------------------------------------------------------*/
//  Integrate \f$\mathcal{N}^s\f$ on each hyperedge
/** \f$\mathcal{N}^s\f$ is computed on both physical and ghost edges
 *  \param[out] a_scrN  Storage for \f$\mathcal(N)^s\f$.  Must have
 *                      SpaceDim components.
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::getScriptN(LevelData<CodimBox<FArrayBox> >& a_scrN)
{
  CH_TIME("LevelGridMetrics::getScriptN");
  const int nOr = CodimBox<FArrayBox> ::getNumOrient(2);
  DataIterator dit = a_scrN.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box box = a_scrN.disjointBoxLayout()[dit()];
      const NewFourthOrderCoordSys *const coordSys = getCoordSys(box);
      CodimBox<FArrayBox>  &scrN = a_scrN[dit];

// Compute scrN on each codimensional 2 geometry in the block including in the
// ghost elements.

      for (int iOr = 0; iOr != nOr; ++iOr)
        {
          // A codimension 2 object is defined by 2 directions.  Get them.
          int dir[2];
          scrN.getDirections(iOr, dir);
          FArrayBox& hyperEdge = scrN.getSequential(iOr);
          Box oBox = hyperEdge.box();    // An orientated box
          coordSys->integrateScriptN(hyperEdge, dir[0], dir[1], oBox);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Average \f$\mathcal{N}^s\f$.
/** This is performed in two steps.  First, we overwrite using a
 *  saved average from a finer level.  Then we store the average for
 *  a coarser level.
 *  \param[in]  a_isBase
 *                      Set to true if this is a base level during a
 *                      regrid in which case don't average to the
 *                      next coarser level.
 *  \param[in]  a_scrN  Computed \f$\mathcal(N)^s\f$ on this level.
 *  \param[out] a_scrN  Possibly overwritten with averages from a
 *                      finer level.
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::averageScriptN(const bool           a_isBase,
                                 LevelData<CodimBox<FArrayBox> >& a_scrN)
{
  CH_TIME("LevelGridMetrics::averageScriptN");
  CH_assert(m_useScriptN);  // Should only be called if using script N
  if (m_scrNAverager)  // Only if has finer
    {
      m_scrNAverager->copyTo(a_scrN);
      delete m_scrNAverager;
      m_scrNAverager = NULL;
    }
  if (!a_isBase && hasCoarser())
    {
      m_coarserLevelGridMetrics->m_scrNAverager =
        new CoarseAverageCodim(
          a_scrN.getBoxes(),
          2,
          SpaceDim,
          m_coarserLevelGridMetrics->refRatio());
      m_coarserLevelGridMetrics->m_scrNAverager->sum(a_scrN);
    }
}

/*--------------------------------------------------------------------*/
//  Average \f$N\f$ directly if \f$\mathcal{N}^s\f$ doesn't exist.
/** This is performed in two steps.  First, we overwrite using a
 *  saved average from a finer level.  Then we store the average for
 *  a coarser level.
 *  \param[in]  a_isBase
 *                      Set to true if this is a base level during a
 *                      regrid in which case don't average to the
 *                      next coarser level.
 *  \param[in]  a_N     Computed \f$N\f$ on this level.
 *  \param[out] a_N     Possibly overwritten with averages from a
 *                      finer level.
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::averageN(const bool          a_isBase,
                           LevelData<FluxBox>& a_N)
{
  CH_TIME("LevelGridMetrics::averageN");
  CH_assert(!m_useScriptN);  // Should only be called if not using script N
  if (m_NAverager)  // Only if has finer
    {
      m_NAverager->copyTo(a_N);
      delete m_NAverager;
      m_NAverager = NULL;

      // MB face register
      if (isMultiblock())
        {
          // Define the coarse side
          m_Nmbreg->defineCoarse(a_N.getBoxes());
          Interval intv = a_N.interval();
          // Update totals
          // 'false' means that N is not directional (although it has vector
          // components, the direction in a face is meaningless)
          m_Nmbreg->reface(a_N, intv, intv, -1., false);
          delete m_Nmbreg;
          m_Nmbreg = NULL;
        }
    }

  if (!a_isBase && hasCoarser())
    {
      m_coarserLevelGridMetrics->m_NAverager =
        new CoarseAverageFace(
          a_N.getBoxes(),
          SpaceDim*SpaceDim,
          m_coarserLevelGridMetrics->refRatio());
      m_coarserLevelGridMetrics->m_NAverager->average(a_N);

      // MB face register
      if (isMultiblock())
        {
          m_coarserLevelGridMetrics->m_Nmbreg = new MultiBlockFaceRegister(
            &(m_coarserLevelGridMetrics->getCoordSys()));
          MultiBlockFaceRegister *const NmbregCr =
            m_coarserLevelGridMetrics->m_Nmbreg;
          NmbregCr->defineFine(a_N.getBoxes(),
                               problemDomain(),
                               m_coarserLevelGridMetrics->refRatio(),
                               a_N.nComp());
          NmbregCr->setToZero();
          Interval intv = a_N.interval();
          // Increment fine values
          for (DataIterator dit = a_N.dataIterator(); dit.ok(); ++dit)
            {
              FluxBox &NFlBx = a_N[dit];
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  // Separate out sides to avoid a change in sign (we aren't
                  // actually doing a divergence)
                  NmbregCr->incrementFine(NFlBx[dir],
                                          -1.,
                                          dit(),
                                          intv,
                                          intv,
                                          dir,
                                          Side::Lo);
                  NmbregCr->incrementFine(NFlBx[dir],
                                          1.,
                                          dit(),
                                          intv,
                                          intv,
                                          dir,
                                          Side::Hi);
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the volume flux
/** \param[out] a_NtX   Volume flux
 *  \param[in]  a_N     Metrics on the faces
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::getVolFlux(LevelData<FluxBox>&        a_NtX,
                             const LevelData<FluxBox>&  a_N)
{
  CH_TIME("LevelGridMetrics::volFlux");
  DataIterator dit = a_N.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box box = a_N.disjointBoxLayout()[dit()];
      const NewFourthOrderCoordSys *const coordSys = getCoordSys(box);
      Box cellBox = resizeWithJReqGhosts(a_N.getBoxes()[dit]);
      coordSys->volFlux(a_NtX[dit], a_N[dit], cellBox);
    }
}


/*--------------------------------------------------------------------*/
//  Compute the divergence of the volume flux to set <J> in m_J
/** \param[in] a_NtX   Volume flux
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::setAvgJfromNtX(LevelData<FluxBox>& a_NtX)
{
  CH_TIME("LevelGridMetrics::setAvgJfromNtX");
  DataIterator dit = m_J.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box box = getBoxes()[dit];
      const NewFourthOrderCoordSys *const coordSys = getCoordSys(box);
      // Get the box we need to compute J on
      Box cellBox = m_J[dit].box();
      // In cubed-sphere, NtX is NON-contiguous.
      coordSys->getAvgJ(m_J[dit], a_NtX[dit], cellBox);
    }
}

/*--------------------------------------------------------------------*/
//  Average the volume flux \f$N^T X\f$.
/** This is performed in two steps.  First, we overwrite using a
 *  saved average from a finer level.  Then we store the average for
 *  a coarser level.
 *  \param[in]  a_isBase
 *                      Set to true if this is a base level during a
 *                      regrid in which case don't average to the
 *                      next coarser level.
 *  \param[in]  a_NtX   Computed \f$N^T X\f$ on this level.
 *  \param[out] a_NtX   Possibly overwritten with averages from a
 *                      finer level.
 *  \param[in]  a_age   During a regrid, averages exists for both the
 *                      new and old metrics.  This specifies which one
 *                      to average.
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::averageVolFlux(const bool          a_isBase,
                                 LevelData<FluxBox>& a_NtX,
                                 const AgeType       a_age)
{
  CH_TIME("LevelGridMetrics::averageVolFlux");
  if (m_NtXAverager[a_age])  // Only if has finer
    {
      m_NtXAverager[a_age]->copyTo(a_NtX);
      delete m_NtXAverager[a_age];
      m_NtXAverager[a_age] = NULL;

      // MB face register
      if (isMultiblock())
        {
          // Define the coarse side
          m_NtXmbreg[a_age]->defineCoarse(a_NtX.getBoxes());
          Interval intv = a_NtX.interval();
          // Update totals
          //**FIXME The false argument should be 'true' in general since the
          //**       volume flux is directional.  But for the
          //**       CubedSphere2DPanelCS coordsys, NtX always matches sign
          //**       across boundaries, even if directions change sign.
          m_NtXmbreg[a_age]->reface(a_NtX, intv, intv, -1., false);
          delete m_NtXmbreg[a_age];
          m_NtXmbreg[a_age] = NULL;
        }
    }

  if (!a_isBase && hasCoarser())
    {
      m_coarserLevelGridMetrics->m_NtXAverager[a_age] =
        new CoarseAverageFace(
          a_NtX.getBoxes(),
          1,
          m_coarserLevelGridMetrics->refRatio());
      m_coarserLevelGridMetrics->m_NtXAverager[a_age]->average(a_NtX);

      // MB face register
      if (isMultiblock())
        {
          m_coarserLevelGridMetrics->m_NtXmbreg[a_age] =
            new MultiBlockFaceRegister(
              &(m_coarserLevelGridMetrics->getCoordSys()));
          MultiBlockFaceRegister *const NtXmbregCr =
            m_coarserLevelGridMetrics->m_NtXmbreg[a_age];
          NtXmbregCr->defineFine(a_NtX.getBoxes(),
                                 problemDomain(),
                                 m_coarserLevelGridMetrics->refRatio(),
                                 a_NtX.nComp());
          NtXmbregCr->setToZero();
          Interval intv = a_NtX.interval();
          // Increment fine values
          for (DataIterator dit = a_NtX.dataIterator(); dit.ok(); ++dit)
            {
              FluxBox &NtXFlBx = a_NtX[dit];
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  // Separate out sides to avoid a change in sign (we aren't
                  // actually doing a divergence)
                  NtXmbregCr->incrementFine(NtXFlBx[dir],
                                            -1.,
                                            dit(),
                                            intv,
                                            intv,
                                            dir,
                                            Side::Lo);
                  NtXmbregCr->incrementFine(NtXFlBx[dir],
                                            1.,
                                            dit(),
                                            intv,
                                            intv,
                                            dir,
                                            Side::Hi);
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute \f$N^T X\f$ on all finer levels and average down
/** The base of the recursion must always have a coarser level
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::volFluxCascade()
{
  CH_TIME("LevelGridMetrics::volFluxCascade");
  // Walk up the hierarchy
  if (hasFiner()) 
    m_finerLevelGridMetrics->volFluxCascade();

  // Compute the volume flux
  LevelData<FluxBox> NtXLev(m_gridsNP, 1, getJReqGhostVect());
  getVolFlux(NtXLev, m_N);
  // Average down from finer and save averages to coarser.  A walk through this
  // cascade is to update the volume flux on the old meshes.
  averageVolFlux(false, NtXLev, AgeOld);
}

/*--------------------------------------------------------------------*/
//  Average the snapback solution flux, \f$N^T F\f$.
/** This is performed in two steps.  First, we overwrite using a
 *  saved average from a finer level.  Then we store the average for
 *  a coarser level.
 *  \param[in]  a_isBase
 *                      Set to true if this is a base level during a
 *                      regrid in which case don't average to the
 *                      next coarser level.
 *  \param[in]  a_NtFDiff
 *                      Computed \f$N^T F\f$ on this level.
 *  \param[out] a_NtFDiff
 *                      Possibly overwritten with averages from a
 *                      finer level.
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::averageSnapbackSolFlux(const bool          a_isBase,
                                         LevelData<FluxBox>& a_NtFDiff)
{
  CH_TIME("LevelGridMetrics::averageSnapbackSolFlux");
  if (m_NtFDiffAverager)  // Only if has finer
    {
      m_NtFDiffAverager->copyTo(a_NtFDiff, m_finerLevelGridMetrics->getBoxes());
      delete m_NtFDiffAverager;
      m_NtFDiffAverager = NULL;
    }
  if (!a_isBase && hasCoarser())
    {
      m_coarserLevelGridMetrics->m_NtFDiffAverager =
        new IntersectionAverageFace(
          a_NtFDiff,
          m_coarserLevelGridMetrics->refRatio());
    }
}

/*--------------------------------------------------------------------*/
//  Store the DisjointBoxLayout and define multiblock structures
/** \param[in]  a_grids Original grid
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::defineGrids(const DisjointBoxLayout& a_grids)
{
  CH_TIME("LevelGridMetrics::defineGrids");
  m_grids = a_grids;
  if (m_coordSys->numBlocks() > 1)
    {
      m_mbgeo.define(m_coordSys, m_grids, m_ghostVect[0], m_spaceOrder);
      m_mbex.define(&m_mbgeo, m_ghostVect[0], m_spaceOrder);
      if (m_haveMultiBlockVectorData)
        {
          m_mbex.defineVector();
        }
    }
}

/*--------------------------------------------------------------------*/
//  Create a non-periodic DisjointBoxLayout
/** \param[in]  a_grids Original grid
 *//*-----------------------------------------------------------------*/

void
LevelGridMetrics::defineGridsNP(const DisjointBoxLayout& a_grids)
{
  CH_TIME("LevelGridMetrics::defineGridsNP");
  if (problemDomain().isPeriodic())
    {
      const bool periodicDir[SpaceDim] =
        {
          D_DECL6(false, false, false, false, false, false)
        };
      const ProblemDomain probDomNP(problemDomain().domainBox(),
                                    periodicDir);
      m_gridsNP.define(a_grids, probDomNP);
    }
  else
    {
      m_gridsNP = a_grids;
    }
}

#include "NamespaceFooter.H"
