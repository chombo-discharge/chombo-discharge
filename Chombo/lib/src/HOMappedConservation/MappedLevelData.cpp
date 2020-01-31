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
 * \file MappedLevelData.cpp
 *
 * \brief Non-inline definitions for classes in MappedLevelData.H
 *
 *//*+*************************************************************************/

#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "MultiBlockLevelExchange.H"
#include "NamespaceHeader.H"


/*******************************************************************************
 *
 * Class MappedLevelData: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/

MappedLevelData::MappedLevelData()
  :
  m_indexJUOld(0),
  m_indexJUNew(1),
  m_allocU(-1),
  m_allocJUNew(-1),
  m_allocJUOld(-1),
  m_haveValidCrFnU(false),
  m_defined(false)
{
}

/*--------------------------------------------------------------------*/

MappedLevelData::~MappedLevelData()
{
  deallocate();
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
/**
 *   \note
 *   <ul>
 *     <li> Call again to redefine
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
MappedLevelData::define(LevelGridMetrics *const           a_LGMP,
                        LevelData<FArrayBox> *const       a_ULevP,
                        LevelData<FArrayBox> *const       a_JUNewLevP,
                        LevelData<FArrayBox> *const       a_JUOldLevP,
                        const LevelData<FArrayBox> *const a_CrULevP,
                        const int                         a_numComp,
                        const Interval&                   a_velocityIntv,
                        const int                         a_numGhostU,
                        const int                         a_numGhostJU,
                        const bool                        a_trimUExchanger,
                        const bool                        a_trimJUExchanger,
                        const bool                        a_useJUExchangerForU)
{
  m_numComp = a_numComp;
  m_velocityIntv = a_velocityIntv;
  m_numGhostU = a_numGhostU;
  m_numGhostJU = a_numGhostJU;
  m_trimUExchanger = a_trimUExchanger;
  m_trimJUExchanger = a_trimJUExchanger;
  m_useJUExchangerForU = a_useJUExchangerForU;

  deallocate();

  m_LGMP = a_LGMP;

  if (a_ULevP == 0)
    {
      m_ULevP = new LevelData<FArrayBox>;
      m_allocU = 0;
    }
  else
    {
      m_ULevP = a_ULevP;
    }

  if (a_JUNewLevP == 0)
    {
      m_JULevP[m_indexJUNew] = new LevelData<FArrayBox>;
      m_allocJUNew = m_indexJUNew;
    }
  else
    {
      m_JULevP[m_indexJUNew] = a_JUNewLevP;
    }

  if (a_JUOldLevP == 0)
    {
      m_JULevP[m_indexJUOld] = new LevelData<FArrayBox>;
      m_allocJUOld = m_indexJUOld;
    }
  else
    {
      m_JULevP[m_indexJUOld] = a_JUOldLevP;
    }

  m_CrULevP = a_CrULevP;

  m_defined = true;
  reshapeU();
  reshapeJUNew();
  reshapeJUOld();

  invalidate();
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::reshapeU()
{
  CH_assert(m_defined);
  const IntVect ghostVectU = m_numGhostU*IntVect::Unit;
  rawU().define(m_LGMP->getBoxes(), m_numComp, ghostVectU);
  m_UExchangeCopier.exchangeDefine(m_LGMP->getBoxes(), ghostVectU);
  if (m_trimUExchanger)
    {
      m_UExchangeCopier.trimEdges(m_LGMP->getBoxes(), ghostVectU);
    }
  m_haveValidCrFnU = false;
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::reshapeJUNew()
{
  CH_assert(m_defined);
  const IntVect ghostVectJU = m_numGhostJU*IntVect::Unit;
  rawJU().define(m_LGMP->getBoxes(), m_numComp, ghostVectJU);
  m_JUExchangeCopier.exchangeDefine(m_LGMP->getBoxes(), ghostVectJU);
  if (m_trimJUExchanger)
    {
      m_JUExchangeCopier.trimEdges(m_LGMP->getBoxes(), ghostVectJU);
    }
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::reshapeJUOld()
{
  CH_assert(m_defined);
  getJUOld().define(m_LGMP->getBoxes(), m_numComp, m_numGhostJU*IntVect::Unit);
}

/*--------------------------------------------------------------------*/
//
/** The result is constant but we do need to set it first.  To do that
 *  we cast away const for the class (it's either this or making
 *  everything mutable)
 *//*-----------------------------------------------------------------*/

const LevelData<FArrayBox>&
MappedLevelData::getU(const int a_numValidGhost,
                      const int a_numInvalidGhost) const
{
  return const_cast<MappedLevelData*>(this)->getU(a_numValidGhost,
                                                  a_numInvalidGhost);
}

/*--------------------------------------------------------------------*/
/** \param[in]  a_numInvalidGhost
 *                      Number of invalid ghost cells to compute \<U\>
 *                      in, by interpolating from a coarser level (if
 *                      one exisits).
 *                      If (a_numInvalidGhost == 0), then \<U\> will
 *                      "only" be determined in the valid cells.
 *                      However, note that \<JU\> will have to be
 *                      interpolated into 1 layer of invalid ghosts
 *                      cells and this will have the side effect of
 *                      filling m_numGhostU layers of invalid ghost
 *                      cells with \<U\>!
 *                      If (a_numInvalidGhost < 0) then only the valid
 *                      cells of <JU> will used (\<JU\> will not be
 *                      computed in invalid ghost cells).  This means
 *                      that \<U\> will not be properly updated in
 *                      valid cells adjacent to the coarse-fine
 *                      interface.  However, sometimes this is useful
 *                      if we are only trying to load a coarsened
 *                      representation of the next finer level (proper
 *                      nesting avoids the cells adjacent to the
 *                      interface with the next coarser level.
 *//*-----------------------------------------------------------------*/

LevelData<FArrayBox>&
MappedLevelData::getU(const int a_numValidGhost,
                      const int a_numInvalidGhost)
{
  CH_assert(m_defined);
  CH_assert(a_numValidGhost >= 0 && a_numValidGhost <= m_numGhostU);
  CH_assert(a_numInvalidGhost >= -1 && a_numInvalidGhost <= m_numGhostU);

//--Fill invalid ghosts (this will also fill required invalid ghosts in JU)

  if (hasCoarser() && (a_numInvalidGhost > m_invalidGhostStatusU))
    {
      fillInvalidGhosts();
    }

//--Fill valid cells

  // Do we need to fill?
  bool fillValidCells = false;
  if (m_validStatusU < 0)
    {
      fillValidCells = true;      // Obvious case
    }
  if (m_validCFLayerStatusU > 0)  // Some layers of valid cells near the CF
    {                             // inteface have not been filled
      if (m_validCFLayerStatusU + a_numInvalidGhost > 0)
                                  // We want a_numInvalidGhost layers filled
        {
          m_validCFLayerStatusU = -std::min(0, a_numInvalidGhost);
          fillValidCells = true;
        }
    }
  else if ((m_validStatusU < 0) && (a_numInvalidGhost < 0))
    {
      // Only if m_validStatusU < 0 do we initially set m_validCFLayerStatusU
      m_validCFLayerStatusU = -a_numInvalidGhost;
    }

  if (fillValidCells)
    {
      // To compute <U>, we generally require 1 ghost cell of <JU> everywhere.
      // But if a_numInvalidGhost == -1, we don't need to find invalid ghost
      // data for <JU>
      int numInvalidJU = std::min(a_numInvalidGhost+1, 1);
      // Note: block boundaries are treated the same as domain boundaries
      m_LGMP->computeValidU(rawU(), getJU(1, numInvalidJU));
      m_validStatusU = 0;
    }

//--Fill valid ghosts

  if (a_numValidGhost > m_validStatusU)
    {
      if (m_useJUExchangerForU && a_numValidGhost <= m_numGhostJU)
        {
          rawU().exchange(m_JUExchangeCopier);
          m_validStatusU = m_numGhostJU;
        }
      else
        {
          rawU().exchange(m_UExchangeCopier);
          m_validStatusU = m_numGhostU;
        }
    }

  // Extra-block ghost cells (which are a type of valid ghost) are tracked
  // separately from valid ghosts because different exchangers could be used
  // for the valid ghosts above.
  if (a_numValidGhost > m_extraBlockStatusU)
    {
      m_LGMP->multiblockExchangeU(rawU(), m_velocityIntv);
      m_extraBlockStatusU = m_numGhostU;
    }

  return rawU();
}

/*--------------------------------------------------------------------*/
//
/** The result is constant but we do need to set it first.  To do that
 *  we cast away const for the class (it's either this or making
 *  everything mutable)
 *//*-----------------------------------------------------------------*/

const LevelData<FArrayBox>&
MappedLevelData::getJU(const int a_numValidGhost,
                       const int a_numInvalidGhost) const
{
  return const_cast<MappedLevelData*>(this)->getJU(a_numValidGhost,
                                                   a_numInvalidGhost);
}

/*--------------------------------------------------------------------*/
/**
 *  \note
 *  <ul>
 *    <li> *** WARNING *** Ghost cells of \<JU\> are _never_ filled
 *         across block boundaries.
 *         Internal notes: If this is ever needed, we can exchange
 *         \<U\> and use derivatives of the interpolating polynomial
 *         to get Grad(U).  One-sided gradients of \<J\> would be used
 *         as for filling invalid ghosts between levels.  The method
 *         gradAvgXipExact to get gradients of Xi from \<Xi\> could be
 *         reused from FourthOrderMappedFineInterpF.ChF.  But note
 *         that for Ax=b, Ax is stored as a vector in the multiblock-
 *         exchange technique.  Hence extra vectors would need to be
 *         stored in each dimension for A(dx/dxi_d).
 *  </ul>
 *//*-----------------------------------------------------------------*/

LevelData<FArrayBox>&
MappedLevelData::getJU(const int a_numValidGhost,
                       const int a_numInvalidGhost)
{
  CH_assert(m_defined);
  CH_assert(a_numValidGhost >= 0 && a_numValidGhost <= m_numGhostJU);
  CH_assert(a_numInvalidGhost >= 0 && a_numInvalidGhost <= m_numGhostJU);


//--Fill the invalid ghosts (NOT valid cells of another box at this level)

  if (hasCoarser() && (a_numInvalidGhost > m_invalidGhostStatusJU))
    {
      fillInvalidGhosts();
    }

//--Fill the valid ghosts (valid cells of another box at this level)

  if (a_numValidGhost > m_validStatusJU)
    {
      rawJU().exchange(m_JUExchangeCopier);
      m_validStatusJU = m_numGhostJU;
    }

  return rawJU();
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::invalidateU(const bool a_this2FnInterpDataInvalid)
{
  m_invalidGhostStatusU = 0;
  m_validStatusU = -1;
  m_validCFLayerStatusU = 0;
  m_extraBlockStatusU = 0;
  // Ensure that the finer level space-interpolator knows that
  // any coarse data is now invalid
  if (a_this2FnInterpDataInvalid)
    {
      m_LGMP->invalidateThis2FnInterpolatorCrFnLevData();
    }
}

/*--------------------------------------------------------------------*/
//  Invalidate all of \<U\>
/** \param[in]  a_this2FnInterpDataInvalid
 *              T - interpolation data for a finer level is invalid
 *                  (default)
 *              F - interpolation data is still valid.  This usually
 *                  is required after advance of level 'l'.
 *                  Interpolation for level 'l+1' still depends on
 *                  the old solution and is not affected
 *//*-----------------------------------------------------------------*/

void
MappedLevelData::invalidateJUGhosts()
{
  m_invalidGhostStatusJU = 0;
  m_validStatusJU = 0;
}

/*--------------------------------------------------------------------*/
//  Invalidate all ghost data and all of \<U\>
/** \param[in]  a_this2FnInterpDataInvalid
 *              T - interpolation data for a finer level is invalid
 *                  (default)
 *              F - interpolation data is still valid.  This usually
 *                  is required after advance of level 'l'.
 *                  Interpolation for level 'l+1' still depends on
 *                  the old solution and is not affected
 *//*-----------------------------------------------------------------*/

void
MappedLevelData::invalidate(const bool a_this2FnInterpDataInvalid)
{
  invalidateU(a_this2FnInterpDataInvalid);
  invalidateJUGhosts();
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::copyJUNewToOld()
{
  CH_assert(m_defined);
  for (DataIterator dit = rawJU().dataIterator(); dit.ok(); ++dit)
    {
      getJUOld()[dit()].copy(rawJU()[dit()]);
    }
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::copyJUOldToNew()
{
  CH_assert(m_defined);
  getJUOld().copyTo(getJUOld().interval(),
                    rawJU(),
                    rawJU().interval());
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::swapJUNewAndOld()
{
  m_indexJUNew = m_indexJUOld;
  m_indexJUOld = !m_indexJUNew;
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::deallocate()
{
  if (m_allocU >= 0)
    {
      delete m_ULevP;
    }
  if (m_allocJUNew >= 0)
    {
      delete m_JULevP[m_allocJUNew];
    }
  if (m_allocJUOld >= 0)
    {
      delete m_JULevP[m_allocJUOld];
    }
  m_allocU = -1;
  m_allocJUNew = -1;
  m_allocJUOld = -1;
}

/*--------------------------------------------------------------------*/

void
MappedLevelData::fillInvalidGhosts()
{
  // Always assume coarsened-fine data <U> has been preset - the interpolator
  // should catch cases where the coarsened-fine data becomes invalid (see
  // invalidateU() in this class)
  m_LGMP->fillFineGhostCells(rawU(), rawJU());
  m_invalidGhostStatusU = m_numGhostU;
  m_invalidGhostStatusJU = m_numGhostJU;
}

#include "NamespaceFooter.H"
