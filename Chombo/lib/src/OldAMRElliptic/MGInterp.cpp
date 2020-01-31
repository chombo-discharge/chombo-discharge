#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "REAL.H"
#include "DataIterator.H"
#include "LevelMGF_F.H"

#include "MGInterp.H"
#include "NamespaceHeader.H"

MGInterp::MGInterp()
  :m_isDefined(false)
{
}

MGInterp::~MGInterp()
{
}

MGInterp::MGInterp(const DisjointBoxLayout& a_fineDomain,
                   int                      a_numcomps,
                   int                      a_refRatio,
                   const Box&               a_problemDomain)
  :m_isDefined(false)
{
  ProblemDomain finePhysdomain(a_problemDomain);

  define(a_fineDomain, a_numcomps, a_refRatio, finePhysdomain);
}

MGInterp::MGInterp(const DisjointBoxLayout& a_fineDomain,
                   int                      a_numcomps,
                   int                      a_refRatio,
                   const ProblemDomain&     a_problemDomain)
  :m_isDefined(false)
{
  define(a_fineDomain, a_numcomps, a_refRatio, a_problemDomain);
}

void MGInterp::define(const DisjointBoxLayout& a_fineDomain,
                      int                      a_numcomps,
                      int                      a_refRatio,
                      const Box&               a_problemDomain)
{
  ProblemDomain finePhysdomain(a_problemDomain);

  define(a_fineDomain, a_numcomps, a_refRatio, finePhysdomain);
}

void MGInterp::define(const DisjointBoxLayout& a_fineDomain,
                      int                      a_numcomps,
                      int                      a_refRatio,
                      const ProblemDomain&     a_problemDomain)
{
  m_refRatio = a_refRatio;
  m_problemDomain = a_problemDomain;
  m_grids = a_fineDomain;

  // create the work array
  DisjointBoxLayout coarsenedFineDomain;

  coarsen(coarsenedFineDomain, a_fineDomain, m_refRatio);

  m_coarsenedFineData.define(coarsenedFineDomain,
                             a_numcomps,
                             IntVect::Zero);

  m_isDefined = true;
}

// interpolate from coarse level to fine level
void MGInterp::interpToFine(LevelData<FArrayBox>&       a_fineData,
                            const LevelData<FArrayBox>& a_coarseData)
{
  CH_assert(m_isDefined);
  a_coarseData.copyTo(a_coarseData.interval(),
                      m_coarsenedFineData,
                      m_coarsenedFineData.interval());

  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& fine = a_fineData[dit()];
      const FArrayBox& coarsenedFine = m_coarsenedFineData[dit()];

      Box crseBox = coarsenedFine.box();
      Box fineBox = m_grids.get(dit());

      fineBox.coarsen(m_refRatio);

      CH_assert(fineBox == crseBox);

      Box nrefbox(IntVect::Zero,
                  (m_refRatio-1)*IntVect::Unit);

      FORT_INTERPMG(CHF_FRA(fine),
                    CHF_FRA(coarsenedFine),
                    CHF_BOX(crseBox),
                    CHF_CONST_INT(m_refRatio),
                    CHF_BOX(nrefbox));
    }
}

bool MGInterp::isDefined() const
{
  return m_isDefined;
}
#include "NamespaceFooter.H"
