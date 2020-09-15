#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TestSolnData.H"

#include "NamespaceHeader.H"

TestSolnData::TestSolnData()
{
  m_isDefined = false;
  m_accum = NULL;
}

TestSolnData::~TestSolnData()
{
}

void
TestSolnData::define(/// box layout at this level
                        const DisjointBoxLayout&    a_layout,
                        /// number of conserved components
                        const int                   a_nComp,
                        /// ghost vector
                        const IntVect&              a_ghostVect
                        )
{
  m_grids = a_layout;
  m_ghostVect = a_ghostVect;
  m_nComp = a_nComp;
  m_data.define(m_grids, m_nComp, m_ghostVect);
  m_isDefined = true;
}

void
TestSolnData::define(const TestSolnData& a_state)
{
  const LevelData<FArrayBox>& srcData = a_state.data();
  m_grids = srcData.disjointBoxLayout();
  m_nComp = srcData.nComp();
  m_ghostVect = srcData.ghostVect();
  m_data.define(m_grids, m_nComp, m_ghostVect);
  m_isDefined = true;
}

/// Constructor that aliases an incoming LevelData<FArrayBox>
void
TestSolnData::aliasData(
    LevelData<FArrayBox>& a_data,
    LevelData<FArrayBox>* a_accum)
{
  Interval aliasInt(0,a_data.nComp()-1);
  aliasLevelData<FArrayBox>(m_data, &a_data, aliasInt);
  m_isDefined = true;

  if (a_accum != NULL)
    m_accum = a_accum;
}

void
TestSolnData::copy(const TestSolnData& a_state)
{
  CH_TIMERS("TestSolnData::copy(TestSolnData)");
  const LevelData<FArrayBox>& srcData = a_state.data();
  DataIterator dit = srcData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& srcDataFab = srcData[dit];
      FArrayBox& dataFab = m_data[dit];
      dataFab.copy(srcDataFab);
    }
}

/*
void
TestSolnData::copy(const TestRhsData& a_value)
{
  CH_TIMERS("TestSolnData::copy(TestRhsData)");
  const LevelData<FArrayBox>& srcData = a_value.data();
  DataIterator dit = srcData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& srcDataFab = srcData[dit];
      FArrayBox& dataFab = m_data[dit];
      dataFab.copy(srcDataFab);
    }
}
*/
  
void
TestSolnData::zero()
{
  CH_TIMERS("TestSolnData::zero");
  DataIterator dit = m_data.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& dataFab = m_data[dit];
      dataFab.setVal(0);
    }
}

void
TestSolnData::increment(const TestRhsData& a_rhs,
                      Real a_factor,
                      bool a_updateFluxReg)
{
  CH_TIMERS("TestSolnData::increment");
  CH_assert(m_isDefined);
  const LevelData<FArrayBox>& rhsData = a_rhs.data();
  DataIterator dit = rhsData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const FArrayBox& rhsDataFab = rhsData[dit];
    FArrayBox& dataFab = m_data[dit];
    dataFab.plus(rhsDataFab, a_factor);

    if (m_accum != NULL)
    {
      FArrayBox& accumFab = (*m_accum)[dit];
      accumFab.plus(rhsDataFab, a_factor);
    }
  }
}

#include "NamespaceFooter.H"
