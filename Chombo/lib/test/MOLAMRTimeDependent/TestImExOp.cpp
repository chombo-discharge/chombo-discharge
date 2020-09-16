#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TestImExOp.H"

#include "NamespaceHeader.H"

const Real TestImExOp::s_cI = -6;
const Real TestImExOp::s_cE = .5;

TestImExOp::TestImExOp()
{
  m_isDefined = false;
}

TestImExOp::~TestImExOp()
{
}

void
TestImExOp::define(const TestSolnData&   a_state,
                    Real a_dt,
                    Real a_dtscale)
{
  m_dt = a_dt;
  m_stage0 = -1;
  m_stage1 = -1;
  m_dtscale = a_dtscale;
  m_isDefined = true;
}

void
TestImExOp::resetDt(Real a_dt)
{
  CH_assert(isDefined());
  m_dt = a_dt;
}

void
TestImExOp::stage0stage1(int a_stage0,
                         int a_stage1)
{
  CH_assert(isDefined());
  m_stage0 = a_stage0;
  m_stage1 = a_stage1;
}


void
TestImExOp::explicitOp(TestRhsData& a_result, Real a_time, 
                       const TestSolnData&  a_state,
                       int a_stage)
{
  CH_TIMERS("TestImExOp::explicitOp");
  CH_assert(isDefined());

  const LevelData<FArrayBox>& stateData = a_state.data();

  const DisjointBoxLayout& layout = stateData.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();

  LevelData<FArrayBox>& resultData = a_result.data();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& stateDataFab = stateData[dit];
      FArrayBox& resultDataFab = resultData[dit];

      // This is just a simple test for the time integration
      Real coef = s_cE + 1.0/(1.0 + a_time);
      resultDataFab.copy(stateDataFab);
      resultDataFab.mult(coef);
    }
}

void
TestImExOp::implicitOp(TestRhsData&             a_result,
                       Real                    a_time,
                       const TestSolnData&  a_state,
                       int a_stage)
{
  CH_TIMERS("TestImExOp::implicitOp");
  CH_assert(isDefined());

  const LevelData<FArrayBox>& stateData = a_state.data();
  const DisjointBoxLayout& grids = stateData.disjointBoxLayout();

  a_result.define(a_state);
  LevelData<FArrayBox>& resultData = a_result.data();

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& stateDataFab = stateData[dit];
      FArrayBox& resultDataFab = resultData[dit];
      resultDataFab.setVal(0.);

      // Just a simple test op for the time integration
      resultDataFab.copy(stateDataFab);
      resultDataFab.mult(s_cI);
    }
}

void
TestImExOp::solve(TestSolnData&   a_soln,
                   const TestRhsData&  a_rhs,
                   Real               a_time,
                   int a_stage)
{
  CH_TIMERS("TestImExOp::solve");
  CH_assert(isDefined());

  LevelData<FArrayBox>& solnData = a_soln.data();
  const DisjointBoxLayout& grids = solnData.disjointBoxLayout();
  const LevelData<FArrayBox>& rhsData = a_rhs.data();

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& solnDataFab = solnData[dit];
      const FArrayBox& rhsDataFab = rhsData[dit];

      // Just a simple test op solve for the time integration
      Real coef = 1.0 / (1.0 - m_dtscale*s_cI*m_dt);
      solnDataFab.copy(rhsDataFab);
      solnDataFab.mult(coef);
    }
}

#include "NamespaceFooter.H"
