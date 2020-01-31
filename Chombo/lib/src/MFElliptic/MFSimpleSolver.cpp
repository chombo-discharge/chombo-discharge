#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFSimpleSolver.H"
#include "NamespaceHeader.H"

MFSimpleSolver::MFSimpleSolver()
{
  m_isDefined = false;

  m_operator = NULL;
  m_homogeneous = true;

  m_numSmooths = 8;

  m_verbosity = 1;
}

MFSimpleSolver::~MFSimpleSolver()
{
}

void MFSimpleSolver::setHomogeneous(bool a_homogeneous)
{
  m_homogeneous = a_homogeneous;
}

void MFSimpleSolver::define(LinearOp<LevelData<MFCellFAB> >* a_operator,
                            bool                             a_homogeneous)
{
  m_isDefined = true;

  m_operator = dynamic_cast <MGLevelOp<LevelData<MFCellFAB> >* > (a_operator);
  if (m_operator == NULL)
  {
    MayDay::Error("MFSimpleSolver::define - operator not a MGLevelOp");
  }

  m_homogeneous = a_homogeneous;
}

void MFSimpleSolver::setNumSmooths(const int& a_numSmooths)
{
  m_numSmooths = a_numSmooths;
}

void MFSimpleSolver::solve(LevelData<MFCellFAB>&       a_phi,
                           const LevelData<MFCellFAB>& a_rhs)
{
  CH_assert(m_isDefined);

  m_operator->relax(a_phi,a_rhs,m_numSmooths);
}
#include "NamespaceFooter.H"
