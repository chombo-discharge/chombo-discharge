#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFSIMPLESOLVER_H_
#define _MFSIMPLESOLVER_H_

#include "LevelData.H"
#include "LinearSolver.H"
#include "MultiGrid.H"

#include "MFCellFAB.H"
#include "NamespaceHeader.H"

class MFSimpleSolver : public LinearSolver<LevelData<MFCellFAB> >
{
public:
  int m_verbosity;

  MFSimpleSolver();

  virtual ~MFSimpleSolver();

  virtual void setHomogeneous(bool a_homogeneous);

  virtual void define(LinearOp<LevelData<MFCellFAB> >* a_operator,
                      bool                             a_homogeneous);

  virtual void setNumSmooths(const int& a_numSmooths);

  virtual void solve(LevelData<MFCellFAB>&       a_phi,
                     const LevelData<MFCellFAB>& a_rhs);

private:
  bool m_isDefined;

  MGLevelOp<LevelData<MFCellFAB> >* m_operator;

  bool m_homogeneous;
  int  m_numSmooths;
};

#include "NamespaceFooter.H"
#endif
