#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GrowInBlock.H"
#include "MBMiscUtil.H"

#include "NamespaceHeader.H"

GrowInBlock::GrowInBlock(
                         LevelGridMetrics* a_levelGridMetricsPtr,
                         int a_ghost)
{
  m_levelGridMetricsPtr = a_levelGridMetricsPtr;
  m_ghost = a_ghost;
}

GrowInBlock::~GrowInBlock()
{
}

Box
GrowInBlock::operator() (const Box& a_inputBox)
{
  const ProblemDomain& domain =
    m_levelGridMetricsPtr->blockDomain(a_inputBox, m_ghost);
  Box bx = grow(a_inputBox, m_ghost);
  bx &= domain;
  return bx;
}

#include "NamespaceFooter.H"
