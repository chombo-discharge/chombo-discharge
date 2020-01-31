#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DFMartin, Sun, May 5, 2002

#include "GSRBSmoother.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

GSRBSmoother::GSRBSmoother()
{
  // set default value for m_numSmooth
  m_numSmooth = 60;
}

GSRBSmoother::~GSRBSmoother()
{
}

BaseBottomSmoother* GSRBSmoother::new_bottomSmoother() const
{
  GSRBSmoother* newsmoother = new GSRBSmoother();

  if (newsmoother == NULL)
    {
      MayDay::Error("Out of Memory in GSRBSmoother::new_bottomSmoother");
    }

  return static_cast<BaseBottomSmoother*>(newsmoother);
}

/***********************/
// True to its name, we do a bunch of GSRB here
/***********************/
void GSRBSmoother::doBottomSmooth(LevelData<FArrayBox>&       a_phi,
                                  const LevelData<FArrayBox>& a_rhs,
                                  LevelOp*                    a_levelopPtr)
{

  for (int n = 0; n < m_numSmooth; n++)
    {
      a_levelopPtr->smooth(a_phi, a_rhs);
    }
}

void GSRBSmoother::setNumSmooth(int a_numSmooth)
{
  m_numSmooth = a_numSmooth;
}

int GSRBSmoother::numSmooth() const
{
  return m_numSmooth;
}
#include "NamespaceFooter.H"
