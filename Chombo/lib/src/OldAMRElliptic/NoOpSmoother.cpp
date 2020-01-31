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

#include "NoOpSmoother.H"
#include "DotProduct.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

NoOpSmoother::NoOpSmoother()
{
}

NoOpSmoother::~NoOpSmoother()
{
}

BaseBottomSmoother* NoOpSmoother::new_bottomSmoother() const
{
  NoOpSmoother* newsmoother = new NoOpSmoother();

  if (newsmoother == NULL)
    {
      MayDay::Error("Out of Memory in NoOpSmoother::new_bottomSmoother");
    }

  return static_cast<BaseBottomSmoother*>(newsmoother);
}

/***********************/
// True to its name, nothing gets done here
/***********************/
void NoOpSmoother::doBottomSmooth(LevelData<FArrayBox>&           a_phi,
                                  const LevelData<FArrayBox>& a_rhs,
                                  LevelOp*                    a_levelopPtr)
{
  // look at me doing nothing here...
}
#include "NamespaceFooter.H"
