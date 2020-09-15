#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOUNDARYAREAREFCRIT_H_
#define _BOUNDARYAREAREFCRIT_H_

#include "WrappedGShop.H"
#include "RealVect.H"
#include "IntVect.H"
#include "IrregNode.H"
#include "NamespaceHeader.H"

///base class for special refinement
/**
   Refines if the nodes boundary area FRACTION is greater  a threshold
 */
class BoundaryAreaRefCrit: public WGSRefinementCriterion
{
public:
  ///
  virtual bool refineHere(const IrregNode& a_node, const IntVect& a_iv, const Real& a_dx) const;

  ///
  BoundaryAreaRefCrit(Real a_thresholdBA)
  {
    m_thresholdBA = a_thresholdBA;
  }

  ///
  virtual ~BoundaryAreaRefCrit()
  {
  }

private:
  Real m_thresholdBA;


  BoundaryAreaRefCrit()
  {
  }
};
#include "NamespaceFooter.H"
#endif
