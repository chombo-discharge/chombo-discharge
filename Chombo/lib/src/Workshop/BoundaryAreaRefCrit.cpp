#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "REAL.H"
#include "RealVect.H"
#include "Box.H"
#include "IntVect.H"
#include "EBISBox.H"
#include "GeometryService.H"
#include "WrappedGShop.H"
#include "LSquares.H"
#include "BaseIF.H"
#include "AMRIO.H"
#include "CutCellMoments.H"
#include "BoundaryAreaRefCrit.H"
#include "BoundaryAreaRefCrit.H"
#include "NamespaceHeader.H"


bool
BoundaryAreaRefCrit::
refineHere(const IrregNode& a_node, const IntVect& a_iv, const Real& a_dx) const
{
  Real bndryArea     = a_node.m_EBMoments[IvSpaceDim::Zero];
  Real bndryAreaFrac = bndryArea;
  for(int idir = 0; idir < SpaceDim-1; idir++)
    {
      bndryAreaFrac /= a_dx;
    }
  bool retval = (bndryAreaFrac > m_thresholdBA);
  return retval;
}


#include "NamespaceFooter.H"
