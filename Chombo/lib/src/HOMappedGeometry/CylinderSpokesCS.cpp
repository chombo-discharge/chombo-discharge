#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CylinderSpokesCS.H"
#include "CylinderSpokesBlockCS.H"

#include "NamespaceHeader.H"

CylinderSpokesCS::CylinderSpokesCS()
{
}

CylinderSpokesCS::~CylinderSpokesCS()
{
  if (m_gotCoordSysVect)
    {
      for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
        {
          delete m_coordSysVect[iblock];
        }
    }
}

void CylinderSpokesCS::define(const ProblemDomain& a_levelDomain,
                                   const RealVect& a_dx)
{
  CylinderCS::define(a_levelDomain, a_dx);
  for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
    {
      const Box& bx = m_mappingBlocks[iblock];
      m_coordSysVect[iblock] =
        new CylinderSpokesBlockCS(iblock, a_dx, bx);
    }
  m_gotCoordSysVect = true;
}

// ----------------------------------------------------------------------------

MultiBlockCoordSys*
CylinderSpokesCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                          const RealVect& a_dx) const
{
  CylinderSpokesCS* coordSysPtr = new CylinderSpokesCS();
  coordSysPtr->define(a_levelDomain, a_dx);
  setAllPhysical(coordSysPtr);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
