#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NewCoordSys.H"
#include "CartesianCS.H"

#include "NamespaceHeader.H"

/// constructor
CartesianCS::CartesianCS(const RealVect& a_dx,
                         const RealVect& a_origin,
                         const RealVect& a_stretch)
  :
  m_origin(a_origin),
  m_stretch(a_stretch)
{
  m_dx = a_dx;
}
/**
   Destructor.
*/
CartesianCS::~CartesianCS()
{
}

/// given coordinate in mapped space, return its location in real space
RealVect
CartesianCS::realCoord(const RealVect& a_Xi) const
{
  // simplest possible mapping
  return m_stretch * a_Xi + m_origin;
}

/// given coordinate in real space, return its location in the mapped space
RealVect CartesianCS::mappedCoord(const RealVect& a_x) const
{
  // simplest possible mapping
  return (a_x - m_origin)/m_stretch;
}

// Calculate the derivative of each coordinate vector
Real CartesianCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  Real value = 0.0;
  if (a_dirX == a_dirXi) value = m_stretch[a_dirXi];
  return value;

}

// -- begin factory implementations ---------------------------

CartesianCSFactory::CartesianCSFactory(const RealVect& a_origin,
                                       const RealVect& a_stretch)
  :
  m_origin(a_origin),
  m_stretch(a_stretch)
{
}


NewCoordSys*
CartesianCSFactory::getCoordSys(const ProblemDomain& a_domain,
                                const RealVect& a_dx) const
{
  CartesianCS* newCSPtr = new CartesianCS(a_dx, m_origin, m_stretch);
  return static_cast< NewCoordSys* >(newCSPtr);
}


#include "NamespaceFooter.H"
