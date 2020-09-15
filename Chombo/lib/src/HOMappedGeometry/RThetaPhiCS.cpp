#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RThetaPhiCS.H"

#include "NamespaceHeader.H"

// Constructor
RThetaPhiCS::RThetaPhiCS(const RealVect& a_dx,
                         const Real      a_stretch,
                         const Real      a_rMin,
                         const Real      a_phiMin)
  :
  m_stretch(a_stretch),
  m_rMin(a_rMin),
  m_phiMin(a_phiMin)
{
  m_dx = a_dx;
}

// Destructor
RThetaPhiCS::~RThetaPhiCS()
{
}

// Given coordinate in mapped space, return its location in real space
RealVect
RThetaPhiCS::realCoord(const RealVect& a_hatXi) const
{
  RealVect Xi = unnormalize(a_hatXi);
  return realCoordUN(Xi);
}

// Given coordinate in real space, return its location in the mapped space
RealVect
RThetaPhiCS::mappedCoord(const RealVect& a_x) const
{
  CH_assert(SpaceDim >= 2 && SpaceDim <= 3);

  RealVect Xi;

  const Real r_xy = sqrt(a_x[0]*a_x[0] + a_x[1]*a_x[1]);
  Xi[1] = acos(a_x[0]/r_xy);
  if (a_x[1] < 0)
    {
      Xi[1] = 2*Pi - Xi[1];
    }

#if (CH_SPACEDIM == 2)
  Xi[0] = r_xy;
#elif (CH_SPACEDIM >= 3)
  Xi[0] = sqrt(r_xy*r_xy + a_x[2]*a_x[2]);
  Xi[2] = acos(a_x[2]/Xi[0]);
#endif

  RealVect hatXi = normalize(Xi);
  return hatXi;
}

// Calculate the derivative of each coordinate vector
Real RThetaPhiCS::dXdXi(const RealVect& a_hatXi, int a_dirX, int a_dirXi) const
{
  CH_assert(a_dirX  >= 0 && a_dirX  < SpaceDim);
  CH_assert(a_dirXi >= 0 && a_dirXi < SpaceDim);

  const RealVect Xi = unnormalize(a_hatXi);
  const RealVect X  = realCoordUN(Xi);

  // Notational convenience
  const Real r     = Xi[0];
  const Real theta = Xi[1];
#if   (CH_SPACEDIM == 2)
  const Real phi   = 0.5*Pi;
#elif (CH_SPACEDIM >= 3)
  const Real phi   = Xi[2];
#elif (CH_SPACEDIM == 1)
  // add this so that code compiles in 1D, not that it will be used anyway
  const Real phi   = 0.;
#endif

  switch (4*a_dirX + a_dirXi)
    {
    case 4*0 + 0:  // dx/dr
      return m_stretch*sin(phi)*cos(theta);
      break;
    case 4*0 + 1:  // dx/dtheta
      return -r*sin(phi)*sin(theta)*2*Pi;
      break;
    case 4*0 + 2:  // dx/dphi
      return r*cos(phi)*(Pi + 2*m_phiMin)*cos(theta);
      break;
    case 4*1 + 0:  // dy/dr
      return m_stretch*sin(phi)*sin(theta);
      break;
    case 4*1 + 1:  // dy/dtheta
      return r*sin(phi)*cos(theta)*2*Pi;
      break;
    case 4*1 + 2:  // dy/dphi
      return r*cos(phi)*(Pi + 2*m_phiMin)*sin(theta);
      break;
    case 4*2 + 0:  // dz/dr
      return m_stretch*cos(theta);
      break;
    case 4*2 + 1:  // dz/dtheta
      return 0.;
      break;
    case 4*2 + 2:  // dz/dphi
      return -r*sin(phi)*(Pi + 2*m_phiMin);
      break;
    }
  return 0;
}

// Normalize the mapped coordinates \f$(0 <= \hat{\xi} <= 1)\f$
RealVect
RThetaPhiCS::normalize(const RealVect & a_Xi) const
{
  RealVect hatXi(
    D_DECL6(a_Xi[0]/m_stretch - m_rMin,
           a_Xi[1]/(2*Pi),
            (a_Xi[2] - m_phiMin)/(Pi - 2*m_phiMin),0.,0.,0.));
  return hatXi;
}

// Unnormalized the mapped coordinates
RealVect
RThetaPhiCS::unnormalize(const RealVect & a_hatXi) const
{
  RealVect Xi(
    D_DECL6(m_stretch*(a_hatXi[0] + m_rMin),
           2*Pi*a_hatXi[1],
            a_hatXi[2]*(Pi - 2*m_phiMin) + m_phiMin,0.,0.,0.));
  return Xi;
}

// Return location in real space of an unnormalized Coord. from mapped space
RealVect
RThetaPhiCS::realCoordUN(const RealVect& a_Xi) const
{
#if (CH_SPACEDIM == 2)
  const Real phi = 0.5*Pi;
#elif (CH_SPACEDIM >= 3)
  const Real phi = a_Xi[2];
#elif (CH_SPACEDIM == 1)
  // add this so that code compiles in 1D, not that it will be used anyway
  const Real phi = 0.;
#endif
  RealVect xyz(
    D_DECL6(a_Xi[0]*sin(phi)*cos(a_Xi[1]),
           a_Xi[0]*sin(phi)*sin(a_Xi[1]),
            a_Xi[0]*cos(phi),0.,0.,0.));
  return xyz;
}

// -- begin factory implementations ---------------------------

RThetaPhiCSFactory::RThetaPhiCSFactory(const Real a_stretch,
                                       const Real a_rMin,
                                       const Real a_phiMin)
  :
  m_stretch(a_stretch),
  m_rMin(a_rMin),
  m_phiMin(a_phiMin)
{
}

NewCoordSys*
RThetaPhiCSFactory::getCoordSys(const ProblemDomain& a_domain,
                                const RealVect&      a_dx) const
{
  RThetaPhiCS* newCSPtr = new RThetaPhiCS(a_dx, m_stretch, m_rMin, m_phiMin);
  return static_cast<NewCoordSys*>(newCSPtr);
}

#include "NamespaceFooter.H"
