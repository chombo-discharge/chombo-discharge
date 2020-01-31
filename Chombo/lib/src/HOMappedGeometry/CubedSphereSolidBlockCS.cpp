#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSphereSolidBlockCS.H"
// #include <cmath>

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

CubedSphereSolidBlockCS::CubedSphereSolidBlockCS(int a_blockID,
                                                 const RealVect& a_dxVect,
                                                 const Box& a_bx)
{
  CH_assert((a_blockID >= 0) && (a_blockID < 7));
  m_blockID = a_blockID;
  m_dx = a_dxVect;
  m_mappedLo = RealVect(surroundingNodes(a_bx).smallEnd()) * m_dx;
  m_signCoord = RealVect::Unit;
  switch (m_blockID)
    {
      // For non-CENTRAL blocks, set xyz[m_rCoord] = G * m_signCoord,
      // then xyz[m_aCoord] *= a and xyz[m_bCoord] *= b.
    case XPOS:
      { // x = G, y = b*G, z = a*G
        m_rCoord = 0;
        m_aCoord = 2;
        m_bCoord = 1;
        break;
      }
    case YPOS:
      { // x = -b*G, y = G, z = a*G
        m_rCoord = 1;
        m_aCoord = 2;
        m_bCoord = 0;
        m_signCoord[0] = -1.;
        break;
      }
    case XNEG:
      { // x = -G, y = -b*G, z = a*G
        m_rCoord = 0;
        m_aCoord = 2;
        m_bCoord = 1;
        m_signCoord[0] = -1.;
        m_signCoord[1] = -1.;
        break;
      }
    case YNEG:
      { // x = b*G, y = -G, z = a*G
        m_rCoord = 1;
        m_aCoord = 2;
        m_bCoord = 0;
        m_signCoord[1] = -1.;
        break;
      }
    case ZPOS:
      { // x = -a*G, y = b*G, z = G
        m_rCoord = 2;
        m_aCoord = 0;
        m_bCoord = 1;
        m_signCoord[0] = -1.;
        break;
      }
    case ZNEG:
      { // x = a*G, y = b*G, z = -G
        m_rCoord = 2;
        m_aCoord = 0;
        m_bCoord = 1;
        m_signCoord[2] = -1.;
        break;
      }
    }
}

CubedSphereSolidBlockCS::~CubedSphereSolidBlockCS()
{
}

RealVect
CubedSphereSolidBlockCS::realCoord(const RealVect& a_Xi) const
{
  // Convert a_Xi to coordinates abc in [-1, 1].
  // Mapped coordinates of this block are m_mappedLo + [0:1].
  RealVect abc = (2. * (a_Xi - m_mappedLo) - RealVect::Unit) * m_signCoord;
  RealVect X;
  if (m_blockID == CENTRAL)
    {
      X = m_r0 * abc;
    }
  else
    {
      // Get the correct components a, b, c from abc.
      Real a = abc[m_aCoord];
      Real b = abc[m_bCoord];
      Real c = abc[m_rCoord];
      Real r = (m_r0 * (1. - c) + m_r1 * (1. + c)) / 2.;
      Real E = 1. + a*a + b*b;
      Real F2 = ((m_r1 - r) + (r - m_r0) * E) / (m_r1 - m_r0);
      if (F2 <= 0.) 
        {
          pout() << "Bad call to realCoord for block " << m_blockID << ":  "
                 << a_Xi << endl;
        }
      Real F = sqrt(F2);
      Real G = r / F;
      X = G * m_signCoord;
      X[m_aCoord] *= a;
      X[m_bCoord] *= b;
    }
  return X;
}


RealVect
CubedSphereSolidBlockCS::mappedCoord(const RealVect& a_x) const
{
  RealVect abc;
  if (m_blockID == CENTRAL)
    {
      abc = a_x / m_r0;
    }
  else
    {
      // xUnsigned contains G, a*G, b*G in some order.
      RealVect xUnsigned = a_x * m_signCoord;
      Real G = xUnsigned[m_rCoord];
      Real a = xUnsigned[m_aCoord] / G;
      Real b = xUnsigned[m_bCoord] / G;
      Real E = 1. + a*a + b*b;
      // Note that contrary to Diener et al. 2007, it is NOT generally
      // true that r = sqrt(x^2 + y^2 + z^2).
      // We need to solve a quadratic equation for r.
      Real aQuadratic = (m_r1 - m_r0);
      Real bQuadratic = (1. - E) * G*G;
      Real cQuadratic = (m_r0 * E - m_r1) * G*G;
      Real discrim = bQuadratic*bQuadratic - 4.*aQuadratic*cQuadratic;
      // Take only the positive root.
      Real r = (-bQuadratic + sqrt(discrim)) / (2.*aQuadratic);
      Real c = (2.*r - m_r1 - m_r0) / (m_r1 - m_r0);
      if (discrim < 0.) 
        {
          pout() << "Bad call to mappedCoord for block " << m_blockID << ":  "
                 << a_x << ", a = " << a << ", b = " << b << ", c = " << c << endl;
        }
      // Set the correct components of abc to match the order in xi.
      abc[m_aCoord] = a;
      abc[m_bCoord] = b;
      abc[m_rCoord] = c;
    }
  // Convert abc in [-1, 1] to mapped coordinates of this block.
  // Mapped coordinates of this block are m_mappedLo + [0:1].
  RealVect xi = m_mappedLo + (m_signCoord * abc + RealVect::Unit)/2.;
  return xi;
}

Real
CubedSphereSolidBlockCS::dXdXi(const RealVect& a_Xi,
                               int a_dirX,
                               int a_dirXi) const
{
  Real retval;
  if (m_blockID == CENTRAL)
    {
      if (a_dirX == a_dirXi)
        {
          retval = 2. * m_r0;
        }
      else
        {
          retval = 0.;
        }
    }
  else
    {
      // RealVect xi = m_mappedLo + (m_signCoord * abc + 1.)/2.;
      RealVect abc = (2.*(a_Xi - m_mappedLo) - RealVect::Unit) * m_signCoord;
      Real a = abc[m_aCoord];
      Real b = abc[m_bCoord];
      Real c = abc[m_rCoord];
      Real r = (m_r0 * (1. - c) + m_r1 * (1. + c)) / 2.;
      Real E = 1 + a*a + b*b;
      Real F2 = ((m_r1 - r) + (r - m_r0) * E) / (m_r1 - m_r0);
      if (F2 <= 0.) 
        {
          pout() << "Bad call to dXdXi for block " << m_blockID << ":  "
                 << a_Xi << endl;
        }
      Real F = sqrt(F2);
      Real G = r / F;
      Real H = 1. / ((m_r0 - m_r1) * (F*F*F));
      if (a_dirXi == m_rCoord)
        { // Find derivative of xyz[a_dirX] with respect to c.
          retval = 1./(F*F*F) + (a*a + b*b)*(m_r0 - r/2.)*H;
          if (a_dirX == m_aCoord) retval *= a;
          if (a_dirX == m_bCoord) retval *= b;
          // Now retval contains dG/dr or d(a*G)/dr or d(b*G)/dr.
          retval *= m_signCoord[a_dirX];
          // Now retval contains dX[a_dirX]/dr.
          Real drdc = (m_r1 - m_r0)/2.;
          retval *= drdc;
          // Now retval contains dX[a_dirX]/dc.
        }
      else
        { // Find derivative of xyz[a_dirX] with respect to a or b.
          retval = r*(r - m_r0) * H;
          if (a_dirX == m_aCoord) retval *= a;
          if (a_dirX == m_bCoord) retval *= b;
          if (a_dirXi == m_aCoord) retval *= a;
          if (a_dirXi == m_bCoord) retval *= b;
          if (a_dirX == a_dirXi) retval += G;
          // Now retval contains one of d{G, a*G, b*G}/d{a, b}.
          retval *= m_signCoord[a_dirX];
          // Now retval contains dX[a_dirX]/a or dX[a_dirX]/b.
        }
      // Now retval contains one of dX[a_dirX]/d{a, b, c}.
      retval *= 2. * m_signCoord[a_dirXi];
      // Now retval contains dX[a_dirX]/dXi[a_dirXi].
    }
  return retval;
}

#include "NamespaceFooter.H"
