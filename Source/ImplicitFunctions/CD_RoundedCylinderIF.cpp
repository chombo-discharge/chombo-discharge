/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RoundedCylinderIF.cpp
  @brief  Implementation of CD_RoundedCylinderIF.H
  @date   Feb. 2021
  @author Robert Marskar
*/

// Chombo includes
#include <PolyGeom.H>
#include <TorusIF.H>
#include <IntersectionIF.H>
#include <TransformIF.H>
#include <PlaneIF.H>
#include <SmoothIntersection.H>
#include <SmoothUnion.H>

// Our includes
#include <CD_CylinderSdf.H>
#include <CD_RoundedCylinderIF.H>
#include <CD_RoundedBoxIF.H>
#include <CD_TorusSdf.H>
#include <CD_NamespaceHeader.H>

RoundedCylinderIF::RoundedCylinderIF(const RealVect a_center1,
                                     const RealVect a_center2,
                                     const Real     a_radius,
                                     const Real     a_curv,
                                     const bool     a_fluidInside)
{
  m_center1     = a_center1;
  m_center2     = a_center2;
  m_length      = (m_center2 - m_center1).vectorLength();
  m_radius      = a_radius;
  m_curv        = a_curv;
  m_fluidInside = a_fluidInside;

  this->makeBaseIF();
}

RoundedCylinderIF::RoundedCylinderIF(const RoundedCylinderIF& a_inputIF)
{
  m_fluidInside = a_inputIF.m_fluidInside;
  m_baseIF      = a_inputIF.m_baseIF;
}

Real
RoundedCylinderIF::value(const RealVect& a_point) const
{
  Real retval = m_baseIF->value(a_point);

  if (m_fluidInside) {
    retval = -retval;
  }

  return retval;
}

BaseIF*
RoundedCylinderIF::newImplicitFunction() const
{
  return (BaseIF*)(new RoundedCylinderIF(*this));
}

void
RoundedCylinderIF::makeBaseIF()
{
#if CH_SPACEDIM == 2
  BaseIF* bif = this->makeBaseIF2D();
#elif CH_SPACEDIM == 3
  BaseIF* bif = this->makeBaseIF3D();
#endif

  // Rotate and translate into place
  TransformIF* transif = new TransformIF(*bif);

  const int      dir = SpaceDim - 1;
  const RealVect up  = BASISREALV(dir);
  if (m_center2[dir] >= m_center1[dir]) {
    transif->rotate(up, m_center2 - m_center1);
    transif->translate(m_center1);
  }
  else {
    transif->rotate(up, m_center1 - m_center2);
    transif->translate(m_center2);
  }

  delete bif;

  m_baseIF = RefCountedPtr<BaseIF>(transif);
}

#if CH_SPACEDIM == 2
BaseIF*
RoundedCylinderIF::makeBaseIF2D()
{
  const RealVect x0 = RealVect::Zero - m_radius * BASISREALV(0);
  const RealVect x1 = RealVect::Zero + m_radius * BASISREALV(0) + m_length * BASISREALV(1);

  return (BaseIF*)(new RoundedBoxIF(x0, x1, m_curv, false));
}
#endif

#if CH_SPACEDIM == 3
BaseIF*
RoundedCylinderIF::makeBaseIF3D()
{

  // TLDR: Construct m_baseIF from a main cylinderk.  on each we put a torus and then a smaller cylinder between everything. Default orientation
  //       is along +z.

  const RealVect up = BASISREALV(SpaceDim - 1);

  const RealVect x0 = RealVect::Zero;
  const RealVect x1 = x0 + m_length * up;

  const RealVect y0 = x0 + m_curv * up;
  const RealVect y1 = x1 - m_curv * up;

  const Real majorRadius = m_radius - m_curv;
  const Real minorRadius = m_curv;

  BaseIF* mainCylinder   = (BaseIF*)(new CylinderSdf(y0, y1, m_radius, false));
  BaseIF* insideCylinder = (BaseIF*)(new CylinderSdf(x0, x1, m_radius - m_curv, false));

  BaseIF* torusBottom = (BaseIF*)(new TorusSdf(y0, majorRadius, minorRadius, false));
  BaseIF* torusTop    = (BaseIF*)(new TorusSdf(y1, majorRadius, minorRadius, false));

  // Make the intersection of these
  Vector<BaseIF*> parts;
  parts.push_back(mainCylinder);
  parts.push_back(insideCylinder);
  parts.push_back(torusBottom);
  parts.push_back(torusTop);

  BaseIF* isectIF = (BaseIF*)(new IntersectionIF(parts));

  for (int i = 0; i < parts.size(); i++) {
    delete parts[i];
  }

  return isectIF;
}
#endif

#include <CD_NamespaceFooter.H>
