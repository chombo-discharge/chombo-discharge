/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RoundedBoxIF.cpp
  @brief  Implementation of CD_RoundedBoxIF.H
  @author Robert Marskar
*/

// Chombo includes
#include <PlaneIF.H>
#include <SmoothUnion.H>
#include <TransformIF.H>

// Our includes
#include <CD_RoundedBoxIF.H>
#include <CD_NamespaceHeader.H>

RoundedBoxIF::RoundedBoxIF(const RealVect a_loCorner,
                           const RealVect a_hiCorner,
                           const Real     a_curv,
                           const bool     a_fluidInside)
{
  m_fluidInside = a_fluidInside;

  const RealVect xyz = a_hiCorner - a_loCorner;

  // Make a slab whose center is at zero and whose widths are given by a_xyz
  Vector<BaseIF*> parts;
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const int      s = sign(sit());
      const RealVect n = s * BASISREALV(dir);
      const RealVect p = n * 0.5 * xyz[dir];

      BaseIF* baseif = (BaseIF*)new PlaneIF(n, p, true);
      parts.push_back(baseif);
    }
  }

  // Do rounded corners.
  BaseIF* bif = (BaseIF*)new SmoothUnion(parts, a_curv);

  // Translate box
  TransformIF* tif = new TransformIF(*bif);
  tif->translate(a_loCorner + 0.5 * xyz);

  delete bif;
  for (int i = 0; i < parts.size(); i++) {
    delete parts[i];
  }

  // Done;
  m_baseIF = RefCountedPtr<BaseIF>(tif);
}

RoundedBoxIF::RoundedBoxIF(const RoundedBoxIF& a_inputIF)
{
  m_fluidInside = a_inputIF.m_fluidInside;
  m_baseIF      = a_inputIF.m_baseIF;
}

RoundedBoxIF::~RoundedBoxIF() {}

Real
RoundedBoxIF::value(const RealVect& a_point) const
{

  // TLDR: m_baseIF designed so that f < 0 outside the box. This means that the fluid is outside. Revert if inside.
  Real retval = m_baseIF->value(a_point);

  if (m_fluidInside) {
    retval = -retval;
  }

  return retval;
}

BaseIF*
RoundedBoxIF::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new RoundedBoxIF(*this));
}

#include <CD_NamespaceFooter.H>
