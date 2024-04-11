/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PerlinSlabSdf.cpp
  @brief  Implementation of CD_PerlinSlabSdf.H
  @author Robert Marskar
*/

// Chombo includes
#include <SmoothUnion.H>
#include <PlaneIF.H>
#include <TransformIF.H>

// Our includes
#include <CD_PerlinSlabSdf.H>
#include <CD_PerlinPlaneSdf.H>
#include <CD_NamespaceHeader.H>

PerlinSlabSdf::PerlinSlabSdf(const RealVect a_ccPoint,
                             const RealVect a_normal,
                             const RealVect a_xyz,
                             const RealVect a_noiseFreq,
                             const int      a_octaves,
                             const Real     a_noiseAmp,
                             const Real     a_persistence,
                             const Real     a_cornerCurv,
                             const bool     a_reseed,
                             const bool     a_fluidInside)
{
  m_fluidInside = a_fluidInside;

  constexpr int up = CH_SPACEDIM - 1;

  // Make a slab whose center is at zero and whose widths are given by a_xyz
  Vector<BaseIF*> parts;
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const int      s = sign(sit());
      const RealVect n = s * BASISREALV(dir);
      const RealVect p = n * 0.5 * a_xyz[dir];

      if (dir == up && sit() == Side::Hi) {
        BaseIF*
          baseif = (BaseIF*)new PerlinPlaneSdf(n, p, true, a_noiseAmp, a_noiseFreq, a_persistence, a_octaves, a_reseed);
        parts.push_back(baseif);
      }
      else {
        BaseIF* baseif = (BaseIF*)new PlaneIF(n, p, true);
        parts.push_back(baseif);
      }
    }
  }

  // Do rounded corners.
  BaseIF* bif = (BaseIF*)new SmoothUnion(parts, a_cornerCurv);

  // // Rotate and translate into place
  TransformIF* tif = new TransformIF(*bif);
  tif->translate(-0.5 * BASISREALV(up) *
                 a_xyz[up]);             // Move so that "top" point is at RealVect::Zero. This is the rotation point
  tif->rotate(BASISREALV(up), a_normal); // Rotate so that +z/+y points along a_normal
  tif->translate(a_ccPoint);             // Translate to point

  // Done. Delete the rest.
  m_baseif = RefCountedPtr<BaseIF>(tif);

  for (int i = 0; i < parts.size(); i++) {
    delete parts[i];
  }
  delete bif;
}

PerlinSlabSdf::PerlinSlabSdf(const PerlinSlabSdf& a_inputIF)
{
  m_baseif      = a_inputIF.m_baseif;
  m_fluidInside = a_inputIF.m_fluidInside;
}

PerlinSlabSdf::~PerlinSlabSdf()
{}

Real
PerlinSlabSdf::value(const RealVect& a_pos) const
{
  Real retval = m_baseif->value(a_pos);

  if (m_fluidInside) {
    retval = -retval;
  }

  return retval;
}

BaseIF*
PerlinSlabSdf::newImplicitFunction() const
{
  return (BaseIF*)new PerlinSlabSdf(*this);
}

#include <CD_NamespaceFooter.H>
