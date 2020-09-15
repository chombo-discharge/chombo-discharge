#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphereArrayIF.H"

#include "NamespaceHeader.H"

SphereArrayIF::SphereArrayIF(const Real&     a_radius,
                             const RealVect& a_firstCenter,
                             const RealVect& a_spacing)
{
  // Remember the parameters
  m_radius = a_radius;
  m_firstCenter = a_firstCenter;
  m_spacing = a_spacing;

  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}

SphereArrayIF::SphereArrayIF(const SphereArrayIF& a_inputIF)
{
  // Remember the parameters
  m_radius      = a_inputIF.m_radius;
  CH_assert(m_radius > 0);
  m_firstCenter = a_inputIF.m_firstCenter;
  m_spacing     = a_inputIF.m_spacing;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_spacing[idir] < 2.25*m_radius)
        {
          MayDay::Abort("SphereArrayIF: spheres too close");
        }
    }
  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}

SphereArrayIF::~SphereArrayIF()
{
}

Real SphereArrayIF::value(const RealVect& a_point) const
{
  //compute which center we are closest to
  RealVect dist = a_point;
  dist -= m_firstCenter;
  dist /= m_spacing;
  //make distance integer numbers in all directions
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      dist[idir] = rint(dist[idir]);
    }
  //now multiply the distance by the spacing and add it to the firstCenter;
  dist *= m_spacing;
  RealVect center = m_firstCenter;
  center += dist;

  // The distance squared for m_center to a_point
  Real distance2;

  // Compute the distance squared
  distance2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    Real cur;
    cur = a_point[idir] - center[idir];

    distance2 += cur*cur;
  }

  // Return the difference between the sqaures (zero on the sphere)
  Real retval = m_radius2 - distance2;

  return retval;
}
                                                       //
BaseIF* SphereArrayIF::newImplicitFunction() const
{
  SphereArrayIF* spherePtr = new SphereArrayIF(m_radius,
                                               m_firstCenter,
                                               m_spacing);

  return static_cast<BaseIF*>(spherePtr);
}


#include "NamespaceFooter.H"
