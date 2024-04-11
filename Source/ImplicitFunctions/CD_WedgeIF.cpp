/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file  CD_WedgeIF.cpp
  @brief Implementation of CD_WedgeIF.H
  @author Robert Marskar
*/

// Chombo includes
#include <IntersectionIF.H>
#include <SphereIF.H>
#include <PlaneIF.H>
#include <TransformIF.H>
#include <UnionIF.H>

// Our includes
#include <CD_WedgeIF.H>
#include <CD_CylinderSdf.H>
#include <CD_SphereSdf.H>
#include <CD_NamespaceHeader.H>

WedgeIF::WedgeIF(const int a_dir, const Real a_angle, const Real a_curv, const RealVect a_point, const bool a_inside)
{

  // By default, constructed with infinite extent along the z-axis
  const RealVect point = a_point - a_curv * RealVect(D_DECL(0, 1, 0));
  const Real     theta = 0.5 * a_angle * M_PI / 180.0; // Half opening angle
  const Real     cosT  = cos(theta);
  const Real     sinT  = sin(theta);

  const RealVect n1 = RealVect(D_DECL(-cosT, sinT, 0.0));
  const RealVect n2 = RealVect(D_DECL(cosT, sinT, 0.0));
  const RealVect p1 = point + n1 * a_curv;
  const RealVect p2 = point + n2 * a_curv;

  // Create the wedge planes as a union. Two of these form the sides, the last cuts off the top.
  Vector<BaseIF*> planes;
  planes.push_back(new PlaneIF(n1, p1, a_inside));
  planes.push_back(new PlaneIF(n2, p2, a_inside));
  BaseIF* base_planes = static_cast<BaseIF*>(new UnionIF(planes));
  BaseIF* cut_plane   = static_cast<BaseIF*>(new PlaneIF(RealVect(D_DECL(0.0, 1.0, 0.0)), p1, a_inside));

  // Union the planes and create the cylinder for rounding
  Vector<BaseIF*> plane_parts;
  plane_parts.push_back(base_planes);
  plane_parts.push_back(cut_plane);
  BaseIF* wedge = static_cast<BaseIF*>(new UnionIF(plane_parts));

#if CH_SPACEDIM == 3
  const RealVect c1        = point - 1.E10 * RealVect(BASISV(2));
  const RealVect c2        = point + 1.E10 * RealVect(BASISV(2));
  BaseIF*        round_cyl = static_cast<BaseIF*>(new CylinderSdf(c1, c2, a_curv, !a_inside));
#else
  BaseIF* round_cyl = static_cast<BaseIF*>(new SphereSdf(point, a_curv, !a_inside));
#endif

  // All parts
  Vector<BaseIF*> parts;
  parts.push_back(wedge);
  parts.push_back(round_cyl);

  // Rotate the wedge
  TransformIF* rot_wedge = new TransformIF(*(new IntersectionIF(parts)));
  if (a_dir == 0)
    rot_wedge->rotate(RealVect(BASISV(1)), RealVect(BASISV(0)), point);
  if (a_dir == 2)
    rot_wedge->rotate(RealVect(BASISV(1)), RealVect(BASISV(2)), point);

  m_baseIF = RefCountedPtr<BaseIF>(static_cast<BaseIF*>(rot_wedge));
}

WedgeIF::WedgeIF(const WedgeIF& a_inputIF)
{
  CH_assert(!a_inputIF.m_baseIF.isNull());
  m_baseIF = a_inputIF.m_baseIF;
}

WedgeIF::~WedgeIF()
{}

Real
WedgeIF::value(const RealVect& a_pos) const
{
  return m_baseIF->value(a_pos);
}

BaseIF*
WedgeIF::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new WedgeIF(*this));
}

#include <CD_NamespaceFooter.H>
