/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Rod.cpp
  @brief  Implementation CD_Rod.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Rod.H>
#include <CD_NamespaceHeader.H>

using Vec3    = EBGeometry::Vec3T<Real>;
using ImpFunc = EBGeometry::ImplicitFunction<Real>;

Rod::Rod() noexcept
{
  ParmParse pp("Rod");

  bool live;
  Real radius;

  Vector<Real> vec1(3);
  Vector<Real> vec2(3);

  pp.get("radius", radius);
  pp.get("live", live);
  pp.getarr("point1", vec1, 0, 3);
  pp.getarr("point2", vec2, 0, 3);

  Vec3 point1;
  Vec3 point2;
  point1 = Vec3(vec1[0], vec1[1], vec1[2]);
  point2 = Vec3(vec2[0], vec2[1], vec2[2]);

  std::shared_ptr<ImpFunc> rod = std::make_shared<EBGeometry::CapsuleSDF<Real>>(point1, point2, radius);

  const auto implicitFunction = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(rod, true));

  m_electrodes.push_back(Electrode(implicitFunction, live));
}

#include <CD_NamespaceFooter.H>
