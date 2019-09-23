/*!
  @file rod_if.cpp
  @brief Implementation of rod_if
  @date Nov. 2017
  @author Robert Marskar
*/

#include "rod_if.H"
#include "new_sphere_if.H"
#include "cylinder_if.H"

#include <SphereIF.H>
#include <IntersectionIF.H>
#include <PolyGeom.H>

rod_if::rod_if(const RealVect& a_center1,
	       const RealVect& a_center2,
	       const Real&     a_radius,
	       const bool&     a_inside){
  

  const RealVect axis    = (a_center2 - a_center1);
  const RealVect axisVec = axis/axis.vectorLength();

  // Find two new centers where can place cylinder edges and spheres
  const RealVect c1 = a_center1 + axisVec*a_radius;
  const RealVect c2 = a_center2 - axisVec*a_radius;

  // Build the cylinder
  Vector<BaseIF*> isects;
  isects.push_back(static_cast<BaseIF*> (new cylinder_if(c1, c2, a_radius, a_inside)));
#if 0
  isects.push_back(static_cast<BaseIF*> (new SphereIF(a_radius, c1, a_inside)));
  isects.push_back(static_cast<BaseIF*> (new SphereIF(a_radius, c2, a_inside)));
#else
  isects.push_back(static_cast<BaseIF*> (new new_sphere_if(c1, a_radius, a_inside)));
  isects.push_back(static_cast<BaseIF*> (new new_sphere_if(c2, a_radius, a_inside)));
#endif

  // Build the rod
  m_baseif = RefCountedPtr<BaseIF>(new IntersectionIF(isects));

  // Delete everything we allocated so far
  for (int i = 0; i < isects.size(); i++){
    delete isects[i];
  }
}

rod_if::rod_if(const rod_if& a_inputIF){
  this->m_baseif  = a_inputIF.m_baseif;
}

Real rod_if::value(const RealVect& a_point) const{
  return m_baseif->value(a_point);
}

BaseIF* rod_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new rod_if(*this));
}
