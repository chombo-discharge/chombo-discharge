/*!
  @file polygon_rod_if.cpp
  @brief Implementation of polygon_rod_if.hpp
  @author Robert Marskar
  @date Aug. 2016
*/

#include <IntersectionIF.H>
#include <SphereIF.H>
#include <PlaneIF.H>
#include <TransformIF.H>
#include <UnionIF.H>

#include "polygon_rod_if.H"
#include "cylinder_if.H"

//
polygon_rod_if::polygon_rod_if(const int&  a_nsides,
			       const Real& a_radius,
			       const Real& a_length,
			       const Real& a_cornerCurv,
			       const bool& a_inside){

  // CSG things
  Vector<BaseIF*> planes;
  Vector<BaseIF*> parts;

  
  // Sides
  const Real dTheta = 2.*M_PI/a_nsides;
  for (int iside = 0; iside < a_nsides; iside++){
    const Real theta  = iside*dTheta;
    const RealVect n1 = RealVect(D_DECL(cos(theta + 0.5*dTheta),sin(theta+0.5*dTheta),0));
    const RealVect n2 = RealVect(D_DECL(cos(theta),sin(theta),0));
    const RealVect p1 = a_radius*RealVect(D_DECL(cos(theta), sin(theta), 0.));
    const RealVect p2 = (a_radius - a_cornerCurv)*RealVect(D_DECL(cos(theta), sin(theta), 0.));
    
    // Sides
    planes.push_back(static_cast<BaseIF*> (new PlaneIF(-n1, p1, a_inside)));

#if AMREX_SPACEDIM == 2

#elif AMREX_SPACEDIM == 3
    const RealVect v  = RealVect(0., 0., 0.5*a_length);
    const RealVect p  = (a_radius - a_cornerCurv)*RealVect(D_DECL(cos(theta), sin(theta), 0.));
    const RealVect c1 = p2 + v;
    const RealVect c2 = p2 - v;
    parts.push_back(static_cast<BaseIF*> (new cylinder_if(c1, c2, a_cornerCurv, a_inside)));
#endif
  }

  // Top/Bottom planes
#if AMREX_SPACEDIM == 3
  const RealVect norm  = BASISV(SpaceDim - 1);
  const RealVect point = 0.5*a_length*norm;
  planes.push_back(static_cast<BaseIF*> (new PlaneIF(-norm,  point, a_inside))); // Top
  planes.push_back(static_cast<BaseIF*> (new PlaneIF( norm, -point, a_inside))); // Bottom
#endif


  // Create rod without rounded edges. 
  BaseIF* unionIF = static_cast<BaseIF*> (new UnionIF(planes));

  // Intersection with cylinders to round edges
  parts.push_back(unionIF);
  m_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts));

  // Release memory
  for (int i = 0; i < planes.size(); i++){
    delete planes[i];
  }
  for (int i = 0; i < parts.size(); i++){
    //    delete parts[i];
  }
  delete unionIF;
}
    
polygon_rod_if::polygon_rod_if(const polygon_rod_if& a_inputIF){
  CH_assert(!a_inputIF.m_baseif.isNull());
  m_baseif = a_inputIF.m_baseif;
}

polygon_rod_if::~polygon_rod_if(){
}

Real polygon_rod_if::value(const RealVect& a_pos) const {
  return m_baseif->value(a_pos);
}

BaseIF* polygon_rod_if::newImplicitFunction() const {
  return static_cast<BaseIF*> (new polygon_rod_if(*this));
}
