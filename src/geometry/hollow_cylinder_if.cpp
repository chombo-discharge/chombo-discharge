/*!
  @file hollow_cylinder_if.cpp
  @brief Implementation of hollow_cylinder_if.H
  @date Jan. 2017
  @author Sigurd Midttun
*/

#include "cylinder_if.H"
#include "hollow_cylinder_if.H"
#include "rounded_cylinder_if.H"

#include <SmoothUnion.H>
#include <UnionIF.H>
#include <IntersectionIF.H>

namespace ChomboDischarge {

  hollow_cylinder_if::hollow_cylinder_if(const RealVect a_center1,
					 const RealVect a_center2,
					 const Real     a_majorRadius,
					 const Real     a_minorRadius,
					 const Real     a_curv,
					 const bool     a_fluidInside){

    // Make the SmoothUnion of this stuff. 
    Vector<BaseIF*> parts;

    RealVect axis = (a_center2 - a_center1);
    axis = axis/axis.vectorLength();

    const RealVect c2 = a_center2 + a_curv*axis;
    const RealVect c1 = a_center1 - a_curv*axis;
  
    BaseIF* bigCylinder   = (BaseIF*) (new rounded_cylinder_if(a_center1, a_center2, a_majorRadius, a_curv,  a_fluidInside));
    BaseIF* smallCylinder = (BaseIF*) (new cylinder_if        (a_center1, a_center2, a_minorRadius,         !a_fluidInside));

    parts.push_back(bigCylinder);
    parts.push_back(smallCylinder);
  
    // Make union
    m_baseif = RefCountedPtr<BaseIF> (new SmoothUnion(parts, 2*a_curv));
  }

  hollow_cylinder_if::hollow_cylinder_if(const hollow_cylinder_if& a_inputIF){
    m_baseif = a_inputIF.m_baseif;
  }

  Real hollow_cylinder_if::value(const RealVect& a_point) const{
    return m_baseif->value(a_point);
  }

  BaseIF* hollow_cylinder_if::newImplicitFunction() const{
    return (BaseIF*) (new hollow_cylinder_if(*this));
  }
}
