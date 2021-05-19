/*!
  @file    profile_cylinder_if.cpp
  @brief   Implementation of profile_cylinder_if.H
  @author  Robert Marskar
*/

#include "profile_cylinder_if.H"
#include <CD_RoundedCylinderIF.H>
#include "cylinder_if.H"
#include "torus_if.H"

#include "TransformIF.H"
#include <SmoothUnion.H>

#include "CD_NamespaceHeader.H"

profile_cylinder_if::profile_cylinder_if(const RealVect  a_endPoint1,
					 const RealVect  a_endPoint2,
					 const Real      a_cylinderRadius,
					 const Real      a_torusMajorRadius,
					 const Real      a_torusMinorRadius,
					 const Real      a_ccDistance,
					 const Real      a_shift,
					 const Real      a_roundingRadius,
					 const int       a_numLeft,
					 const int       a_numRight,
					 const bool      a_fluidInside){
  if(SpaceDim != 3){
    MayDay::Abort("profile_cylinder_if::profile_cylinder_if - this is a 3D object!");
  }
  

  const Real length     = (a_endPoint2 - a_endPoint1).vectorLength();
  const RealVect zhat   = RealVect(D_DECL(0., 0., 1));
  const RealVect center = 0.5*length*zhat;

  Vector<BaseIF*> parts;

  BaseIF* cyl   = (BaseIF*) new RoundedCylinderIF(RealVect::Zero, length*zhat, a_cylinderRadius, a_roundingRadius, a_fluidInside);
  BaseIF* torus = (BaseIF*) new torus_if(center, a_torusMajorRadius, a_torusMinorRadius,  !a_fluidInside);

  // "Left" profiles
  for (int ileft = 0; ileft < a_numLeft; ileft++){
    TransformIF* transif = new TransformIF(*torus);

    const RealVect shift = (-ileft + 0.5)*a_ccDistance*zhat + a_shift*zhat;

    transif->translate(shift);
    parts.push_back(transif);
  }

  // "Right" profiles
  for (int iright = 0; iright < a_numRight; iright++){
    TransformIF* transif = new TransformIF(*torus);

    const RealVect shift = (iright + 0.5)*a_ccDistance*zhat + a_shift*zhat;

    transif->translate(shift);
    parts.push_back(transif);
  }

  parts.push_back(cyl);

  // Make the union
  BaseIF* base = (BaseIF*) (new SmoothUnion(parts, a_roundingRadius));

  // Translate it
  TransformIF* transif = new TransformIF(*base);

  if(a_endPoint2[2] >= a_endPoint1[2]){
    transif->rotate(zhat, a_endPoint2-a_endPoint1);
    transif->translate(a_endPoint1);
  }
  else{
    transif->rotate(zhat, a_endPoint1-a_endPoint2);
    transif->translate(a_endPoint2);
  }

  // Ok, we're done. Set m_baseif and clear up memory.
  m_baseif = RefCountedPtr<BaseIF> (transif);

  for (int i = 0; i < parts.size(); i++){
    delete parts[i];
  }
  delete base;
}

profile_cylinder_if::profile_cylinder_if(const profile_cylinder_if& a_inputIF){
  m_baseif = a_inputIF.m_baseif;
}
  
profile_cylinder_if::~profile_cylinder_if(){

}

Real profile_cylinder_if::value(const RealVect& a_pos) const{
  return m_baseif->value(a_pos);
}

BaseIF* profile_cylinder_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new profile_cylinder_if(*this));
}
#include "CD_NamespaceFooter.H"
