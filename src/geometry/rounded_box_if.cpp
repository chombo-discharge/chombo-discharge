/*!
  @file rounded_box_if.cpp
  @brief Implementation of rounded_box_if.H
  @date Nov. 2017
  @author Robert Marskar
*/

#include "rounded_box_if.H"
#include "box_if.H"
#include "rod_if.H"
#include "new_sphere_if.H"

#include <IntersectionIF.H>
#include <UnionIF.H>
#include <SphereIF.H>

//
rounded_box_if::rounded_box_if(const RealVect& a_loCorner,
			       const RealVect& a_hiCorner,
			       const Real&     a_curv,
			       const bool&     a_inside){
 
  CH_assert(a_hiCorner > a_loCorner);
  CH_assert(a_hiCorner - a_loCorner > 2.0*a_curv*RealVect::Unit);
  
  // Remember the parameters
  m_loCorner = a_loCorner;
  m_hiCorner = a_hiCorner;
  m_curv     = a_curv;
  m_inside   = a_inside;

  // Build the box in 2D
#if CH_SPACEDIM == 2
  buildBox2D();
#elif CH_SPACEDIM == 3
  buildBox3D();
#endif
}

//
rounded_box_if::rounded_box_if(const rounded_box_if& a_inputIF){
   
  // Copy everything
  m_hiCorner   = a_inputIF.m_hiCorner;
  m_loCorner   = a_inputIF.m_loCorner;
  m_curv       = a_inputIF.m_curv;
  m_inside     = a_inputIF.m_inside;
  m_baseif     = a_inputIF.m_baseif;
}

#if CH_SPACEDIM == 2
void rounded_box_if::buildBox2D(){

  const RealVect base    = m_curv*RealVect::Unit;
  const RealVect loLeft  = m_loCorner;
  const RealVect hiRight = m_hiCorner;
  const RealVect hiLeft  = RealVect(m_loCorner[0], m_hiCorner[1]);
  const RealVect loRight = RealVect(m_hiCorner[0], m_loCorner[1]);

  // We will subtract four boxes from each of the corners of a parent box by taking the union of their outsides
  Vector<BaseIF*> unions;
  unions.push_back(static_cast<BaseIF*>  (new box_if(m_loCorner,     m_hiCorner,      m_inside)));
  unions.push_back(static_cast<BaseIF*>  (new box_if(loLeft  - base, loLeft  + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*>  (new box_if(hiRight - base, hiRight + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*>  (new box_if(hiLeft  - base, hiLeft  + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*>  (new box_if(loRight - base, loRight + base, !m_inside)));
  BaseIF* unionIF = static_cast<BaseIF*> (new UnionIF(unions));

  
  // // Then, we will add in sphere at every corner
  Vector<BaseIF*> isects;
  const RealVect c1 = m_curv*RealVect( 1, 1);
  const RealVect c2 = m_curv*RealVect(-1,-1);
  const RealVect c3 = m_curv*RealVect( 1,-1);
  const RealVect c4 = m_curv*RealVect(-1, 1);
  isects.push_back(unionIF);
  isects.push_back(static_cast<BaseIF*> (new new_sphere_if(loLeft  + c1, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new new_sphere_if(hiRight + c2, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new new_sphere_if(hiLeft  + c3, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new new_sphere_if(loRight + c4, m_curv, m_inside)));

  // Now build our box
  m_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(isects));
  
  // Delete the stuff where we built the union
  for (int i = 0; i < unions.size(); i++){
    delete unions[i];
  }
  for (int i = 0; i < isects.size(); i++){
    delete isects[i];
  }
}
#endif

//
#if CH_SPACEDIM == 3
void rounded_box_if::buildBox3D(){

  // TLDR; A box in 3D has 12 edges; we build the 3D box by subtracting a base box on each edge, and then replace it by
  //       a cylinder. The 8 corners are then filled with spheres.
  //

  // This is the thickness of the boxes we subtract
  const RealVect base    = m_curv*RealVect::Unit;

  // Corners of the box
  const RealVect c1 = RealVect(m_loCorner[0], m_loCorner[1], m_loCorner[2]);
  const RealVect c2 = RealVect(m_loCorner[0], m_hiCorner[1], m_loCorner[2]);
  const RealVect c3 = RealVect(m_hiCorner[1], m_loCorner[1], m_loCorner[2]);
  const RealVect c4 = RealVect(m_hiCorner[1], m_hiCorner[1], m_loCorner[2]);
  const RealVect c5 = RealVect(m_loCorner[0], m_loCorner[1], m_hiCorner[2]);
  const RealVect c6 = RealVect(m_loCorner[0], m_hiCorner[1], m_hiCorner[2]);
  const RealVect c7 = RealVect(m_hiCorner[1], m_loCorner[1], m_hiCorner[2]);
  const RealVect c8 = RealVect(m_hiCorner[1], m_hiCorner[1], m_hiCorner[2]);


  // We will now subtract a base box from each edge
  Vector<BaseIF*> unions;
  unions.push_back(static_cast<BaseIF*> (new box_if(m_loCorner, m_hiCorner, m_inside))); // Base box
  unions.push_back(static_cast<BaseIF*> (new box_if(c1 - base, c2 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c1 - base, c5 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c1 - base, c3 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c2 - base, c4 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c2 - base, c6 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c3 - base, c4 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c3 - base, c7 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c4 - base, c8 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c5 - base, c6 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c5 - base, c7 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c6 - base, c8 + base, !m_inside)));
  unions.push_back(static_cast<BaseIF*> (new box_if(c7 - base, c8 + base, !m_inside)));
  BaseIF* unionIF = static_cast<BaseIF*> (new UnionIF(unions));

  // Add in a rod_if on each edge
  Vector<BaseIF*> isects;
  isects.push_back(unionIF);
  const RealVect c12 = m_curv*RealVect( 1, 0, 1);
  const RealVect c13 = m_curv*RealVect( 0, 1, 1);
  const RealVect c15 = m_curv*RealVect( 1, 1, 0);
  const RealVect c24 = m_curv*RealVect( 0,-1, 1);
  const RealVect c26 = m_curv*RealVect( 1,-1, 0);
  const RealVect c34 = m_curv*RealVect(-1, 0, 1);
  const RealVect c37 = m_curv*RealVect(-1, 1, 0);
  const RealVect c48 = m_curv*RealVect(-1,-1, 0);
  const RealVect c56 = m_curv*RealVect( 1, 0,-1);
  const RealVect c57 = m_curv*RealVect( 0, 1,-1);
  const RealVect c68 = m_curv*RealVect( 0,-1,-1);
  const RealVect c78 = m_curv*RealVect(-1, 0,-1);
  isects.push_back(static_cast<BaseIF*> (new rod_if(c1 + c12, c2 + c12, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c1 + c13, c3 + c13, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c1 + c15, c5 + c15, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c2 + c24, c4 + c24, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c2 + c26, c6 + c26, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c3 + c34, c4 + c34, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c3 + c37, c7 + c37, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c4 + c48, c8 + c48, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c5 + c56, c6 + c56, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c5 + c57, c7 + c57, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c6 + c68, c8 + c68, m_curv, m_inside)));
  isects.push_back(static_cast<BaseIF*> (new rod_if(c7 + c78, c8 + c78, m_curv, m_inside)));

  //
  m_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(isects));

  // Delete all the stuff we just created
  for (int i = 0; i < unions.size(); i++){
    delete unions[i];
  }
  for (int i = 0; i < isects.size(); i++){
    delete isects[i];
  }
}
#endif

//
rounded_box_if::~rounded_box_if(){
}

//
Real rounded_box_if::value(const RealVect& a_point) const{
  return m_baseif->value(a_point);
}

//
BaseIF* rounded_box_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new rounded_box_if(*this));
}

