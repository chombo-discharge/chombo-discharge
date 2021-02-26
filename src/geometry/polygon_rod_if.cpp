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
#include <SmoothIntersection.H>
#include <SmoothUnion.H>

#include "polygon_rod_if.H"
#include "cylinder_if.H"

polygon_rod_if::polygon_rod_if(const RealVect a_endPoint1,
			       const RealVect a_endPoint2,
			       const Real     a_radius,
			       const Real     a_cornerCurv,
			       const int      a_numSides,
			       const bool     a_fluidInside){

  if(SpaceDim != 3){
    MayDay::Abort("polygon_rod_if::polygon_rod_if - this is a 3D object!");
  }

  const Real length = (a_endPoint2 - a_endPoint1).vectorLength();
  const Real dTheta = 2.*M_PI/a_numSides;  // Internal angle
  const Real alpha  = M_PI/(a_numSides);   // Internal half angle 
  const Real beta   = 0.5*M_PI - alpha;    // External half angle
  
  const Real a = sin(alpha);             // Triangle factor, opposite internal angle
  const Real b = sin(beta);              // Triangle factor, opposite external angle
  const Real r = a_cornerCurv;           // Radius of curvature. Shortcut. 
  const Real R = a_radius - r + r/b;     // Radius accounts for curvature
  const Real c = R - r/b;                // Sphere center for rounding
  
  const RealVect zhat = RealVect(D_DECL(0., 0., 1.0));
  
  // Add the "sides" of the cylinder. 
  Vector<BaseIF*> planes;
  for (int iside = 0; iside < a_numSides; iside++){
    const Real theta = iside*dTheta;
    const RealVect n =   RealVect(D_DECL(cos(theta + 0.5*dTheta),sin(theta+0.5*dTheta), 0));
    const RealVect p = R*RealVect(D_DECL(cos(theta),             sin(theta),            0.));

    planes.push_back(static_cast<BaseIF*> (new PlaneIF(-n, p, a_fluidInside)));
  }

  // Add cuts above/below.
  planes.push_back((BaseIF*) (new PlaneIF(-zhat, length*zhat,    a_fluidInside)));
  planes.push_back((BaseIF*) (new PlaneIF( zhat, RealVect::Zero, a_fluidInside)));

  // Make a smooth union of those planes.
  BaseIF* isect = (BaseIF*) new SmoothUnion(planes, a_cornerCurv);

  // Do a transform, translating the rod into its specified place. 
  TransformIF* transif = new TransformIF(*isect);
  
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

  for (int i = 0; i < planes.size(); i++) {
    delete planes[i];
  }
}

polygon_rod_if::polygon_rod_if(const int&  a_nsides,
			       const Real& a_radius,
			       const Real& a_length,
			       const Real& a_curvrad,
			       const bool& a_inside){

  if(SpaceDim != 3){
    MayDay::Abort("polygon_rod_if::polygon_rod_if - this is a 3D object!");
  }
  Vector<BaseIF*> sides;
  Vector<BaseIF*> cut_sides;
  Vector<BaseIF*> corners;
  Vector<BaseIF*> parts;

  const Real dTheta = 2.*M_PI/a_nsides;  // Internal angle
  const Real alpha  = M_PI/(a_nsides);   // Internal half angle 
  const Real beta   = 0.5*M_PI - alpha;  // External half angle
  
  const Real a = sin(alpha);             // Triangle factor, opposite internal angle
  const Real b = sin(beta);              // Triangle factor, opposite external angle
  const Real r = a_curvrad;              // 
  const Real R = a_radius - r + r/b;     // Radius accounts for curvature
  const Real c = R - r/b;                // Sphere center for rounding

  
  // Add sides
  for (int iside = 0; iside < a_nsides; iside++){
    const Real theta = iside*dTheta;
    const RealVect n = RealVect(D_DECL(cos(theta + 0.5*dTheta),sin(theta+0.5*dTheta),0));
    const RealVect p = R*RealVect(D_DECL(cos(theta), sin(theta), 0.));

    sides.push_back(static_cast<BaseIF*> (new PlaneIF(-n, p, a_inside)));
  }

  // Cut above and below
  const RealVect norm  = BASISV(SpaceDim - 1);
  const RealVect point = 0.5*a_length*norm;
  sides.push_back(static_cast<BaseIF*> (new PlaneIF(-norm,  point, a_inside))); 
  sides.push_back(static_cast<BaseIF*> (new PlaneIF( norm, -point, a_inside))); 


  // Cut corners with planes
  for (int iside = 0; iside < a_nsides; iside++){
    const Real theta = iside*dTheta; // Plane angle
    const RealVect n = RealVect(D_DECL(cos(theta), sin(theta), 0.));
    const RealVect p = n*(c + b*r);
    cut_sides.push_back(static_cast<BaseIF*> (new PlaneIF(-n, p, a_inside)));
  }

  // Rounded corners
  for (int iside = 0; iside < a_nsides; iside++){
    const Real theta  = iside*dTheta;
    const RealVect n  = RealVect(D_DECL(cos(theta), sin(theta), 0.));
    const RealVect p  = c*n;
    const RealVect c1 = p - 0.5*a_length*RealVect(D_DECL(0.0, 0.0, 1.0));
    const RealVect c2 = p + 0.5*a_length*RealVect(D_DECL(0.0, 0.0, 1.0));

    corners.push_back(static_cast<BaseIF*> (new cylinder_if(c1, c2, r, a_inside)));
  }


  BaseIF* baserod = static_cast<BaseIF*> (new UnionIF(sides));      // Rod without rounded edges
  BaseIF* cutgeom = static_cast<BaseIF*> (new UnionIF(cut_sides));  // Planes that cut corners of rod
  
  // Cut rod at its corners
  parts.push_back(baserod);
  parts.push_back(cutgeom);
  BaseIF* cutrod = static_cast<BaseIF*> (new UnionIF(parts));

  // Add in cylinders on corners
  corners.push_back(cutrod);
  m_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(corners));

  // Clear up memory
  for (int i = 0; i < sides.size(); i++){
    delete sides[i];
  }
  for (int i = 0; i < parts.size(); i++){
    //    delete parts[i];
  }
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
