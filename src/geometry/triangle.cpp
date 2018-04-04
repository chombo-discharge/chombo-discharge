/*!
  @brief  triangle.cpp
  @brief  Implementation of triangle.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "triangle.H"
#include "PolyGeom.H"

#define USE_3D_METHOD 1
#define DEV_VERB      1

triangle::triangle(){
#if CH_SPACEDIM==2
  MayDay::Abort("triangle::triangle - this is a purely 3D class");
#endif
}

triangle::triangle(const RealVect a_x1, const RealVect a_x2, const RealVect a_x3, const RealVect a_normal){
  this->define(a_x1, a_x2, a_x3, a_normal);
}

triangle::~triangle(){

}

void triangle::define(const RealVect a_x1, const RealVect a_x2, const RealVect a_x3, const RealVect a_normal){
  m_x1       = a_x1;
  m_x2       = a_x2;
  m_x3       = a_x3;
  m_normal   = a_normal;
  m_normal   = m_normal/m_normal.vectorLength();
  m_centroid = (m_x1 + m_x2 + m_x3)/3.0;

  const RealVect x21 = m_x2 - m_x1;
  const RealVect x31 = m_x3 - m_x1;
  const RealVect x32 = m_x3 - m_x2;

  m_v1 =  x21/x21.vectorLength() + x31/x31.vectorLength();
  m_v2 =  x32/x32.vectorLength() - x21/x21.vectorLength();
  m_v3 = -x31/x31.vectorLength() - x32/x32.vectorLength();

  // Create bounding sphere
}

Real triangle::signed_distance(const RealVect a_x0) const {
#if USE_3D_METHOD
  return signed_distance3D(a_x0);
#else
  return signed_distance2D(a_x0);
#endif
}

RealVect triangle::centroid() const{
  return m_centroid;
}

Real triangle::signed_distance3D(const RealVect a_x0) const{

  Real retval = 1.234567E89; // Return value
  
  // Vector from m_x1 to the projection of x0 onto the triangle plane
  const RealVect xp    = a_x0 - PolyGeom::dot(m_normal, a_x0 - m_x1)*m_normal;

  // Transient things for outside test
  const Real f1 = PolyGeom::dot(PolyGeom::cross(m_v1, m_x1-xp), m_normal);
  const Real f2 = PolyGeom::dot(PolyGeom::cross(m_v2, m_x2-xp), m_normal);
  const Real f3 = PolyGeom::dot(PolyGeom::cross(m_v3, m_x3-xp), m_normal);


  // Outside test
  bool outside;

  
  if(f1 > 0.0 && f2 < 0.0){ // 1-2 edge test
    if(PolyGeom::dot(PolyGeom::cross(m_x1-xp, m_x2-xp), m_normal) < 0.0){
      outside = true;
    }
  }
 else if(f1 < 0.0 && f2 > 0.0){ // 1-2 edge test
   if(PolyGeom::dot(PolyGeom::cross(m_x2-xp, m_x1-xp), m_normal) < 0.0){
      outside = true;
    }
  }
  else if(f1 > 0.0 && f3 < 0.0){ // 1-3 edge test
    if(PolyGeom::dot(PolyGeom::cross(m_x1-xp, m_x3-xp), m_normal) < 0.0){
      outside = true;
    }
  }
  else if(f1 < 0.0 && f3 > 0.0){ // 1-3 edge test
    if(PolyGeom::dot(PolyGeom::cross(m_x3-xp, m_x1-xp), m_normal) < 0.0){
      outside = true;
    }
  }
  else if(f2 > 0.0 && f3 < 0.0){ // 2-3 edge test
    if(PolyGeom::dot(PolyGeom::cross(m_x2-xp, m_x3-xp), m_normal) < 0.0){
      outside = true;
    }
  }
  else if(f2 < 0.0 && f3 > 0.0){ // 2-3 edge test
    if(PolyGeom::dot(PolyGeom::cross(m_x3-xp, m_x2-xp), m_normal) < 0.0){
      outside = true;
    }
  }
  else {
    outside = false;
  }

  
  if(outside){

    Real ret12;
    Real ret13;
    Real ret23;

    { // Find distance from outside point to edge/vertex on the x1-x2 edge
      const RealVect R   = PolyGeom::cross(PolyGeom::cross(m_x2-xp, m_x1-xp), m_x2-m_x1);
      const RealVect r   = R/R.vectorLength();
      const RealVect xpp = xp + PolyGeom::dot(m_x1-xp, r)*r;
      const Real t       = (xpp - m_x1).vectorLength()/(m_x2 - m_x1).vectorLength();

      if(t <= 0.0){ // Closest to x1
	ret12 = PolyGeom::dot(a_x0-m_x1, m_normal);
      }
      else if(t >= 1.0){ // Closest to x2
	ret12 = PolyGeom::dot(a_x0-m_x2, m_normal);
      }
      else{ // Lies on the line connecting x1 and x2
	MayDay::Abort("triangle::signed_distance - stop");
      }
    }

    // Shortest distance to edge/vertex
    retval = ret12;
    retval = Min(retval, ret13);
    retval = Min(retval, ret23);
  }
  else{
    retval = PolyGeom::dot(a_x0 - m_x1, m_normal);
  }


  
#if DEV_VERB
  pout() << "x0 = " << a_x0 << endl;
  pout() << "x1 = " << m_x1 << endl;
  pout() << "x2 = " << m_x2 << endl;
  pout() << "x3 = " << m_x3 << endl;
  pout() << "xp  = " << xp << endl;
  pout() << "f1 = " << f1 << endl;
  pout() << "f2 = " << f2 << endl;
  pout() << "f3 = " << f3 << endl;
  pout() << "outside = " << outside << endl;
  pout() << "dist  = " << retval << endl;
#endif
  
  return retval;
}

Real triangle::signed_distance2D(const RealVect a_x0) const{
  MayDay::Warning("triangle::signed_distance2D - not implemented");

  return 0.;
}
