/*!
  @file   real_box.cpp
  @brief  Implementation of real_box.H
  @author Robert Marskar
  @date   June 2018
*/

#include "real_box.H"

real_box::real_box(){
  m_lo = RealVect::Zero;
  m_hi = RealVect::Zero;
}


real_box::real_box(const RealVect a_lo, const RealVect a_hi){
  m_lo = a_lo;
  m_hi = a_hi;
}

real_box::real_box(const Box a_box, const RealVect a_origin, const Real a_dx){

  const IntVect lo = a_box.smallEnd();
  const IntVect hi = a_box.bigEnd();

  m_lo = a_origin + a_dx*RealVect(lo);
  m_hi = a_origin + a_dx*RealVect(hi);
}

real_box::~real_box(){

}

RealVect real_box::get_lo() const {
  return m_lo;
}

RealVect real_box::get_hi() const {
  return m_hi;
}

Vector<RealVect> real_box::get_corners() const {

  Vector<RealVect> corners(pow(2, SpaceDim));

  const RealVect DD = m_hi - m_lo;
  const RealVect DX = RealVect(D_DECL(DD[0], 0.0,   0.0));
  const RealVect DY = RealVect(D_DECL(0.0 ,  DD[1], 0.0));
#if CH_SPACEDIM==3
  const RealVect DZ = RealVect(D_DECL(0.0 ,  0.0,   DD[2]));
#endif

  corners[0] = m_lo;
  corners[1] = m_lo + DX;
  corners[2] = m_lo + DY;
  corners[3] = m_lo + DX + DY; // = m_hi
#if CH_SPACEDIM==3
  corners[4] = m_lo + DZ;
  corners[5] = m_lo + DZ + DX;
  corners[6] = m_lo + DZ + DY;
  corners[7] = m_lo + DX + DX + DY; // = m_hi
#endif

  return corners;
}

bool real_box::intersect(const real_box& a_box) const {

  //  bool ret = false;

  const RealVect LO = a_box.get_lo();
  const RealVect HI = a_box.get_hi();

  if (   LO[0] < m_hi[0] && HI[0] > m_lo[0]  //  Input x-edges either to left or righ
      && LO[1] < m_hi[1] && HI[1] > m_lo[1]
#if CH_SPACEDIM==3
      && LO[2] < m_hi[2] && HI[2] > m_lo[2]
#endif
	 ){
    return true;
  }
  else{
    return false;
  }
      
    
  // const Vector<RealVect> corners = a_box.get_corners();
  
  // for (int i = 0; i < corners.size(); i++){
  //   if(is_point_inside(corners[i])){
  //     ret = true;
  //   }
  // }

  // return ret;
}


bool real_box::is_point_inside(const RealVect a_point) const {

  bool ret = false;

  if(   a_point[0] >= m_lo[0] && a_point[0] <= m_hi[0] 
     && a_point[1] >= m_lo[1] && a_point[1] <= m_hi[1]
#if CH_SPACEDIM==3
     && a_point[2] >= m_lo[2] && a_point[2] <= m_hi[2]
#endif
     ){
    ret = true;
  }

  return ret;
}

bool real_box::is_box_inside(const real_box& a_box) const {

  Vector<RealVect> corners = a_box.get_corners();

  bool inside = true;
  return true;
  for (int i = 0; i < corners.size(); i++){

    // Check if any of the corners of the input box lies inside this box. If one of the
    // corners of the input box is outside, a_box is NOT completely contained by this box.
    const bool is_corner_inside = this->is_point_inside(corners[i]);

    if(!is_corner_inside){
      inside = false;
      break;
    }
  }

  return inside;
}
