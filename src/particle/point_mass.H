/*!
  @file   point_mass.H
  @brief  Declaration of a simple class for storing a point with mass 
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _POINTMASS_H_
#define _POINTMASS_H_

#include "RealVect.H"

/*!
  @brief  Declaration of a simple class for storing a point with mass 
*/
class point_mass{
public:
  point_mass();
  point_mass(const RealVect a_pos, const Real a_mass);
  point_mass(const std::vector<point_mass>& a_point_masses);
  ~point_mass();

  inline Real operator[](int a_dir) const{
    return m_pos[a_dir];
  }
  
  inline const RealVect pos() const{
    return m_pos;
  }
  
  inline const Real mass() const{
    return m_mass;
  }
  
  inline bool can_split() const{
    return m_mass*(1.0 + 1.E-6) >= 2.0;
  }
  
  RealVect m_pos;   // Position of this point_mass
  Real     m_mass;  // How many points does this this point_mass represent
};

#endif
