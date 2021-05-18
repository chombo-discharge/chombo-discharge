/*!
  @file   point_mass.H
  @brief  Declaration of a simple class for storing a point with mass 
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _POINTMASS_H_
#define _POINTMASS_H_

#include "RealVect.H"

#include "CD_NamespaceHeader.H"
  
/*!
  @brief  Declaration of a simple class for storing a point with mass 
*/
class point_mass{
public:
  point_mass();
  point_mass(const RealVect a_pos, const Real a_mass, const Real a_energy);
  point_mass(const std::vector<point_mass>& a_point_masses);
  ~point_mass();

  inline
  void define(const RealVect& a_pos, const Real& a_mass, const Real& a_energy){
    m_pos    = a_pos;
    m_mass   = a_mass;
    m_energy = a_energy;
  }
  
  inline
  Real& operator[](int a_dir) {
    return m_pos[a_dir];
  }

  inline
  const Real& operator[](int a_dir) const {
    return m_pos[a_dir];
  }

  inline
  RealVect& pos(){
    return m_pos;
  }
  
  inline
  const RealVect& pos() const{
    return m_pos;
  }

  inline
  Real& mass(){
    return m_mass;
  }
  
  inline
  const Real& mass() const{
    return m_mass;
  }

  inline
  Real& energy(){
    return m_energy;
  }
  
  inline
  const Real& energy() const{
    return m_energy;
  }
  
  inline
  bool can_split() const {
    return m_mass*(1.0 + 1.E-6) >= 2.0;
  }
  
  RealVect m_pos;   // Position of this point_mass
  Real     m_mass;  // How many points does this this point_mass represent
  Real     m_energy;
};
#include "CD_NamespaceFooter.H"

#endif
