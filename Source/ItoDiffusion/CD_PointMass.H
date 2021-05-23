/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PointMass.H
  @brief  Declaration of a simple class for storing a point with mass 
  @author Robert Marskar
*/

#ifndef CD_PointMass_H
#define CD_PointMass_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Declaration of a simple class for storing a point with mass. This is not a computational particle. 
*/
class PointMass{
public:
  PointMass();
  PointMass(const RealVect a_pos, const Real a_mass, const Real a_energy);
  PointMass(const std::vector<PointMass>& a_PointMasses);
  ~PointMass();

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
  bool canSplit() const {
    return m_mass*(1.0 + 1.E-6) >= 2.0;
  }
  
  RealVect m_pos;   // Position of this PointMass
  Real     m_mass;  // How many points does this this PointMass represent
  Real     m_energy;
};

#include <CD_NamespaceFooter.H>

#endif
