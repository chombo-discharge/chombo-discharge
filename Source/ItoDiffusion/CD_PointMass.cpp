/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PointMass.cpp
  @brief  Implementation of CD_PointMass.H
  @author Robert Marskar
*/

// Our includes
#include <CD_PointMass.H>
#include <CD_NamespaceHeader.H>
  
PointMass::PointMass(){

}

PointMass::PointMass(const RealVect a_pos, const Real a_mass, const Real a_energy){
  m_pos    = a_pos;
  m_mass   = a_mass;
  m_energy = a_energy;
}

PointMass::PointMass(const std::vector<PointMass>& a_PointMasses){
  m_pos    = RealVect::Zero;
  m_mass   = 0.0;
  m_energy = 0.0;

  for (int i = 0; i < a_PointMasses.size(); i++){
    const RealVect& p = a_PointMasses[i].m_pos;
    const Real& m     = a_PointMasses[i].m_mass;
    const Real& e     = a_PointMasses[i].m_energy;
    
    m_mass   += m;
    m_pos    += m*p;
    m_energy += m*e;
  }

  m_pos    = m_pos/m_mass;
  m_energy = m_energy/m_mass;
}

PointMass::~PointMass(){

}

#include <CD_NamespaceFooter.H>
