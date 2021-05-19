/*!
  @file   air2_species.cpp
  @brief  Implementation of air2_species.H
  @author Robert Marskar
  @date   Feb. 2018
  @todo   Could definitely see ways to cut down on the typic 
*/

#include "air2_species.H"
#include "units.H" 

#include <ParmParse.H>
#include <PolyGeom.H>

air2::electron::electron(){
  m_name      = "electron density";
  m_unit      = "m-3";
  m_chargeNumber    = -1;
  m_isDiffusive = true;
  m_isMobile    = true;

  ParmParse pp("air2");
  Vector<Real> vec(SpaceDim);
  pp.get("initial_ionization", m_initial_ionization);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius", m_seed_radius);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_position = RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air2::electron::~electron(){

}

air2::positive_ion::positive_ion() {
  m_name      = "positive_ion";
  m_unit      = "m-3";
  m_chargeNumber    = 1;
  m_isDiffusive = false;
  m_isMobile    = false;

  ParmParse pp("air2");
  Vector<Real> vec(SpaceDim);
  pp.get("initial_ionization", m_initial_ionization);
  pp.get("seed_density",       m_seed_density);
  pp.get("seed_radius",        m_seed_radius);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_position = RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air2::positive_ion::~positive_ion(){

}

air2::Photon_one::Photon_one(){
  m_name   = "Photon_one";
  m_A      = 1.12E-4; 
  m_lambda = 4.15E-2;
  m_pO2    = 0.2*units::s_atm2pascal;
  
  ParmParse pp("air2");
  pp.query("Photon1_A_coeff",      m_A);
  pp.query("Photon1_lambda_coeff", m_lambda);
}

air2::Photon_one::~Photon_one(){

}

air2::Photon_two::Photon_two(){
  m_name   = "Photon_two";
  m_A      = 2.88E-3;
  m_lambda = 1.09E-1;
  m_pO2    = 0.2*units::s_atm2pascal;
  
  ParmParse pp("air2");
  pp.query("Photon2_A_coeff",      m_A);
  pp.query("Photon2_lambda_coeff", m_lambda);
}

air2::Photon_two::~Photon_two(){

}

air2::Photon_three::Photon_three(){
  m_name   = "Photon_three";
  m_A      = 2.76E-1;
  m_lambda = 6.69E-1;
  m_pO2    = 0.2*units::s_atm2pascal;
  
  ParmParse pp("air2");
  pp.query("Photon3_A_coeff",      m_A);
  pp.query("Photon3_lambda_coeff", m_lambda);
}

air2::Photon_three::~Photon_three(){

}

Real air2::electron::initialData(const RealVect a_pos, const Real a_time) const {
  const RealVect p = (a_pos - m_seed_position)/(m_seed_radius);
  return m_initial_ionization + m_seed_density*exp(-PolyGeom::dot(p, p));
}

Real air2::positive_ion::initialData(const RealVect a_pos, const Real a_time) const {
  const RealVect p = (a_pos - m_seed_position)/(m_seed_radius);
  return m_initial_ionization + m_seed_density*exp(-PolyGeom::dot(p, p));
}

Real air2::Photon_one::getKappa(const RealVect a_pos) const { return m_lambda*m_pO2/(sqrt(3.0));}
Real air2::Photon_one::get_lambda() const { return m_lambda;}
Real air2::Photon_one::get_A() const {return m_A;}
Real air2::Photon_one::get_pO2() const {return m_pO2;}

Real air2::Photon_two::getKappa(const RealVect a_pos) const { return m_lambda*m_pO2/(sqrt(3.0));}
Real air2::Photon_two::get_lambda() const { return m_lambda;}
Real air2::Photon_two::get_A() const {return m_A;}
Real air2::Photon_two::get_pO2() const {return m_pO2;}

Real air2::Photon_three::getKappa(const RealVect a_pos) const { return m_lambda*m_pO2/(sqrt(3.0));}
Real air2::Photon_three::get_lambda() const { return m_lambda;}
Real air2::Photon_three::get_A() const {return m_A;}
#include "CD_NamespaceFooter.H"
