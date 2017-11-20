/*!
  @file cdr_solver.cpp
  @brief Implementation of cdr_solver.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "cdr_solver.H"

#include <EBAMRDataOps.H>

cdr_solver::cdr_solver(){
  
}

cdr_solver::cdr_solver(const RefCountedPtr<computational_geometry>& a_compgeom){
  this->set_computational_geometry(a_compgeom);
}

cdr_solver::~cdr_solver(){

}

void cdr_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  m_compgeom = a_compgeom;
}

void cdr_solver::advance(const Real& a_dt){
  this->advance(m_state, a_dt);
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const Real& a_dt){
  this->compute_rhs(a_rhs, m_state, a_dt);
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const EBAMRCellData& a_state, const Real& a_dt){

  // Allocate temporary storage
  this->allocate_temporaries();

  // Delete temporary storage
  this->delete_temporaries();
}

void cdr_solver::allocate_temporaries(){

}

void cdr_solver::delete_temporaries(){

}

void cdr_solver::conservative_divergence(EBAMRCellData&       a_cons_div,
					 const EBAMRFluxData& a_face_vel,
					 const EBAMRFluxData& a_face_state){

}

void cdr_solver::nonconservative_divergence(EBAMRCellData& a_noncons_div, const EBAMRCellData& a_cons_div){

}

											    

void cdr_solver::hybrid_divergence(EBAMRCellData&       a_hybrid_div,
				   EBAMRIVData&         a_mass_diff,
				   const EBAMRCellData& a_cons_div,
				   const EBAMRCellData& a_noncons_div){

}


void cdr_solver::increment_flux_registers(const EBAMRCellData& a_face_state, const EBAMRCellData& a_velo_face){

}


void cdr_solver::coarse_fine_increment(const EBAMRIVData& m_mass_diff){

}


void cdr_solver::hyperbolic_redistribution(EBAMRCellData&       a_del_vel_rho,
					   const EBAMRIVData&   a_mass_diff,
					   const EBAMRCellData& a_redist_weights){

}


void cdr_solver::reflux_redist_dance(EBAMRCellData& a_del_vel_rho, const EBAMRCellData& a_redist_weights){

}


void cdr_solver::set_weights(const EBAMRCellData& a_redist_weights){

}


void cdr_solver::increment_redist_flux(){

}

void cdr_solver::coarse_fine_redistribution(EBAMRCellData& a_state){

}

void cdr_solver::reflux(EBAMRCellData& a_state){

}

EBAMRCellData& cdr_solver::get_state(){
  return m_state;
}

EBAMRCellData& cdr_solver::get_source(){
  return m_source;
}

EBAMRFluxData& cdr_solver::get_velo_face(){
  return m_velo_face;
}

EBAMRIVData& cdr_solver::get_velo_eb(){
  return m_velo_eb;
}

EBAMRFluxData& cdr_solver::get_diffco_face(){
  return m_diffco_face;
}

EBAMRIVData& cdr_solver::get_diffco_eb(){
  return m_diffco_eb;
}
