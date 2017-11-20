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

void cdr_solver::advance(EBAMRCellData& a_state, const Real& a_dt){

  // Compute the right-hand side
  EBAMRCellData rhs;
  this->allocate(rhs);

  EBAMRDataOps::incr(a_state, rhs, a_dt);

  this->deallocate(rhs);
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const Real& a_dt){
  this->compute_rhs(a_rhs, m_state, a_dt);
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const EBAMRCellData& a_state, const Real& a_dt){

  EBAMRDataOps::setVal(a_rhs, 0.0);

  // Advective derivative
  EBAMRCellData advective_derivative;
  this->allocate(advective_derivative);
  this->compute_advective_derivative(advective_derivative, a_state);
  EBAMRDataOps::incr(a_rhs, advective_derivative, -1.0);
  this->deallocate(advective_derivative);

  // Diffusion term
  if(this->is_diffusive()){
    EBAMRCellData diffusion_term;
    this->allocate(diffusion_term);
    this->compute_diffusion_term(diffusion_term, a_state);
    EBAMRDataOps::incr(a_rhs, diffusion_term, 1.0);
  }

  // Source term
  EBAMRDataOps::incr(a_rhs, m_source, 1.0);
}

void cdr_solver::compute_advective_derivative(EBAMRCellData& a_adv_deriv, const EBAMRCellData& a_state){

  EBAMRDataOps::setVal(a_adv_deriv, 0.0); // Reset

  EBAMRFluxData face_state; 
  EBAMRIVData   noncons_div;
  EBAMRIVData   mass_diff;

  this->allocate(face_state);
  this->allocate(noncons_div);
  this->allocate(mass_diff);

  // Compute the advective derivative
  this->extrapolate_to_faces(face_state, a_state);                     // Do the face extrapolation
  this->conservative_divergence(a_adv_deriv, face_state, m_velo_face); // a_adv_deriv holds the conservative divergence
  this->nonconservative_divergence(noncons_div, a_adv_deriv);          // noncons_div holds the non-conservative divergence
  this->hybrid_divergence(a_adv_deriv, mass_diff, noncons_div);        // a_adv_deriv is now the hybrid divergence. 
  this->increment_flux_register(face_state, m_velo_face);              // Increment flux registers
  this->hyperbolic_redistribution(a_adv_deriv, mass_diff, a_state);    // Redistribute, use a_state as weights
  

  if(true){ // This should come from amr_mesh
    this->coarse_fine_increment(mass_diff);        // Increment the coarse-fine redistribution objects
    this->increment_redist_flux();                 // Increment flux registers with the redistribution stuff
    this->coarse_fine_redistribution(a_adv_deriv); // Redistribute
  }

  this->reflux(a_adv_deriv); // Reflux

  this->deallocate(face_state);
  this->deallocate(noncons_div);
  this->deallocate(mass_diff);
}

void cdr_solver::compute_diffusion_term(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state){

}

void cdr_solver::allocate(EBAMRCellData& a_data){

}

void cdr_solver::allocate(EBAMRFluxData& a_data){

}

void cdr_solver::allocate(EBAMRIVData& a_data){

}

void cdr_solver::deallocate(EBAMRCellData& a_data){

}

void cdr_solver::deallocate(EBAMRFluxData& a_data){

}

void cdr_solver::deallocate(EBAMRIVData& a_data){

}

void cdr_solver::conservative_divergence(EBAMRCellData&       a_cons_div,
					 const EBAMRFluxData& a_face_vel,
					 const EBAMRFluxData& a_face_state){

}

void cdr_solver::nonconservative_divergence(EBAMRIVData& a_noncons_div, const EBAMRCellData& a_cons_div){

}

void cdr_solver::hybrid_divergence(EBAMRCellData&     a_hybrid_div,
				   EBAMRIVData&       a_mass_diff,
				   const EBAMRIVData& a_noncons_div){

}

void cdr_solver::increment_flux_register(const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_velo_face){

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

Real cdr_solver::compute_dt(){

}

bool cdr_solver::is_diffusive(){
  return m_diffusive;
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
