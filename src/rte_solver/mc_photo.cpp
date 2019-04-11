/*!
  @file   mc_photo.cpp
  @brief  Implementation of mc_photo.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "mc_photo.H"
#include "data_ops.H"

mc_photo::mc_photo(){

}

mc_photo::~mc_photo(){

}

bool mc_photo::advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi){
  data_ops::set_value(a_state, 1.2345);
}
  
void mc_photo::allocate_internals(){
  const int ncomp = 1;
  m_amr->allocate(m_state,  m_phase, ncomp); // This is the deposited 
  m_amr->allocate(m_source, m_phase, ncomp);
}
  
void mc_photo::cache_state(){

}

void mc_photo::deallocate_internals(){
  m_amr->deallocate(m_state);
  m_amr->deallocate(m_source);
}

void mc_photo::regrid(const int a_old_finest_level, const int a_new_finest_level){

}

void mc_photo::compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state){
  data_ops::set_value(a_ebflux, 0.0);
}


void mc_photo::compute_domain_flux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_state){
  data_ops::set_value(a_domainflux, 0.0);
}

void mc_photo::compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state){
  MayDay::Abort("mc_photo::compute_flux - I don't think this should ever be called.");
}

void mc_photo::compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state){
  MayDay::Abort("mc_photo::compute_density - I don't think this should ever be called.");
}


void mc_photo::write_plot_file(){

}

int mc_photo::query_ghost() const {
  return 3;
}
