/*!
  @file rte_solver.cpp
  @brief Implementation of rte_solver.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "rte_solver.H"
#include "MFAliasFactory.H"


rte_solver::rte_solver(){
}

rte_solver::rte_solver(const RefCountedPtr<computational_geometry> a_compgeom){
  this->set_computational_geometry(a_compgeom);
}

rte_solver::~rte_solver(){
}

void rte_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  m_compgeom = a_compgeom;
}

void rte_solver::set_time(const Real a_time) {
  m_time = a_time;
}

Real rte_solver::get_time() const{
  return m_time;
}

EBAMRCellData& rte_solver::get_isotropic_state(){
  return m_isotropic_state;
}

EBAMRCellData& rte_solver::get_source(){
  return m_source;
}

EBAMRFluxData& rte_solver::get_kappa(){
  return m_kappa;
}

EBAMRIVData& rte_solver::get_kappa_eb(){
  return m_kappa_eb;
}


