/*!
  @file rte_solver.cpp
  @brief Implementation of rte_solver.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "rte_solver.H"
#include "data_ops.H"

rte_solver::rte_solver(){

  this->set_phase();
}

rte_solver::~rte_solver(){
}

bool rte_solver::is_stationary(){
  return m_stationary;
}

bool rte_solver::advance(const Real a_dt, const bool a_zerophi){
  CH_TIME("rte_solver::advance(dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(dt)" << endl;
  }

  const bool converged = this->advance(a_dt, m_state, a_zerophi);

  m_time += a_dt;
  m_step++;

  return converged;
}

bool rte_solver::advance(const Real a_dt, EBAMRCellData& a_state, const bool a_zerophi){
  CH_TIME("rte_solver::advance(dt, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(dt, state)" << endl;
  }

  const bool converged = this->advance(a_dt, a_state, m_source, a_zerophi);

  m_time += a_dt;
  m_step++;

  return converged;
}

void rte_solver::set_photon_group(const RefCountedPtr<photon_group> a_photon_group){
  CH_TIME("rte_solver::set_photon_group");
  if(m_verbosity > 5){
    pout() << m_name + "::set_photon_group" << endl;
  }

  m_photon_group = a_photon_group;
  m_name = m_photon_group->get_name();
}

void rte_solver::set_phase(const phase::which_phase a_phase){
  CH_TIME("rte_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void rte_solver::sanity_check(){
  CH_TIME("rte_solver::sanity_check");
  if(m_verbosity > 5){
    pout() << m_name + "::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_physdom.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_photon_group.isNull());
  CH_assert(!m_ebis.isNull());
}

void rte_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("rte_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  
  m_compgeom = a_compgeom;

  const RefCountedPtr<mfis> mfis = m_compgeom->get_mfis();
  
  this->set_ebis(mfis->get_ebis(m_phase));
}

void rte_solver::set_ebis(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("rte_solver::set_ebis");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebis" << endl;
  }

  m_ebis = a_ebis;
}

void rte_solver::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("rte_solver::set_physical_domain");
  if(m_verbosity > 5){
    pout() << m_name + "::set_physical_domain" << endl;
  }

  m_physdom = a_physdom;
}

void rte_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("rte_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }

  m_amr = a_amr;
}

void rte_solver::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("rte_solver::set_time");
  if(m_verbosity > 5){
    pout() << m_name + "::set_time" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void rte_solver::set_stationary(const bool a_stationary) {
  CH_TIME("rte_solver::set_stationary");
  if(m_verbosity > 5){
    pout() << m_name + "::set_stationry" << endl;
  }
  
  m_stationary = a_stationary;
}

void rte_solver::set_verbosity(const int a_verbosity){
  CH_TIME("rte_solver::set_verbosity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
  
  m_verbosity = a_verbosity;
}

void rte_solver::set_source(const EBAMRCellData& a_source){
  CH_TIME("rte_solver::set_source(ebamrcell)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source(ebamrcell)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_source[lvl]->copyTo(*m_source[lvl]);
  }

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);
}

void rte_solver::set_source(const Real a_source){
  CH_TIME("rte_solver::set_source(constant)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source(constant)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (int comp = 0; comp < m_source[lvl]->nComp(); comp++){
      data_ops::set_value(*m_source[lvl], a_source, comp);
    }
  }

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);
}

void rte_solver::initial_data() {
  CH_TIME("rte_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }
  
  data_ops::set_value(m_state, 0.0);
}

Real rte_solver::get_time() const{
  CH_TIME("rte_solver::get_time");
  if(m_verbosity > 5){
    pout() << m_name + "::get_time" << endl;
  }
  
  return m_time;
}

phase::which_phase rte_solver::get_phase(){
  CH_TIME("rte_solver::get_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::get_phase" << endl;
  }

  return m_phase;
}

EBAMRCellData& rte_solver::get_state(){
  return m_state;
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


