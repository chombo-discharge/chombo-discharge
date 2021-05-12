/*!
  @file rte_solver.cpp
  @brief Implementation of rte_solver.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "rte_solver.H"
#include "data_ops.H"

#include <ParmParse.H>

#include "CD_NamespaceHeader.H"

rte_solver::rte_solver(){
  m_name       = "rte_solver";
  m_class_name = "rte_solver";
}

rte_solver::~rte_solver(){
  
}

std::string rte_solver::get_name(){
  return m_name;
}

const std::string rte_solver::get_Realm() const {
  return m_Realm;
}

Vector<std::string> rte_solver::get_plotvar_names() const {
  CH_TIME("rte_solver::get_plotvar_names");
  if(m_verbosity > 5){
    pout() << m_name + "::get_plotvar_names" << endl;
  }
  
  Vector<std::string> names(0);
  
  if(m_plot_phi) names.push_back(m_name + " phi");
  if(m_plot_src) names.push_back(m_name + " source");

  return names;
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

  return converged;
}

bool rte_solver::advance(const Real a_dt, EBAMRCellData& a_state, const bool a_zerophi){
  CH_TIME("rte_solver::advance(dt, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(dt, state)" << endl;
  }

  const bool converged = this->advance(a_dt, a_state, m_source, a_zerophi);

  return converged;
}

void rte_solver::set_Realm(const std::string a_Realm){
  CH_TIME("rte_solver::set_Realm");
  if(m_verbosity > 5){
    pout() << m_name + "::set_Realm" << endl;
  }

  m_Realm = a_Realm;
}

void rte_solver::set_rte_species(const RefCountedPtr<rte_species> a_rte_species){
  CH_TIME("rte_solver::set_rte_species");
  if(m_verbosity > 5){
    pout() << m_name + "::set_rte_species" << endl;
  }

  m_rte_species = a_rte_species;
  m_name = m_rte_species->get_name();
}

void rte_solver::set_phase(const phase::which_phase a_phase){
  CH_TIME("rte_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void rte_solver::sanityCheck(){
  CH_TIME("rte_solver::sanityCheck");
  if(m_verbosity > 5){
    pout() << m_name + "::sanityCheck" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_rte_species.isNull());
  CH_assert(!m_ebis.isNull());
}

void rte_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("rte_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  
  m_compgeom = a_compgeom;

  const RefCountedPtr<mfis> mfis = m_compgeom->get_mfis();
  
  this->set_ebis(mfis->getEBIndexSpace(m_phase));
}

void rte_solver::set_ebis(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("rte_solver::set_ebis");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebis" << endl;
  }

  m_ebis = a_ebis;
}

void rte_solver::set_amr(const RefCountedPtr<AmrMesh>& a_amr){
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

  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
  

}

void rte_solver::set_source(const EBAMRCellData& a_source){
  CH_TIME("rte_solver::set_source(ebamrcell)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source(ebamrcell)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_source[lvl]->localCopyTo(*m_source[lvl]);
  }

  m_amr->averageDown(m_source, m_Realm, m_phase);
  m_amr->interpGhost(m_source, m_Realm, m_phase);
}

void rte_solver::set_source(const Real a_source){
  CH_TIME("rte_solver::set_source(constant)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source(constant)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (int comp = 0; comp < m_source[lvl]->nComp(); comp++){
      data_ops::set_value(*m_source[lvl], a_source, comp);
    }
  }

  m_amr->averageDown(m_source, m_Realm, m_phase);
  m_amr->interpGhost(m_source, m_Realm, m_phase);
}

void rte_solver::set_plot_variables(){
  CH_TIME("rte_solver::set_plot_variables");
  if(m_verbosity > 5){
    pout() << m_name + "::set_plot_variables" << endl;
  }

  m_plot_phi = false;
  m_plot_src = false;

  ParmParse pp("rte_solver");
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi") m_plot_phi = true;
    else if(str[i] == "src") m_plot_src = true;
  }
}

int rte_solver::get_num_plotvars() const {
  CH_TIME("rte_solver::get_num_plotvars");
  if(m_verbosity > 5){
    pout() << m_name + "::get_num_plotvars" << endl;
  }

  int num_output = 0;

  if(m_plot_phi) num_output = num_output + 1;
  if(m_plot_src) num_output = num_output + 1;

  return num_output;
}

void rte_solver::initial_data() {
  CH_TIME("rte_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }
  
  data_ops::set_value(m_state, 0.0);
}

void rte_solver::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("rte_solver::write_plot_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_data" << endl;
  }

  if(m_plot_phi) write_data(a_output, a_comp, m_state,  true);
  if(m_plot_src) write_data(a_output, a_comp, m_source, false);
}

void rte_solver::write_data(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp){
  CH_TIME("rte_solver::write_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_data" << endl;
  }

  const int comp = 0;
  const int ncomp = a_data[0]->nComp();

  const Interval src_interv(0, ncomp-1);
  const Interval dst_interv(a_comp, a_comp + ncomp - 1);

  // Copy data onto scratch
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_Realm, m_phase, ncomp);
  data_ops::copy(scratch, a_data);

  // Interp if we should
  if(a_interp){
    m_amr->interpToCentroids(scratch, m_Realm, phase::gas);
  }

  m_amr->averageDown(scratch, m_Realm, m_phase);
  m_amr->interpGhost(scratch, m_Realm, m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(a_output.get_Realm() == m_Realm){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else{
      scratch[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }

  data_ops::set_covered_value(a_output, a_comp, 0.0);

  a_comp += ncomp;
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

RefCountedPtr<rte_species>& rte_solver::get_species(){
  return m_rte_species;
}
#include "CD_NamespaceFooter.H"
