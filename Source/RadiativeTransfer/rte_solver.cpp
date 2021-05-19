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
  m_className = "rte_solver";
}

rte_solver::~rte_solver(){
  
}

std::string rte_solver::getName(){
  return m_name;
}

const std::string rte_solver::getRealm() const {
  return m_realm;
}

Vector<std::string> rte_solver::getPlotVariableNames() const {
  CH_TIME("rte_solver::getPlotVariableNames");
  if(m_verbosity > 5){
    pout() << m_name + "::getPlotVariableNames" << endl;
  }
  
  Vector<std::string> names(0);
  
  if(m_plotPhi) names.push_back(m_name + " phi");
  if(m_plotSource) names.push_back(m_name + " source");

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

  const bool converged = this->advance(a_dt, m_phi, a_zerophi);

  return converged;
}

bool rte_solver::advance(const Real a_dt, EBAMRCellData& a_phi, const bool a_zerophi){
  CH_TIME("rte_solver::advance(dt, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(dt, state)" << endl;
  }

  const bool converged = this->advance(a_dt, a_phi, m_source, a_zerophi);

  return converged;
}

void rte_solver::setRealm(const std::string a_realm){
  CH_TIME("rte_solver::setRealm");
  if(m_verbosity > 5){
    pout() << m_name + "::setRealm" << endl;
  }

  m_realm = a_realm;
}

void rte_solver::set_rte_species(const RefCountedPtr<rte_species> a_rte_species){
  CH_TIME("rte_solver::set_rte_species");
  if(m_verbosity > 5){
    pout() << m_name + "::set_rte_species" << endl;
  }

  m_rte_species = a_rte_species;
  m_name = m_rte_species->getName();
}

void rte_solver::setPhase(const phase::which_phase a_phase){
  CH_TIME("rte_solver::setPhase");
  if(m_verbosity > 5){
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

void rte_solver::sanityCheck(){
  CH_TIME("rte_solver::sanityCheck");
  if(m_verbosity > 5){
    pout() << m_name + "::sanityCheck" << endl;
  }

  CH_assert(!m_computationalGeometry.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_rte_species.isNull());
  CH_assert(!m_ebis.isNull());
}

void rte_solver::setComputationalGeometry(const RefCountedPtr<computational_geometry> a_computationalGeometry){
  CH_TIME("rte_solver::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << m_name + "::setComputationalGeometry" << endl;
  }
  
  m_computationalGeometry = a_computationalGeometry;

  const RefCountedPtr<mfis> mfis = m_computationalGeometry->get_mfis();
  
  this->setEbIndexSpace(mfis->getEBIndexSpace(m_phase));
}

void rte_solver::setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("rte_solver::setEbIndexSpace");
  if(m_verbosity > 5){
    pout() << m_name + "::setEbIndexSpace" << endl;
  }

  m_ebis = a_ebis;
}

void rte_solver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("rte_solver::setAmr");
  if(m_verbosity > 5){
    pout() << m_name + "::setAmr" << endl;
  }

  m_amr = a_amr;
}

void rte_solver::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("rte_solver::setTime");
  if(m_verbosity > 5){
    pout() << m_name + "::setTime" << endl;
  }

  m_timeStep = a_step;
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

void rte_solver::setVerbosity(const int a_verbosity){
  CH_TIME("rte_solver::setVerbosity");

  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << m_name + "::setVerbosity" << endl;
  }
  

}

void rte_solver::setSource(const EBAMRCellData& a_source){
  CH_TIME("rte_solver::setSource(ebamrcell)");
  if(m_verbosity > 5){
    pout() << m_name + "::setSource(ebamrcell)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_source[lvl]->localCopyTo(*m_source[lvl]);
  }

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void rte_solver::setSource(const Real a_source){
  CH_TIME("rte_solver::setSource(constant)");
  if(m_verbosity > 5){
    pout() << m_name + "::setSource(constant)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (int comp = 0; comp < m_source[lvl]->nComp(); comp++){
      data_ops::set_value(*m_source[lvl], a_source, comp);
    }
  }

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

int rte_solver::getNumberOfPlotVariables() const {
  CH_TIME("rte_solver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int num_output = 0;

  if(m_plotPhi) num_output = num_output + 1;
  if(m_plotSource) num_output = num_output + 1;

  return num_output;
}

void rte_solver::initialData() {
  CH_TIME("rte_solver::initialData");
  if(m_verbosity > 5){
    pout() << m_name + "::initialData" << endl;
  }
  
  data_ops::set_value(m_phi, 0.0);
}

void rte_solver::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("rte_solver::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  if(m_plotPhi) writeData(a_output, a_comp, m_phi,  true);
  if(m_plotSource) writeData(a_output, a_comp, m_source, false);
}

void rte_solver::writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp){
  CH_TIME("rte_solver::writeData");
  if(m_verbosity > 5){
    pout() << m_name + "::writeData" << endl;
  }

  const int comp = 0;
  const int ncomp = a_data[0]->nComp();

  const Interval src_interv(0, ncomp-1);
  const Interval dst_interv(a_comp, a_comp + ncomp - 1);

  // Copy data onto scratch
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, ncomp);
  data_ops::copy(scratch, a_data);

  // Interp if we should
  if(a_interp){
    m_amr->interpToCentroids(scratch, m_realm, phase::gas);
  }

  m_amr->averageDown(scratch, m_realm, m_phase);
  m_amr->interpGhost(scratch, m_realm, m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(a_output.getRealm() == m_realm){
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

EBAMRCellData& rte_solver::getPhi(){
  return m_phi;
}

EBAMRCellData& rte_solver::getSource(){
  return m_source;
}

EBAMRFluxData& rte_solver::get_kappa(){
  return m_kappa;
}

EBAMRIVData& rte_solver::get_kappa_eb(){
  return m_kappa_eb;
}

RefCountedPtr<rte_species>& rte_solver::getSpecies(){
  return m_rte_species;
}
#include "CD_NamespaceFooter.H"
