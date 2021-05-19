/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RtSolver.cpp
  @brief  Implementation of CD_RtSolver.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_RtSolver.H>
#include <data_ops.H>
#include <CD_NamespaceHeader.H>

RtSolver::RtSolver(){
  m_name      = "RtSolver";
  m_className = "RtSolver";
}

RtSolver::~RtSolver(){
  
}

std::string RtSolver::getName(){
  return m_name;
}

const std::string RtSolver::getRealm() const {
  return m_realm;
}

Vector<std::string> RtSolver::getPlotVariableNames() const {
  CH_TIME("RtSolver::getPlotVariableNames");
  if(m_verbosity > 5){
    pout() << m_name + "::getPlotVariableNames" << endl;
  }
  
  Vector<std::string> names(0);
  
  if(m_plotPhi) names.push_back(m_name + " phi");
  if(m_plotSource) names.push_back(m_name + " source");

  return names;
}

bool RtSolver::isStationary(){
  return m_stationary;
}

bool RtSolver::advance(const Real a_dt, const bool a_zerophi){
  CH_TIME("RtSolver::advance(dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(dt)" << endl;
  }

  const bool converged = this->advance(a_dt, m_phi, a_zerophi);

  return converged;
}

bool RtSolver::advance(const Real a_dt, EBAMRCellData& a_phi, const bool a_zerophi){
  CH_TIME("RtSolver::advance(dt, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(dt, state)" << endl;
  }

  const bool converged = this->advance(a_dt, a_phi, m_source, a_zerophi);

  return converged;
}

void RtSolver::setRealm(const std::string a_realm){
  CH_TIME("RtSolver::setRealm");
  if(m_verbosity > 5){
    pout() << m_name + "::setRealm" << endl;
  }

  m_realm = a_realm;
}

void RtSolver::setRtSpecies(const RefCountedPtr<RtSpecies> a_RtSpecies){
  CH_TIME("RtSolver::setRtSpecies");
  if(m_verbosity > 5){
    pout() << m_name + "::setRtSpecies" << endl;
  }

  m_RtSpecies = a_RtSpecies;
  m_name = m_RtSpecies->getName();
}

void RtSolver::setPhase(const phase::which_phase a_phase){
  CH_TIME("RtSolver::setPhase");
  if(m_verbosity > 5){
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

void RtSolver::sanityCheck(){
  CH_TIME("RtSolver::sanityCheck");
  if(m_verbosity > 5){
    pout() << m_name + "::sanityCheck" << endl;
  }

  CH_assert(!m_computationalGeometry.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_RtSpecies.isNull());
  CH_assert(!m_ebis.isNull());
}

void RtSolver::setComputationalGeometry(const RefCountedPtr<computational_geometry> a_computationalGeometry){
  CH_TIME("RtSolver::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << m_name + "::setComputationalGeometry" << endl;
  }
  
  m_computationalGeometry = a_computationalGeometry;

  const RefCountedPtr<mfis> mfis = m_computationalGeometry->get_mfis();
  
  this->setEbIndexSpace(mfis->getEBIndexSpace(m_phase));
}

void RtSolver::setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("RtSolver::setEbIndexSpace");
  if(m_verbosity > 5){
    pout() << m_name + "::setEbIndexSpace" << endl;
  }

  m_ebis = a_ebis;
}

void RtSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("RtSolver::setAmr");
  if(m_verbosity > 5){
    pout() << m_name + "::setAmr" << endl;
  }

  m_amr = a_amr;
}

void RtSolver::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("RtSolver::setTime");
  if(m_verbosity > 5){
    pout() << m_name + "::setTime" << endl;
  }

  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void RtSolver::setStationary(const bool a_stationary) {
  CH_TIME("RtSolver::setStationary");
  if(m_verbosity > 5){
    pout() << m_name + "::set_stationry" << endl;
  }
  
  m_stationary = a_stationary;
}

void RtSolver::setVerbosity(const int a_verbosity){
  CH_TIME("RtSolver::setVerbosity");

  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << m_name + "::setVerbosity" << endl;
  }
  

}

void RtSolver::setSource(const EBAMRCellData& a_source){
  CH_TIME("RtSolver::setSource(ebamrcell)");
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

void RtSolver::setSource(const Real a_source){
  CH_TIME("RtSolver::setSource(constant)");
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

int RtSolver::getNumberOfPlotVariables() const {
  CH_TIME("RtSolver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int num_output = 0;

  if(m_plotPhi) num_output = num_output + 1;
  if(m_plotSource) num_output = num_output + 1;

  return num_output;
}

void RtSolver::initialData() {
  CH_TIME("RtSolver::initialData");
  if(m_verbosity > 5){
    pout() << m_name + "::initialData" << endl;
  }
  
  data_ops::set_value(m_phi, 0.0);
}

void RtSolver::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("RtSolver::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  if(m_plotPhi) writeData(a_output, a_comp, m_phi,  true);
  if(m_plotSource) writeData(a_output, a_comp, m_source, false);
}

void RtSolver::writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp){
  CH_TIME("RtSolver::writeData");
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

Real RtSolver::getTime() const{
  CH_TIME("RtSolver::getTime");
  if(m_verbosity > 5){
    pout() << m_name + "::getTime" << endl;
  }
  
  return m_time;
}

phase::which_phase RtSolver::getPhase(){
  CH_TIME("RtSolver::getPhase");
  if(m_verbosity > 5){
    pout() << m_name + "::getPhase" << endl;
  }

  return m_phase;
}

EBAMRCellData& RtSolver::getPhi(){
  return m_phi;
}

EBAMRCellData& RtSolver::getSource(){
  return m_source;
}

EBAMRFluxData& RtSolver::getKappa(){
  return m_kappa;
}

EBAMRIVData& RtSolver::getKappaEb(){
  return m_kappa_eb;
}

RefCountedPtr<RtSpecies>& RtSolver::getSpecies(){
  return m_RtSpecies;
}

#include <CD_NamespaceFooter.H>
