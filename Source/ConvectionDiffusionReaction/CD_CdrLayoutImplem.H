/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrLayoutImplem.H
  @brief  Implementation of CD_CdrLayout.H
  @author Robert Marskar
*/

#ifndef CD_CdrLayoutImplem_H
#define CD_CdrLayoutImplem_H

// Our includes
#include <CD_CdrLayout.H>
#include <CD_CdrIterator.H>
#include <units.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

template <class T>
CdrLayout<T>::CdrLayout(){

}

template <class T>
CdrLayout<T>::CdrLayout(const Vector<RefCountedPtr<CdrSpecies> >& a_species){
  m_species = a_species;
  m_solvers.resize(0);
}

template <class T>
CdrLayout<T>::~CdrLayout(){

}

template <class T>
phase::which_phase CdrLayout<T>::getPhase() const{
  return m_phase;
}

template <class T>
CdrIterator<T> CdrLayout<T>::iterator() {
  return CdrIterator<T>(*this);
}

template <class T>
const std::string CdrLayout<T>::getRealm() const{
  return m_realm;
}

template <class T>
void CdrLayout<T>::setRealm(const std::string a_realm){
  m_realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->setRealm(a_realm);
  }
}

template <class T>
void CdrLayout<T>::addSolver(RefCountedPtr<T> a_solver){
  m_solvers.push_back(a_solver);
}

template <class T>
void CdrLayout<T>::parseOptions(){
  CH_TIME("CdrLayout<T>::parseOptions");
  if(m_verbosity > 6){
    pout() << "CdrLayout<T>::parseOptions" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->parseOptions();
  }
}

template <class T>
void CdrLayout<T>::parseRuntimeOptions(){
  CH_TIME("CdrLayout<T>::parseRuntimeOptions");
  if(m_verbosity > 6){
    pout() << "CdrLayout<T>::parseRuntimeOptions" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->parseRuntimeOptions();
  }
}

template <class T>
void CdrLayout<T>::allocateInternals(){
  CH_TIME("CdrLayout<T>::allocateInternals");
  if(m_verbosity > 6){
    pout() << "CdrLayout<T>::allocateInternals" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->allocateInternals();
  }
}

template <class T>
void CdrLayout<T>::deallocateInternals(){
  CH_TIME("CdrLayout<T>::deallocateInternals");
  if(m_verbosity > 6){
    pout() << "CdrLayout<T>::deallocateInternals" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->deallocateInternals();
  }
}

template <class T>
void CdrLayout<T>::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("CdrLayout<T>::preRegrid");
  if(m_verbosity > 6){
    pout() << "CdrLayout<T>::preRegrid" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->preRegrid(a_lbase, a_oldFinestLevel);
  }
}

template <class T>
void CdrLayout<T>::initialData(){
  CH_TIME("CdrLayout<T>::initialData");
  if(m_verbosity > 6){
    pout() << "CdrLayout<T>::initialData" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->initialData();
  }
}

template <class T>
void CdrLayout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("CdrLayout<T>::regrid");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::regrid" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void CdrLayout<T>::registerOperators(){
  CH_TIME("CdrLayout<T>::registerOperators");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::registerOperators" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->registerOperators();
  }
}

template <class T>
void CdrLayout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("CdrLayout<T>::setAmr");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setAmr" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setAmr(a_amr);
  }
}

template <class T>
void CdrLayout<T>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry){
  CH_TIME("CdrLayout<T>::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setComputationalGeometry" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void CdrLayout<T>::setPhase(phase::which_phase a_phase){
  CH_TIME("CdrLayout<T>::setPhase");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setPhase" << endl;
  }

  m_phase = a_phase;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setPhase(a_phase);
  }
}

template <class T>
void CdrLayout<T>::setVerbosity(const int a_verbosity){
  CH_TIME("CdrLayout<T>::setVerbosity");
  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setVerbosity" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setVerbosity(a_verbosity);
  }
}

template <class T>
void CdrLayout<T>::set_bc(CdrBc a_bc){
  CH_TIME("CdrLayout<T>::set_bc");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::set_bc" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setDomainBc(a_bc);
  }
}

template <class T>
void CdrLayout<T>::sanityCheck(){
  CH_TIME("CdrLayout<T>::sanityCheck");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::sanityCheck" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->sanityCheck();
  }
}

template <class T>
void CdrLayout<T>::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("CdrLayout<T>::setTime");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setTime" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setTime(a_step, a_time, a_dt);
  }
}

template <class T>
void CdrLayout<T>::setVelocity(const EBAMRCellData& a_velo){
  CH_TIME("CdrLayout<T>::setVelocity(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setVelocity(ebamrcell)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setVelocity(a_velo);
  }
}

template <class T>
void CdrLayout<T>::setVelocity(const RealVect a_velo){
  CH_TIME("CdrLayout<T>::setVelocity(constant)");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setVelocity(constant)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setVelocity(a_velo);
  }
}

template <class T>
void CdrLayout<T>::setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient, const EBAMRIVData& a_ebDiffusionCoefficient){
  CH_TIME("CdrLayout<T>::setDiffusionCoefficient(ebamrflux)");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setDiffusionCoefficient(ebamrflux)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setDiffusionCoefficient(a_diffusionCoefficient, a_ebDiffusionCoefficient);
  }
}

template <class T>
void CdrLayout<T>::setDiffusionCoefficient(const Real a_diffusionCoefficient){
  CH_TIME("CdrLayout<T>::setDiffusionCoefficient(constant)");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setDiffusionCoefficient(constant)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setDiffusionCoefficient(a_diffusionCoefficient);
  }
}

template <class T>
void CdrLayout<T>::setSource(const EBAMRCellData& a_source){
  CH_TIME("CdrLayout<T>::setSource(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setSource(ebamrcell)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setSource(a_source);
  }
}

template <class T>
void CdrLayout<T>::setSource(const Real a_source){
  CH_TIME("CdrLayout<T>::setSource(constant)");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setSource(constant)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setSource(a_source);
  }
}

template <class T>
void CdrLayout<T>::setEbFlux(const EBAMRIVData& a_ebFlux){
  CH_TIME("CdrLayout<T>::setSource(ebamriv)");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setSource(ebamriv)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setEbFlux(a_ebFlux);
  }
}

template <class T>
void CdrLayout<T>::setEbFlux(const Real a_ebFlux){
  CH_TIME("CdrLayout<T>::setSource(constant)");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::setSource(constant)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setEbFlux(a_ebFlux);
  }
}

template <class T>
void CdrLayout<T>::writePlotFile(){
  CH_TIME("CdrLayout<T>::writePlotFile");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::writePlotFile" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->writePlotFile();
  }
}

template <class T>
Real CdrLayout<T>::computeAdvectionDt(){
  CH_TIME("CdrLayout<T>::computeAdvectionDt");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::computeAdvectionDt" << endl;
  }

  Real dt = 1.E99;
  
  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const Real this_dt = solver->computeAdvectionDt();

    dt = Min(dt, this_dt);
  }

  return dt;
}

template <class T>
Real CdrLayout<T>::computeDiffusionDt(){
  CH_TIME("CdrLayout<T>::computeDiffusionDt");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::computeDiffusionDt" << endl;
  }

  Real dt = 1.E99;
  
  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const Real this_dt = solver->computeDiffusionDt();

    dt = Min(dt, this_dt);
  }

  return dt;
}

template <class T>
Real CdrLayout<T>::computeAdvectionDiffusionDt(){
  CH_TIME("CdrLayout<T>::computeAdvectionDiffusionDt");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::computeAdvectionDiffusionDt" << endl;
  }

  Real dt = 1.E99;
  
  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const Real this_dt = solver->computeAdvectionDiffusionDt();

    dt = Min(dt, this_dt);
  }

  return dt;
}

template <class T>
Real CdrLayout<T>::computeSourceDt(const Real a_tolerance, const bool a_elec_only){
  CH_TIME("CdrLayout<T>::computeSourceDt");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::computeSourceDt" << endl;
  }

  const int comp = 0;
  
  Real max = 0.0;
  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    Real cur_max, cur_min;
    RefCountedPtr<CdrSolver>& solver = solver_it();
    DataOps::getMaxMin(cur_max, cur_min, solver->getPhi(), comp);
    max = Max(max, cur_max);
    if(a_elec_only){
      break;
    }
  }

  Real dt = 1.E99;
  
  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const Real this_dt = solver->computeSourceDt(max, a_tolerance);

    dt = Min(dt, this_dt);
    if(a_elec_only){
      break;
    }
  }

  return dt;
}

template <class T>
Real CdrLayout<T>::computeTotalCharge(){
  CH_TIME("CdrLayout<T>::computeTotalCharge");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::computeTotalCharge" << endl;
  }

  Real Q = 0.;
  
  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const Real N = solver->computeMass();

    Q += N*units::s_Qe;
  }

  return Q;
}

template <class T>
Vector<Real> CdrLayout<T>::computeMass(){
  CH_TIME("CdrLayout<T>::computeMass");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::computeMass" << endl;
  }

  Vector<Real> mass;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    mass.push_back(solver->computeMass());
  }

  return mass;
}

template <class T>
Vector<Real> CdrLayout<T>::computeCharge(){
  CH_TIME("CdrLayout<T>::computeCharge");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::computeCharge" << endl;
  }

  Vector<Real> charge;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    charge.push_back(solver->computeCharge()*units::s_Qe);
  }

  return charge;
}

template <class T>
Vector<std::string> CdrLayout<T>::getNames() {
  CH_TIME("CdrLayout<T>::getNames");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::getNames" << endl;
  }

  Vector<std::string> names;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    names.push_back(solver->getName());
  }

  return names;
}

template <class T>
Vector<RefCountedPtr<T> >& CdrLayout<T>::getSolvers(){
  return m_solvers;
}

template <class T>
Vector<RefCountedPtr<CdrSpecies> >& CdrLayout<T>::getSpecies(){
  return m_species;
}

template <class T>
Vector<EBAMRCellData*> CdrLayout<T>::getPhis(){
  CH_TIME("CdrLayout<T>::getPhis");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::getPhis" << endl;
  }

  Vector<EBAMRCellData*> states;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    states.push_back(&(solver->getPhi()));
  }

  return states;
}

template <class T>
Vector<EBAMRCellData*> CdrLayout<T>::getSources(){
  CH_TIME("CdrLayout<T>::getSources");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::getSources" << endl;
  }

  Vector<EBAMRCellData*> sources;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    sources.push_back(&(solver->getSource()));
  }

  return sources;
}

template <class T>
Vector<EBAMRCellData*> CdrLayout<T>::getVelocities(){
  CH_TIME("CdrLayout<T>::getVelocities");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::getVelocities" << endl;
  }

  Vector<EBAMRCellData*> velocities;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    velocities.push_back(&(solver->getCellCenteredVelocity()));
  }

  return velocities;
}

template <class T>
Vector<EBAMRFluxData*> CdrLayout<T>::getFaceCenteredDiffusionCoefficient(){
  CH_TIME("CdrLayout<T>::getFaceCenteredDiffusionCoefficient");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::getFaceCenteredDiffusionCoefficient" << endl;
  }

  Vector<EBAMRFluxData*> diffco;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    diffco.push_back(&(solver->getFaceCenteredDiffusionCoefficient()));
  }

  return diffco;
}

template <class T>
Vector<EBAMRIVData*> CdrLayout<T>::getEbCenteredDiffusionCoefficient(){
  CH_TIME("CdrLayout<T>::getEbCenteredDiffusionCoefficient");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::getEbCenteredDiffusionCoefficient" << endl;
  }

  Vector<EBAMRIVData*> diffco_eb;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    diffco_eb.push_back(&(solver->getEbCenteredDiffusionCoefficient()));
  }

  return diffco_eb;
}

template <class T>
Vector<EBAMRIVData*> CdrLayout<T>::getEbFlux(){
  CH_TIME("CdrLayout<T>::getEbFlux");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::getEbFlux" << endl;
  }

  Vector<EBAMRIVData*> ebflux;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    ebflux.push_back(&(solver->getEbFlux()));
  }

  return ebflux;
}

template <class T>
Vector<EBAMRIFData*> CdrLayout<T>::getDomainFlux(){
  CH_TIME("CdrLayout<T>::getDomainFlux");
  if(m_verbosity > 5){
    pout() << "CdrLayout<T>::getDomainFlux" << endl;
  }

  Vector<EBAMRIFData*> domainflux;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    domainflux.push_back(&(solver->getDomainFlux()));
  }

  return domainflux;
}

template <class T>
bool CdrLayout<T>::getMassRedist(){
  return m_mass_redist;
}


template <class T, class S>
CdrFactory<T, S>::CdrFactory(){
}

template <class T, class S>
CdrFactory<T, S>::~CdrFactory(){
}


template <class T, class S>
RefCountedPtr<CdrLayout<T> > CdrFactory<T, S>::newLayout(const Vector<RefCountedPtr<CdrSpecies> > a_species) const{

  // Build cdr layout
  auto cdr = RefCountedPtr<CdrLayout<T> > (new CdrLayout<T>(a_species));

  // Cast solvers and instantiate them
  for (int i = 0; i < a_species.size(); i++){
    auto solver = RefCountedPtr<T> (static_cast<T*> (new S()));
    solver->setSpecies(a_species[i]);
    solver->setPhase(phase::gas);
    solver->setVerbosity(-1);
    cdr->addSolver(solver);
  }

  return cdr;
}

#include <CD_NamespaceFooter.H>

#endif
