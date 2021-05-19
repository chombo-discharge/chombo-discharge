/*!
  @file   cdr_layoutI.H
  @brief  Implementation of cdr_layout.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CDR_LAYOUTI_H
#define _CDR_LAYOUTI_H

#include "cdr_layout.H"
#include <CD_CdrIterator.H>
#include "units.H"

#include "CD_NamespaceHeader.H"

template <class T>
cdr_layout<T>::cdr_layout(){

}

template <class T>
cdr_layout<T>::cdr_layout(const Vector<RefCountedPtr<CdrSpecies> >& a_species){
  m_species = a_species;
  m_solvers.resize(0);
}

template <class T>
cdr_layout<T>::~cdr_layout(){

}

template <class T>
phase::which_phase cdr_layout<T>::getPhase() const{
  return m_phase;
}

template <class T>
CdrIterator<T> cdr_layout<T>::iterator() {
  return CdrIterator<T>(*this);
}

template <class T>
const std::string cdr_layout<T>::getRealm() const{
  return m_realm;
}

template <class T>
void cdr_layout<T>::setRealm(const std::string a_realm){
  m_realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->setRealm(a_realm);
  }
}

template <class T>
void cdr_layout<T>::addSolver(RefCountedPtr<T> a_solver){
  m_solvers.push_back(a_solver);
}

template <class T>
void cdr_layout<T>::parseOptions(){
  CH_TIME("cdr_layout<T>::parseOptions");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::parseOptions" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->parseOptions();
  }
}

template <class T>
void cdr_layout<T>::parseRuntimeOptions(){
  CH_TIME("cdr_layout<T>::parseRuntimeOptions");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::parseRuntimeOptions" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->parseRuntimeOptions();
  }
}

template <class T>
void cdr_layout<T>::allocateInternals(){
  CH_TIME("cdr_layout<T>::allocateInternals");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::allocateInternals" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->allocateInternals();
  }
}

template <class T>
void cdr_layout<T>::deallocateInternals(){
  CH_TIME("cdr_layout<T>::deallocateInternals");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::deallocateInternals" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->deallocateInternals();
  }
}

template <class T>
void cdr_layout<T>::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("cdr_layout<T>::preRegrid");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::preRegrid" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->preRegrid(a_lbase, a_oldFinestLevel);
  }
}

template <class T>
void cdr_layout<T>::initialData(){
  CH_TIME("cdr_layout<T>::initialData");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::initialData" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->initialData();
  }
}

template <class T>
void cdr_layout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("cdr_layout<T>::regrid");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::regrid" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void cdr_layout<T>::registerOperators(){
  CH_TIME("cdr_layout<T>::registerOperators");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::registerOperators" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->registerOperators();
  }
}

template <class T>
void cdr_layout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("cdr_layout<T>::setAmr");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setAmr" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setAmr(a_amr);
  }
}

template <class T>
void cdr_layout<T>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry){
  CH_TIME("cdr_layout<T>::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setComputationalGeometry" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void cdr_layout<T>::setPhase(phase::which_phase a_phase){
  CH_TIME("cdr_layout<T>::setPhase");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setPhase" << endl;
  }

  m_phase = a_phase;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setPhase(a_phase);
  }
}

template <class T>
void cdr_layout<T>::setVerbosity(const int a_verbosity){
  CH_TIME("cdr_layout<T>::setVerbosity");
  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setVerbosity" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setVerbosity(a_verbosity);
  }
}

template <class T>
void cdr_layout<T>::set_bc(CdrBc a_bc){
  CH_TIME("cdr_layout<T>::set_bc");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_bc" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setDomainBc(a_bc);
  }
}

template <class T>
void cdr_layout<T>::sanityCheck(){
  CH_TIME("cdr_layout<T>::sanityCheck");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::sanityCheck" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->sanityCheck();
  }
}

template <class T>
void cdr_layout<T>::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("cdr_layout<T>::setTime");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setTime" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setTime(a_step, a_time, a_dt);
  }
}

template <class T>
void cdr_layout<T>::setVelocity(const EBAMRCellData& a_velo){
  CH_TIME("cdr_layout<T>::setVelocity(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setVelocity(ebamrcell)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setVelocity(a_velo);
  }
}

template <class T>
void cdr_layout<T>::setVelocity(const RealVect a_velo){
  CH_TIME("cdr_layout<T>::setVelocity(constant)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setVelocity(constant)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setVelocity(a_velo);
  }
}

template <class T>
void cdr_layout<T>::setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient, const EBAMRIVData& a_ebDiffusionCoefficient){
  CH_TIME("cdr_layout<T>::setDiffusionCoefficient(ebamrflux)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setDiffusionCoefficient(ebamrflux)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setDiffusionCoefficient(a_diffusionCoefficient, a_ebDiffusionCoefficient);
  }
}

template <class T>
void cdr_layout<T>::setDiffusionCoefficient(const Real a_diffusionCoefficient){
  CH_TIME("cdr_layout<T>::setDiffusionCoefficient(constant)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setDiffusionCoefficient(constant)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setDiffusionCoefficient(a_diffusionCoefficient);
  }
}

template <class T>
void cdr_layout<T>::setSource(const EBAMRCellData& a_source){
  CH_TIME("cdr_layout<T>::setSource(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setSource(ebamrcell)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setSource(a_source);
  }
}

template <class T>
void cdr_layout<T>::setSource(const Real a_source){
  CH_TIME("cdr_layout<T>::setSource(constant)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setSource(constant)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setSource(a_source);
  }
}

template <class T>
void cdr_layout<T>::setEbFlux(const EBAMRIVData& a_ebFlux){
  CH_TIME("cdr_layout<T>::setSource(ebamriv)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setSource(ebamriv)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setEbFlux(a_ebFlux);
  }
}

template <class T>
void cdr_layout<T>::setEbFlux(const Real a_ebFlux){
  CH_TIME("cdr_layout<T>::setSource(constant)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setSource(constant)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->setEbFlux(a_ebFlux);
  }
}

template <class T>
void cdr_layout<T>::writePlotFile(){
  CH_TIME("cdr_layout<T>::writePlotFile");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::writePlotFile" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->writePlotFile();
  }
}

template <class T>
Real cdr_layout<T>::computeAdvectionDt(){
  CH_TIME("cdr_layout<T>::computeAdvectionDt");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::computeAdvectionDt" << endl;
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
Real cdr_layout<T>::computeDiffusionDt(){
  CH_TIME("cdr_layout<T>::computeDiffusionDt");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::computeDiffusionDt" << endl;
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
Real cdr_layout<T>::computeAdvectionDiffusionDt(){
  CH_TIME("cdr_layout<T>::computeAdvectionDiffusionDt");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::computeAdvectionDiffusionDt" << endl;
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
Real cdr_layout<T>::computeSourceDt(const Real a_tolerance, const bool a_elec_only){
  CH_TIME("cdr_layout<T>::computeSourceDt");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::computeSourceDt" << endl;
  }

  const int comp = 0;
  
  Real max = 0.0;
  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    Real cur_max, cur_min;
    RefCountedPtr<CdrSolver>& solver = solver_it();
    data_ops::get_max_min(cur_max, cur_min, solver->getPhi(), comp);
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
Real cdr_layout<T>::compute_total_charge(){
  CH_TIME("cdr_layout<T>::compute_total_charge");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::compute_total_charge" << endl;
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
Vector<Real> cdr_layout<T>::computeMass(){
  CH_TIME("cdr_layout<T>::computeMass");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::computeMass" << endl;
  }

  Vector<Real> mass;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    mass.push_back(solver->computeMass());
  }

  return mass;
}

template <class T>
Vector<Real> cdr_layout<T>::computeCharge(){
  CH_TIME("cdr_layout<T>::computeCharge");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::computeCharge" << endl;
  }

  Vector<Real> charge;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    charge.push_back(solver->computeCharge()*units::s_Qe);
  }

  return charge;
}

template <class T>
Vector<std::string> cdr_layout<T>::getNames() {
  CH_TIME("cdr_layout<T>::getNames");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::getNames" << endl;
  }

  Vector<std::string> names;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    names.push_back(solver->getName());
  }

  return names;
}

template <class T>
Vector<RefCountedPtr<T> >& cdr_layout<T>::getSolvers(){
  return m_solvers;
}

template <class T>
Vector<RefCountedPtr<CdrSpecies> >& cdr_layout<T>::getSpecies(){
  return m_species;
}

template <class T>
Vector<EBAMRCellData*> cdr_layout<T>::getPhis(){
  CH_TIME("cdr_layout<T>::getPhis");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::getPhis" << endl;
  }

  Vector<EBAMRCellData*> states;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    states.push_back(&(solver->getPhi()));
  }

  return states;
}

template <class T>
Vector<EBAMRCellData*> cdr_layout<T>::getSources(){
  CH_TIME("cdr_layout<T>::getSources");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::getSources" << endl;
  }

  Vector<EBAMRCellData*> sources;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    sources.push_back(&(solver->getSource()));
  }

  return sources;
}

template <class T>
Vector<EBAMRCellData*> cdr_layout<T>::get_velocities(){
  CH_TIME("cdr_layout<T>::get_velocities");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::get_velocities" << endl;
  }

  Vector<EBAMRCellData*> velocities;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    velocities.push_back(&(solver->getCellCenteredVelocity()));
  }

  return velocities;
}

template <class T>
Vector<EBAMRFluxData*> cdr_layout<T>::getFaceCenteredDiffusionCoefficient(){
  CH_TIME("cdr_layout<T>::getFaceCenteredDiffusionCoefficient");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::getFaceCenteredDiffusionCoefficient" << endl;
  }

  Vector<EBAMRFluxData*> diffco;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    diffco.push_back(&(solver->getFaceCenteredDiffusionCoefficient()));
  }

  return diffco;
}

template <class T>
Vector<EBAMRIVData*> cdr_layout<T>::getEbCenteredDiffusionCoefficient(){
  CH_TIME("cdr_layout<T>::getEbCenteredDiffusionCoefficient");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::getEbCenteredDiffusionCoefficient" << endl;
  }

  Vector<EBAMRIVData*> diffco_eb;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    diffco_eb.push_back(&(solver->getEbCenteredDiffusionCoefficient()));
  }

  return diffco_eb;
}

template <class T>
Vector<EBAMRIVData*> cdr_layout<T>::getEbFlux(){
  CH_TIME("cdr_layout<T>::getEbFlux");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::getEbFlux" << endl;
  }

  Vector<EBAMRIVData*> ebflux;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    ebflux.push_back(&(solver->getEbFlux()));
  }

  return ebflux;
}

template <class T>
Vector<EBAMRIFData*> cdr_layout<T>::getDomainFlux(){
  CH_TIME("cdr_layout<T>::getDomainFlux");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::getDomainFlux" << endl;
  }

  Vector<EBAMRIFData*> domainflux;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    domainflux.push_back(&(solver->getDomainFlux()));
  }

  return domainflux;
}

template <class T>
bool cdr_layout<T>::get_mass_redist(){
  return m_mass_redist;
}


template <class T, class S>
cdr_factory<T, S>::cdr_factory(){
}

template <class T, class S>
cdr_factory<T, S>::~cdr_factory(){
}


template <class T, class S>
RefCountedPtr<cdr_layout<T> > cdr_factory<T, S>::newLayout(const Vector<RefCountedPtr<CdrSpecies> > a_species) const{

  // Build cdr layout
  auto cdr = RefCountedPtr<cdr_layout<T> > (new cdr_layout<T>(a_species));

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
#include "CD_NamespaceFooter.H"

#endif
