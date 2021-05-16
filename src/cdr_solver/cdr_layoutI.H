/*!
  @file   cdr_layoutI.H
  @brief  Implementation of cdr_layout.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CDR_LAYOUTI_H
#define _CDR_LAYOUTI_H

#include "cdr_layout.H"
#include "cdr_iterator.H"
#include "units.H"

#include "CD_NamespaceHeader.H"

template <class T>
cdr_layout<T>::cdr_layout(){

}

template <class T>
cdr_layout<T>::cdr_layout(const Vector<RefCountedPtr<cdr_species> >& a_species){
  m_species = a_species;
  m_solvers.resize(0);
}

template <class T>
cdr_layout<T>::~cdr_layout(){

}

template <class T>
phase::which_phase cdr_layout<T>::get_phase() const{
  return m_phase;
}

template <class T>
cdr_iterator<T> cdr_layout<T>::iterator() {
  return cdr_iterator<T>(*this);
}

template <class T>
const std::string cdr_layout<T>::get_Realm() const{
  return m_Realm;
}

template <class T>
void cdr_layout<T>::set_Realm(const std::string a_Realm){
  m_Realm = a_Realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->set_Realm(a_Realm);
  }
}

template <class T>
void cdr_layout<T>::add_solver(RefCountedPtr<T> a_solver){
  m_solvers.push_back(a_solver);
}

template <class T>
void cdr_layout<T>::parseOptions(){
  CH_TIME("cdr_layout<T>::parseOptions");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::parseOptions" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->parseOptions();
  }
}

template <class T>
void cdr_layout<T>::parseRuntimeOptions(){
  CH_TIME("cdr_layout<T>::parseRuntimeOptions");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::parseRuntimeOptions" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->parseRuntimeOptions();
  }
}

template <class T>
void cdr_layout<T>::allocateInternals(){
  CH_TIME("cdr_layout<T>::allocateInternals");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::allocateInternals" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->allocateInternals();
  }
}

template <class T>
void cdr_layout<T>::deallocateInternals(){
  CH_TIME("cdr_layout<T>::deallocateInternals");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::deallocateInternals" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->deallocateInternals();
  }
}

template <class T>
void cdr_layout<T>::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("cdr_layout<T>::preRegrid");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::preRegrid" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->preRegrid(a_lbase, a_oldFinestLevel);
  }
}

template <class T>
void cdr_layout<T>::initialData(){
  CH_TIME("cdr_layout<T>::initialData");
  if(m_verbosity > 6){
    pout() << "cdr_layout<T>::initialData" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->initialData();
  }
}

template <class T>
void cdr_layout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("cdr_layout<T>::regrid");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::regrid" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void cdr_layout<T>::registerOperators(){
  CH_TIME("cdr_layout<T>::registerOperators");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::registerOperators" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->registerOperators();
  }
}

template <class T>
void cdr_layout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("cdr_layout<T>::setAmr");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setAmr" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->setAmr(a_amr);
  }
}

template <class T>
void cdr_layout<T>::setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_computationalGeometry){
  CH_TIME("cdr_layout<T>::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::setComputationalGeometry" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void cdr_layout<T>::set_phase(phase::which_phase a_phase){
  CH_TIME("cdr_layout<T>::set_phase");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_phase" << endl;
  }

  m_phase = a_phase;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_phase(a_phase);
  }
}

template <class T>
void cdr_layout<T>::set_verbosity(const int a_verbosity){
  CH_TIME("cdr_layout<T>::set_verbosity");
  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_verbosity" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_verbosity(a_verbosity);
  }
}

template <class T>
void cdr_layout<T>::set_bc(cdr_bc::which_bc a_bc){
  CH_TIME("cdr_layout<T>::set_bc");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_bc" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_domain_bc(a_bc);
  }
}

template <class T>
void cdr_layout<T>::sanityCheck(){
  CH_TIME("cdr_layout<T>::sanityCheck");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::sanityCheck" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->sanityCheck();
  }
}

template <class T>
void cdr_layout<T>::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("cdr_layout<T>::set_time");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_time" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_time(a_step, a_time, a_dt);
  }
}

template <class T>
void cdr_layout<T>::set_velocity(const EBAMRCellData& a_velo){
  CH_TIME("cdr_layout<T>::set_velocity(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_velocity(ebamrcell)" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_velocity(a_velo);
  }
}

template <class T>
void cdr_layout<T>::set_velocity(const RealVect a_velo){
  CH_TIME("cdr_layout<T>::set_velocity(constant)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_velocity(constant)" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_velocity(a_velo);
  }
}

template <class T>
void cdr_layout<T>::set_diffco(const EBAMRFluxData& a_diffco, const EBAMRIVData& a_diffco_eb){
  CH_TIME("cdr_layout<T>::set_diffco(ebamrflux)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_diffco(ebamrflux)" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_diffco(a_diffco, a_diffco_eb);
  }
}

template <class T>
void cdr_layout<T>::set_diffco(const Real a_diffco){
  CH_TIME("cdr_layout<T>::set_diffco(constant)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_diffco(constant)" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_diffco(a_diffco);
  }
}

template <class T>
void cdr_layout<T>::set_source(const EBAMRCellData& a_source){
  CH_TIME("cdr_layout<T>::set_source(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_source(ebamrcell)" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_source(a_source);
  }
}

template <class T>
void cdr_layout<T>::set_source(const Real a_source){
  CH_TIME("cdr_layout<T>::set_source(constant)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_source(constant)" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_source(a_source);
  }
}

template <class T>
void cdr_layout<T>::set_ebflux(const EBAMRIVData& a_ebflux){
  CH_TIME("cdr_layout<T>::set_source(ebamriv)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_source(ebamriv)" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_ebflux(a_ebflux);
  }
}

template <class T>
void cdr_layout<T>::set_ebflux(const Real a_ebflux){
  CH_TIME("cdr_layout<T>::set_source(constant)");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::set_source(constant)" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_ebflux(a_ebflux);
  }
}

template <class T>
void cdr_layout<T>::writePlotFile(){
  CH_TIME("cdr_layout<T>::writePlotFile");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::writePlotFile" << endl;
  }

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->writePlotFile();
  }
}

template <class T>
Real cdr_layout<T>::compute_advection_dt(){
  CH_TIME("cdr_layout<T>::compute_advection_dt");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::compute_advection_dt" << endl;
  }

  Real dt = 1.E99;
  
  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const Real this_dt = solver->compute_advection_dt();

    dt = Min(dt, this_dt);
  }

  return dt;
}

template <class T>
Real cdr_layout<T>::compute_diffusion_dt(){
  CH_TIME("cdr_layout<T>::compute_diffusion_dt");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::compute_diffusion_dt" << endl;
  }

  Real dt = 1.E99;
  
  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const Real this_dt = solver->compute_diffusion_dt();

    dt = Min(dt, this_dt);
  }

  return dt;
}

template <class T>
Real cdr_layout<T>::compute_advection_diffusion_dt(){
  CH_TIME("cdr_layout<T>::compute_advection_diffusion_dt");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::compute_advection_diffusion_dt" << endl;
  }

  Real dt = 1.E99;
  
  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const Real this_dt = solver->compute_advection_diffusion_dt();

    dt = Min(dt, this_dt);
  }

  return dt;
}

template <class T>
Real cdr_layout<T>::compute_source_dt(const Real a_tolerance, const bool a_elec_only){
  CH_TIME("cdr_layout<T>::compute_source_dt");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::compute_source_dt" << endl;
  }

  const int comp = 0;
  
  Real max = 0.0;
  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    Real cur_max, cur_min;
    RefCountedPtr<cdr_solver>& solver = solver_it();
    data_ops::get_max_min(cur_max, cur_min, solver->get_state(), comp);
    max = Max(max, cur_max);
    if(a_elec_only){
      break;
    }
  }

  Real dt = 1.E99;
  
  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const Real this_dt = solver->compute_source_dt(max, a_tolerance);

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
  
  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const Real N = solver->compute_mass();

    Q += N*units::s_Qe;
  }

  return Q;
}

template <class T>
Vector<Real> cdr_layout<T>::compute_mass(){
  CH_TIME("cdr_layout<T>::compute_mass");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::compute_mass" << endl;
  }

  Vector<Real> mass;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    mass.push_back(solver->compute_mass());
  }

  return mass;
}

template <class T>
Vector<Real> cdr_layout<T>::compute_charge(){
  CH_TIME("cdr_layout<T>::compute_charge");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::compute_charge" << endl;
  }

  Vector<Real> charge;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    charge.push_back(solver->compute_charge()*units::s_Qe);
  }

  return charge;
}

template <class T>
Vector<std::string> cdr_layout<T>::get_names() {
  CH_TIME("cdr_layout<T>::get_names");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::get_names" << endl;
  }

  Vector<std::string> names;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    names.push_back(solver->get_name());
  }

  return names;
}

template <class T>
Vector<RefCountedPtr<T> >& cdr_layout<T>::get_solvers(){
  return m_solvers;
}

template <class T>
Vector<RefCountedPtr<cdr_species> >& cdr_layout<T>::get_species(){
  return m_species;
}

template <class T>
Vector<EBAMRCellData*> cdr_layout<T>::get_states(){
  CH_TIME("cdr_layout<T>::get_states");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::get_states" << endl;
  }

  Vector<EBAMRCellData*> states;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    states.push_back(&(solver->get_state()));
  }

  return states;
}

template <class T>
Vector<EBAMRCellData*> cdr_layout<T>::get_sources(){
  CH_TIME("cdr_layout<T>::get_sources");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::get_sources" << endl;
  }

  Vector<EBAMRCellData*> sources;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    sources.push_back(&(solver->get_source()));
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

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    velocities.push_back(&(solver->get_velo_cell()));
  }

  return velocities;
}

template <class T>
Vector<EBAMRFluxData*> cdr_layout<T>::get_diffco_face(){
  CH_TIME("cdr_layout<T>::get_diffco_face");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::get_diffco_face" << endl;
  }

  Vector<EBAMRFluxData*> diffco;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    diffco.push_back(&(solver->get_diffco_face()));
  }

  return diffco;
}

template <class T>
Vector<EBAMRIVData*> cdr_layout<T>::get_diffco_eb(){
  CH_TIME("cdr_layout<T>::get_diffco_eb");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::get_diffco_eb" << endl;
  }

  Vector<EBAMRIVData*> diffco_eb;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    diffco_eb.push_back(&(solver->get_diffco_eb()));
  }

  return diffco_eb;
}

template <class T>
Vector<EBAMRIVData*> cdr_layout<T>::get_ebflux(){
  CH_TIME("cdr_layout<T>::get_ebflux");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::get_ebflux" << endl;
  }

  Vector<EBAMRIVData*> ebflux;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    ebflux.push_back(&(solver->get_ebflux()));
  }

  return ebflux;
}

template <class T>
Vector<EBAMRIFData*> cdr_layout<T>::get_domainflux(){
  CH_TIME("cdr_layout<T>::get_domainflux");
  if(m_verbosity > 5){
    pout() << "cdr_layout<T>::get_domainflux" << endl;
  }

  Vector<EBAMRIFData*> domainflux;

  for (cdr_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    domainflux.push_back(&(solver->get_domainflux()));
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
RefCountedPtr<cdr_layout<T> > cdr_factory<T, S>::new_layout(const Vector<RefCountedPtr<cdr_species> > a_species) const{

  // Build cdr layout
  auto cdr = RefCountedPtr<cdr_layout<T> > (new cdr_layout<T>(a_species));

  // Cast solvers and instantiate them
  for (int i = 0; i < a_species.size(); i++){
    auto solver = RefCountedPtr<T> (static_cast<T*> (new S()));
    solver->set_species(a_species[i]);
    solver->set_phase(phase::gas);
    solver->set_verbosity(-1);
    cdr->add_solver(solver);
  }

  return cdr;
}
#include "CD_NamespaceFooter.H"

#endif
