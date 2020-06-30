/*!
  @file   mf_realm.cpp
  @brief  Implementation of mf_realm.H
  @author Robert Marskar
  @date   July 2020
*/

#include "mf_realm.H"

mf_realm::mf_realm(){
  m_defined = false;
  m_verbosity = 10; 

  // Just empty points until define() is called
  m_realms[phase::gas]   = RefCountedPtr<realm> (new realm());
  m_realms[phase::solid] = RefCountedPtr<realm> (new realm());
}

mf_realm::~mf_realm(){

}

void mf_realm::define(const Vector<DisjointBoxLayout>& a_grids,
		      const Vector<ProblemDomain>& a_domains,
		      const Vector<int>& a_ref_rat,
		      const Vector<Real>& a_dx,
		      const int a_finest_level,
		      const int a_ebghost,
		      const int a_num_ghost,
		      const int a_redist_rad,
		      const bool a_ebcf,
		      const stencil_type::which_type a_centroid_stencil,
		      const stencil_type::which_type a_eb_stencil,
		      const RefCountedPtr<mfis>& a_mfis){
  CH_TIME("mf_realm::define");
  if(m_verbosity > 5){
    pout() << "mf_realm::define" << endl;
  }

  m_ref_ratios = a_ref_rat;
  m_dx = a_dx;
  m_grids = a_grids;
  m_domains = a_domains;
  m_mfis = a_mfis;
  m_finest_level = a_finest_level;

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  if(!ebis_gas.isNull()){
    m_realms[phase::gas]->define(a_grids, a_domains, a_ref_rat, a_dx, a_finest_level, a_ebghost, a_num_ghost, a_redist_rad,
				 a_centroid_stencil, a_eb_stencil, a_ebcf, ebis_gas);
  }

  if(!ebis_sol.isNull()){
    m_realms[phase::solid]->define(a_grids, a_domains, a_ref_rat, a_dx, a_finest_level, a_ebghost, a_num_ghost, a_redist_rad,
				   a_centroid_stencil, a_eb_stencil, a_ebcf, ebis_sol);
  }
}

void mf_realm::set_grids(const Vector<DisjointBoxLayout>& a_grids, const int a_finest_level){
  CH_TIME("mf_realm::set_grids");
  if(m_verbosity > 5){
    pout() << "mf_realm::set_grids" << endl;
  }

  for (auto& r : m_realms){
    r.second->set_grids(a_grids, a_finest_level);
  }
}

void mf_realm::regrid_base(const int a_lmin, const int a_lmax, const int a_hardcap){
  CH_TIME("mf_realm::regrid_base");
  if(m_verbosity > 5){
    pout() << "mf_realm::regrid_base" << endl;
  }
  
  for (auto& r : m_realms){
    r.second->regrid_base(a_lmin, a_lmax, a_hardcap);
  }
}

void mf_realm::regrid_operators(const int a_lmin, const int a_lmax, const int a_regsize){
  CH_TIME("mf_realm::regrid_operators");
  if(m_verbosity > 5){
    pout() << "mf_realm::regrid_operators" << endl;
  }
  
  for (auto& r : m_realms){
    r.second->regrid_operators(a_lmin, a_lmax, a_regsize);
  }
}

void mf_realm::define_mflevelgrid(const int a_lmin){
  CH_TIME("mf_realm::define_mflevelgrid");
  if(m_verbosity > 5){
    pout() << "mf_realm::define_mflevelgrid" << endl;
  }

  m_mflg.resize(1 + m_finest_level);

  realm& gas = this->get_realm(phase::gas);
  realm& sol = this->get_realm(phase::solid);

  const RefCountedPtr<EBIndexSpace>& ebis_gas = gas.get_ebis();
  const RefCountedPtr<EBIndexSpace>& ebis_sol = sol.get_ebis();

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBLevelGrid> eblgs;

    if(!ebis_gas.isNull()) eblgs.push_back(*(gas.get_eblg()[lvl]));
    if(!ebis_sol.isNull()) eblgs.push_back(*(sol.get_eblg()[lvl]));

    m_mflg[lvl] = RefCountedPtr<MFLevelGrid> (new MFLevelGrid(m_mfis, eblgs));
  }
}

void mf_realm::register_operator(const std::string a_operator, const phase::which_phase a_phase){
  CH_TIME("mf_realm::register_operator");
  if(m_verbosity > 5){
    pout() << "mf_realm::register_operator" << endl;
  }
  
  m_realms[a_phase]->register_operator(a_operator);
}

bool mf_realm::query_operator(const std::string a_operator, phase::which_phase a_phase){
  CH_TIME("mf_realm::query_operator");
  if(m_verbosity > 5){
    pout() << "mf_realm::query_operator" << endl;
  }
  
  return m_realms[a_phase]->query_operator(a_operator);
}

realm& mf_realm::get_realm(const phase::which_phase a_phase){
  return *m_realms[a_phase];
}

const Vector<int>& mf_realm::get_ref_rat() {
  return m_ref_ratios;
}

const Vector<Real>& mf_realm::get_dx() {
  return m_dx;
}

const Vector<DisjointBoxLayout>& mf_realm::get_grids() {
  return m_grids;
}

const Vector<ProblemDomain>& mf_realm::get_domains() {
  return m_domains;
}

const Vector<RefCountedPtr<MFLevelGrid> >& mf_realm::get_mflg(){
  return m_mflg;
}

const RefCountedPtr<EBIndexSpace>& mf_realm::get_ebis(const phase::which_phase a_phase) {

  return m_realms[a_phase]->get_ebis();
}

const Vector<EBISLayout>& mf_realm::get_ebisl(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_ebisl();
}

const Vector<RefCountedPtr<EBLevelGrid> >& mf_realm::get_eblg(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_eblg();
}

const Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& mf_realm::get_neighbors(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_neighbors();
}

const Vector<RefCountedPtr<LayoutData<VoFIterator> > >& mf_realm::get_vofit(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_vofit();
}

const irreg_amr_stencil<centroid_interp>& mf_realm::get_centroid_interp_stencils(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_centroid_interp_stencils();
}

const irreg_amr_stencil<eb_centroid_interp>& mf_realm::get_eb_centroid_interp_stencils(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_eb_centroid_interp_stencils();
}

const irreg_amr_stencil<noncons_div>& mf_realm::get_noncons_div_stencils(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_noncons_div_stencils();
}


Vector<RefCountedPtr<ebcoarseaverage> >& mf_realm::get_coarave(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_coarave();
}

Vector<RefCountedPtr<EBGhostCloud> >& mf_realm::get_ghostcloud(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_ghostcloud();
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& mf_realm::get_quadcfi(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_quadcfi();
}

Vector<RefCountedPtr<EBQuadCFInterp> >& mf_realm::get_old_quadcfi(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_old_quadcfi();
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& mf_realm::get_fillpatch(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_fillpatch();
}

Vector<RefCountedPtr<EBPWLFineInterp> >& mf_realm::get_eb_pwl_interp(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_eb_pwl_interp();
}

Vector<RefCountedPtr<EBMGInterp> >& mf_realm::get_eb_mg_interp(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_eb_mg_interp();
}

Vector<RefCountedPtr<EBFluxRegister> >&  mf_realm::get_flux_reg(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_flux_reg();
}

Vector<RefCountedPtr<EBLevelRedist> >&  mf_realm::get_level_redist(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_level_redist();
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  mf_realm::get_coar_to_fine_redist(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_coar_to_fine_redist();
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  mf_realm::get_coar_to_coar_redist(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_coar_to_coar_redist();
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  mf_realm::get_fine_to_coar_redist(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_fine_to_coar_redist();
}

Vector<RefCountedPtr<Copier> >& mf_realm::get_copier(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_copier();
}

Vector<RefCountedPtr<Copier> >& mf_realm::get_reverse_copier(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_reverse_copier();
}
