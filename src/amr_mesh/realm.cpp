/*!
  @file   realm.cpp
  @brief  Implementation of realm.H
  @author Robert Marskar
  @date   July 2020
*/

#include "realm.H"

const std::string realm::primal = "primal";

template <class T>
void realm::copy(Vector<RefCountedPtr<LevelData<T> > >& a_dst, const Vector<RefCountedPtr<LevelData<T> > >& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->copyTo(a_dst[lvl]);
    }
  }
}

realm::realm(){
  m_defined = false;
  m_verbosity = -1;

  // Just empty points until define() is called
  m_realms.emplace(phase::gas,   RefCountedPtr<phase_realm> (new phase_realm()));
  m_realms.emplace(phase::solid, RefCountedPtr<phase_realm> (new phase_realm()));
}

realm::~realm(){

}

void realm::define(const Vector<DisjointBoxLayout>& a_grids,
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
  CH_TIME("realm::define");
  if(m_verbosity > 5){
    pout() << "realm::define" << endl;
  }



  m_ref_ratios = a_ref_rat;
  m_dx = a_dx;
  m_grids = a_grids;
  m_domains = a_domains;

  m_mfis = a_mfis;
  m_finest_level = a_finest_level;
  
  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  m_realms[phase::gas]->define(a_grids, a_domains, a_ref_rat, a_dx, a_finest_level, a_ebghost, a_num_ghost, a_redist_rad,
			       a_centroid_stencil, a_eb_stencil, a_ebcf, ebis_gas);

  m_realms[phase::solid]->define(a_grids, a_domains, a_ref_rat, a_dx, a_finest_level, a_ebghost, a_num_ghost, a_redist_rad,
				 a_centroid_stencil, a_eb_stencil, a_ebcf, ebis_sol);
}

void realm::set_grids(const Vector<DisjointBoxLayout>& a_grids, const int a_finest_level){
  CH_TIME("realm::set_grids");
  if(m_verbosity > 5){
    pout() << "realm::set_grids" << endl;
  }

  for (auto& r : m_realms){
    r.second->set_grids(a_grids, a_finest_level);
  }
}

void realm::regrid_base(const int a_lmin){
  CH_TIME("realm::regrid_base");
  if(m_verbosity > 5){
    pout() << "realm::regrid_base" << endl;
  }
  
  for (auto& r : m_realms){
    r.second->regrid_base(a_lmin);
  }
  this->define_mflevelgrid(a_lmin);
}

void realm::regrid_operators(const int a_lmin, const int a_lmax, const int a_regsize){
  CH_TIME("realm::regrid_operators");
  if(m_verbosity > 5){
    pout() << "realm::regrid_operators" << endl;
  }

  for (auto& r : m_realms){
    r.second->regrid_operators(a_lmin, a_lmax, a_regsize);
  }
}

void realm::define_mflevelgrid(const int a_lmin){
  CH_TIME("realm::define_mflevelgrid");
  if(m_verbosity > 5){
    pout() << "realm::define_mflevelgrid" << endl;
  }

  m_mflg.resize(1 + m_finest_level);

  phase_realm& gas = this->get_realm(phase::gas);
  phase_realm& sol = this->get_realm(phase::solid);

  const RefCountedPtr<EBIndexSpace>& ebis_gas = gas.get_ebis();
  const RefCountedPtr<EBIndexSpace>& ebis_sol = sol.get_ebis();

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBLevelGrid> eblgs;

    if(!ebis_gas.isNull()) eblgs.push_back(*(gas.get_eblg()[lvl]));
    if(!ebis_sol.isNull()) eblgs.push_back(*(sol.get_eblg()[lvl]));

    m_mflg[lvl] = RefCountedPtr<MFLevelGrid> (new MFLevelGrid(m_mfis, eblgs));
  }
}

void realm::register_operator(const std::string a_operator, const phase::which_phase a_phase){
  CH_TIME("realm::register_operator");
  if(m_verbosity > 5){
    pout() << "realm::register_operator" << endl;
  }
  
  m_realms[a_phase]->register_operator(a_operator);
}

bool realm::query_operator(const std::string a_operator, phase::which_phase a_phase){
  CH_TIME("realm::query_operator");
  if(m_verbosity > 5){
    pout() << "realm::query_operator" << endl;
  }
  
  return m_realms[a_phase]->query_operator(a_operator);
}

phase_realm& realm::get_realm(const phase::which_phase a_phase){
  return *m_realms[a_phase];
}

Vector<int>& realm::get_ref_rat() {
  return m_ref_ratios;
}

Vector<Real>& realm::get_dx() {
  return m_dx;
}

Vector<DisjointBoxLayout>& realm::get_grids() {
  return m_grids;
}

Vector<ProblemDomain>& realm::get_domains() {
  return m_domains;
}

Vector<RefCountedPtr<MFLevelGrid> >& realm::get_mflg(){
  return m_mflg;
}

const RefCountedPtr<EBIndexSpace>& realm::get_ebis(const phase::which_phase a_phase) {

  return m_realms[a_phase]->get_ebis();
}

Vector<EBISLayout>& realm::get_ebisl(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_ebisl();
}

Vector<RefCountedPtr<EBLevelGrid> >& realm::get_eblg(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_eblg();
}

Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& realm::get_neighbors(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_neighbors();
}

Vector<RefCountedPtr<LayoutData<VoFIterator> > >& realm::get_vofit(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_vofit();
}

irreg_amr_stencil<centroid_interp>& realm::get_centroid_interp_stencils(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_centroid_interp_stencils();
}

irreg_amr_stencil<eb_centroid_interp>& realm::get_eb_centroid_interp_stencils(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_eb_centroid_interp_stencils();
}

irreg_amr_stencil<noncons_div>& realm::get_noncons_div_stencils(const phase::which_phase a_phase) {
  return m_realms[a_phase]->get_noncons_div_stencils();
}

Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& realm::get_gradsten(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_gradsten();
}

Vector<RefCountedPtr<ebcoarseaverage> >& realm::get_coarave(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_coarave();
}

Vector<RefCountedPtr<EBGhostCloud> >& realm::get_ghostcloud(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_ghostcloud();
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& realm::get_quadcfi(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_quadcfi();
}

Vector<RefCountedPtr<EBQuadCFInterp> >& realm::get_old_quadcfi(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_old_quadcfi();
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& realm::get_fillpatch(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_fillpatch();
}

Vector<RefCountedPtr<EBPWLFineInterp> >& realm::get_eb_pwl_interp(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_eb_pwl_interp();
}

Vector<RefCountedPtr<EBMGInterp> >& realm::get_eb_mg_interp(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_eb_mg_interp();
}

Vector<RefCountedPtr<EBFluxRegister> >&  realm::get_flux_reg(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_flux_reg();
}

Vector<RefCountedPtr<EBLevelRedist> >&  realm::get_level_redist(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_level_redist();
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  realm::get_coar_to_fine_redist(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_coar_to_fine_redist();
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  realm::get_coar_to_coar_redist(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_coar_to_coar_redist();
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  realm::get_fine_to_coar_redist(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_fine_to_coar_redist();
}

Vector<RefCountedPtr<Copier> >& realm::get_copier(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_copier();
}

Vector<RefCountedPtr<Copier> >& realm::get_reverse_copier(const phase::which_phase a_phase){
  return m_realms[a_phase]->get_reverse_copier();
}
