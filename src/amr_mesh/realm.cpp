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
		      const stencil_type a_centroid_stencil,
		      const stencil_type a_eb_stencil,
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
  this->define_masks(a_lmin);
}

void realm::define_masks(const int a_lmin){
  CH_TIME("realm::define_masks");
  if(m_verbosity > 5){
    pout() << "realm::define_masks" << endl;
  }

  // Regrid all masks
  this->define_halo_masks(a_lmin);
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

void realm::define_halo_masks(const int a_lmin){
  CH_TIME("realm::define_halo_masks");
  if(m_verbosity > 5){
    pout() << "realm::define_halo_masks" << endl;
  }

  // Loop through all masks and do something about the halo masks only. 
  for (auto& m : m_masks){

    // Get mask identifier and buffer. 
    const std::string which_mask = m.first.first;
    const int buffer             = m.first.second;

    if(which_mask == s_particle_halo){
      if(buffer < 0) MayDay::Abort("realm::define_halo_masks -- cannot have buffer < 0!");

      AMRMask& mask = m.second;

      mask.resize(1 + m_finest_level);

      for (int lvl = 0; lvl < m_finest_level; lvl++){
	const DisjointBoxLayout& gridsCoar = m_grids[lvl];
	const DisjointBoxLayout& gridsFine = m_grids[lvl+1];

	const ProblemDomain& domainCoar = m_domains[lvl];
	const ProblemDomain& domainFine = m_domains[lvl+1];

	const int ncomp = 1;

	mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >(new LevelData<BaseFab<bool> >(gridsCoar, ncomp, IntVect::Zero));

	this->define_halo_mask(*mask[lvl],
			       domainCoar,
			       domainFine,
			       gridsCoar,
			       gridsFine,
			       buffer,
			       m_ref_ratios[lvl]);
      }
    }
  }
}

void realm::define_halo_mask(LevelData<BaseFab<bool> >& a_coarMask,
			     const ProblemDomain&       a_domainCoar,
			     const ProblemDomain&       a_domainFine,
			     const DisjointBoxLayout&   a_gridsCoar,
			     const DisjointBoxLayout&   a_gridsFine,
			     const int                  a_buffer,
			     const int                  a_ref_rat){
  CH_TIME("realm::define_halo_masks");
  if(m_verbosity > 5){
    pout() << "realm::define_halo_masks" << endl;
  }

  // TLDR: This function creates a GLOBAL view of the CFIVS, a_buffer cells wide on the coarse side of the refinement boundary. If this is not performant,
  //       there are probably better ways to do this. 
  
  const int comp  = 0;
  const int ncomp = 1;

  // First, reset the mask.
  for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
    a_coarMask[dit()].setVal(false);
  }

  // Ok, we need a particle halo. 
  if(a_buffer > 0){

    // Coarsen the fine grid and make a mask on the coarsened fine grid
    DisjointBoxLayout gridsCoFi;
    coarsen(gridsCoFi, a_gridsFine, a_ref_rat);

    IntVectSet halo;                 // Global CFIVS - should find a better way to do this.
    NeighborIterator nit(gridsCoFi); // Neighbor iterator

    // Go through the cofi grid and set the halo to true
    for (DataIterator dit(gridsCoFi); dit.ok(); ++dit){
      const Box coFiBox = gridsCoFi[dit()];

      // Make IntVect set consisting of only ghost cells (a_buffer) for each box. 
      Box grownBox = grow(coFiBox, a_buffer);
      grownBox    &= a_domainCoar;

      // Subtract non-ghosted box. 
      IntVectSet myHalo(grownBox);
      myHalo -= coFiBox;

      // Subtract non-ghosted neighboring boxes. 
      for (nit.begin(dit()); nit.ok(); ++nit){
	const Box neighborBox = gridsCoFi[nit()];
	myHalo -= neighborBox;
      }

      // Add to halo. 
      halo |= myHalo;
    }

    // TLDR: In the above, we found the coarse-grid cells surrounding the fine level, viewed from the fine grids.
    //       Below, we create that view from the coarse grid. We use a BoxLayoutData<FArrayBox> on the coarsened fine grid,
    //       whose "ghost cells" can added to the _actual_ coarse grid. We then loop through those cells and set the mask. 

    // Get some box definitions. 
    Vector<Box> coFiBoxes = gridsCoFi.boxArray();
    Vector<int> coFiProcs = gridsCoFi.procIDs();

    Vector<Box> coarBoxes = a_gridsCoar.boxArray();
    Vector<Box> coarPrioc = a_gridsCoar.boxArray();

    // Grow the coarsened fine boxes
    for (int i = 0; i < coFiBoxes.size(); i++){
      coFiBoxes[i].grow(a_buffer);
      coFiBoxes[i] &= a_domainCoar;
    }

    // Define boxlayouts
    BoxLayout blCoFi(coFiBoxes, coFiProcs);
    BoxLayout blCoar(a_gridsCoar.boxArray(), a_gridsCoar.procIDs());

    BoxLayoutData<FArrayBox> coFiMask(blCoFi, ncomp);
    BoxLayoutData<FArrayBox> coarMask(blCoar, ncomp);
  
    // Reset masks
    for (DataIterator dit(blCoFi); dit.ok(); ++dit){
      coFiMask[dit()].setVal(0.0);
    }
    for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
      coarMask[dit()].setVal(0.0);
    }

    // Run through the halo and set the halo cells to 1 on the coarsened fine grids (these are essentially "ghost cells")
    for (DataIterator dit(blCoFi); dit.ok(); ++dit){
      IntVectSet curHalo = halo & blCoFi[dit()];
      for (IVSIterator ivsit(curHalo); ivsit.ok(); ++ivsit){
	coFiMask[dit()](ivsit(), comp) = 1.0;
      }
    }

    // Add the result to the coarse grid. 
    const Interval interv(0,0);
    coFiMask.addTo(interv, coarMask, interv, a_domainCoar.domainBox());

    // Run through the grids and make the boolean mask
    for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
      const Box box = a_gridsCoar[dit()];

      BaseFab<bool>& mask     = a_coarMask[dit()];
      const FArrayBox& blMask = coarMask[dit()];

      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();
	if(blMask(iv, comp) > 0.0) mask(iv, comp) = true;
      }
    }
  }
}

void realm::register_operator(const std::string a_operator, const phase::which_phase a_phase){
  CH_TIME("realm::register_operator(operator, phase)");
  if(m_verbosity > 5){
    pout() << "realm::register_operator(operator, phase)" << endl;
  }
  
  m_realms[a_phase]->register_operator(a_operator);
}

bool realm::query_operator(const std::string a_operator, const phase::which_phase a_phase) {
  CH_TIME("realm::query_operator");
  if(m_verbosity > 5){
    pout() << "realm::query_operator" << endl;
  }
  
  return m_realms[a_phase]->query_operator(a_operator);
}


void realm::register_levelset(const phase::which_phase a_phase){
  CH_TIME("realm::register_levelset");
  if(m_verbosity > 5){
    pout() << "realm::register_levelset" << endl;
  }

  //  m_masks.emplace(a_std::pair<string, int>(a_mask, a_buffer), AMRMask());
}
								

void realm::query_levelset(const phase::which_phase a_phase){
  CH_TIME("realm::query_levelset");
  if(m_verbosity > 5){
    pout() << "realm::query_levelset" << endl;
  }
}

void realm::register_mask(const std::string a_mask, const int a_buffer){
  CH_TIME("realm::register_mask(mask, buffer)");
  if(m_verbosity > 5){
    pout() << "realm::register_mask(mask, buffer)" << endl;
  }

  m_masks.emplace(std::pair<string, int>(a_mask, a_buffer), AMRMask());
}

bool realm::query_mask(const std::string a_mask, const int a_buffer) const {
  CH_TIME("realm::query_mask(mask, buffer)");
  if(m_verbosity > 5){
    pout() << "realm::query_mask(mask, buffer)" << endl;
  }

  bool ret;

  if(m_masks.count(std::pair<string, int>(a_mask, a_buffer) ) > 0){
    ret = true;
  }
  else{
    ret = false;
  }

  return ret;
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

AMRMask& realm::get_mask(const std::string a_mask, const int a_buffer){

  if(!this->query_mask(a_mask, a_buffer)){
    std::string str = "realm::get_mask - could not find mask '" + a_mask + "'";
    MayDay::Abort(str.c_str());
  }

  return m_masks.at(std::pair<std::string, int>(a_mask, a_buffer));
}
