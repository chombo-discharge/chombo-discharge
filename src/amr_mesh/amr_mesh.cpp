/*!
  @file amr_mesh.cpp
  @brief Implementation of amr_mesh.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "amr_mesh.H"
#include "mfalias.H"

#include <BRMeshRefine.H>
#include <EBEllipticLoadBalance.H>
#include <EBLevelDataOps.H>
#include <MFLevelDataOps.H>
#include <EBArith.H>

amr_mesh::amr_mesh(){

  // Default stuff will crash.
  this->set_verbosity(10);
  this->set_coarsest_num_cells(IntVect::Zero);
  this->set_max_amr_depth(-1);
  this->set_refinement_ratio(-2);
  this->set_blocking_factor(-8);
  this->set_max_box_size(-32);
  this->set_buffer_size(0);
  this->set_ebcf(true);
  this->set_fill_ratio(-1.0);
  this->set_redist_rad(-1);
  this->set_num_ghost(-2);
  this->set_eb_ghost(-4);
  this->set_irreg_sten_order(-1);
  this->set_irreg_sten_radius(-1);
}

amr_mesh::~amr_mesh(){
  
}

template<typename T> void amr_mesh::deallocate(Vector<T*>& a_data){
  CH_TIME("amr_mesh::deallocate");
  if(m_verbosity > 5){
    pout() << "amr_mesh::deallocate" << endl;
  }
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    delete a_data[lvl];
  }
}

template<typename T> void amr_mesh::alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T> >& a_data){
  CH_TIME("amr_mesh::alias");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias" << endl;
  }

  a_alias.resize(a_data.size());
  
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    a_alias[lvl] = &(*a_data[lvl]);
  }
}

void amr_mesh::allocate(EBAMRCellData& a_data, phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(cell)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(cell)" << endl;
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    EBCellFactory fact(m_ebisl[a_phase][lvl]);

    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB> >
      (new LevelData<EBCellFAB>(m_grids[lvl], a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::allocate(EBAMRFluxData& a_data, phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(flux)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(flux)" << endl;
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    EBFluxFactory fact(m_ebisl[a_phase][lvl]);

    a_data[lvl] = RefCountedPtr<LevelData<EBFluxFAB> >
      (new LevelData<EBFluxFAB>(m_grids[lvl], a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::allocate(EBAMRIVData& a_data, phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(baseiv)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(baseiv)" << endl;
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    LayoutData<IntVectSet> irreg_sets(m_grids[lvl]);
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      Box box = m_grids[lvl].get(dit());
      box.grow(ghost);
      box &= m_domains[lvl];

      irreg_sets[dit()] = m_ebisl[a_phase][lvl][dit()].getIrregIVS(box);
    }

    BaseIVFactory<Real> fact(m_ebisl[a_phase][lvl], irreg_sets);

    a_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> >(m_grids[lvl], a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::allocate(MFAMRCellData& a_data, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mf cell)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mf cell)" << endl;
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = a_ncomp;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(phase::num_phases);
    Vector<int>        comps(phase::num_phases, a_ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    
    MFCellFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFCellFAB> >
      (new LevelData<MFCellFAB>(m_grids[lvl], ignored, ghost*IntVect::Unit, factory));

    MFLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
  
}

void amr_mesh::allocate(MFAMRFluxData& a_data, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mf flux)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mf flux)" << endl;
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = 99;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(phase::num_phases);
    Vector<int>        comps(phase::num_phases, a_ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    
    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFFluxFAB> >
      (new LevelData<MFFluxFAB>(m_grids[lvl], ignored, ghost*IntVect::Unit, factory));
  }
}
  
void amr_mesh::allocate(MFAMRIVData& a_data, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mf baseivfab)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mf baseivfab)" << endl;
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = 99;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(phase::num_phases);
    Vector<int>        comps(phase::num_phases, a_ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    
    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFBaseIVFAB> >
      (new LevelData<MFBaseIVFAB>(m_grids[lvl], ignored, ghost*IntVect::Unit, factory));
  }
}

void amr_mesh::allocate_interface(EBAMRIVData& a_data, phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate_interface(baseiv)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate_interface(baseiv)" << endl;
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const IntVectSet interface_ivs = m_mfis->interface_region(m_domains[lvl]);
    LayoutData<IntVectSet> irreg_sets(m_grids[lvl]);
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      Box box = m_grids[lvl].get(dit());
      box.grow(ghost);
      box &= m_domains[lvl];

      irreg_sets[dit()] = interface_ivs & box; 
    }

    BaseIVFactory<Real> fact(m_ebisl[a_phase][lvl], irreg_sets);

    a_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> >(m_grids[lvl], a_ncomp, ghost*IntVect::Unit, fact));
  }
}

void amr_mesh::set_mfis(const RefCountedPtr<mfis>& a_mfis){
  CH_TIME("amr_mesh::set_mfis");
  if(m_verbosity > 3){
    pout() << "amr_mesh::set_mfis" << endl;
  }

  m_mfis = a_mfis;
}

void amr_mesh::build_domains(){

  const int nlevels = 1 + m_max_amr_depth;
  
  m_domains.resize(nlevels);
  m_dx.resize(nlevels);
  m_ref_ratios.resize(nlevels, m_ref_ratio);
  m_grids.resize(nlevels);

  m_eblg.resize(phase::num_phases);
  m_ebisl.resize(phase::num_phases);
  m_coarave.resize(phase::num_phases);
  m_quadcfi.resize(phase::num_phases);
  m_old_quadcfi.resize(phase::num_phases);
  m_flux_reg.resize(phase::num_phases);
  m_level_redist.resize(phase::num_phases);
  m_centroid_interp.resize(phase::num_phases);
  m_eb_centroid_interp.resize(phase::num_phases);
  m_coar_to_fine_redist.resize(phase::num_phases);
  m_coar_to_coar_redist.resize(phase::num_phases);
  m_fine_to_coar_redist.resize(phase::num_phases);

  m_eblg[phase::gas].resize(nlevels);
  m_ebisl[phase::gas].resize(nlevels);
  m_coarave[phase::gas].resize(nlevels);
  m_quadcfi[phase::gas].resize(nlevels);
  m_flux_reg[phase::gas].resize(nlevels);
  m_old_quadcfi[phase::gas].resize(nlevels);
  m_level_redist[phase::gas].resize(nlevels);
  m_coar_to_fine_redist[phase::gas].resize(nlevels);
  m_coar_to_coar_redist[phase::gas].resize(nlevels);
  m_fine_to_coar_redist[phase::gas].resize(nlevels);

  m_eblg[phase::solid].resize(nlevels);
  m_ebisl[phase::solid].resize(nlevels);
  m_coarave[phase::solid].resize(nlevels);
  m_quadcfi[phase::solid].resize(nlevels);
  m_flux_reg[phase::solid].resize(nlevels);
  m_old_quadcfi[phase::solid].resize(nlevels);
  m_level_redist[phase::solid].resize(nlevels);
  m_coar_to_fine_redist[phase::solid].resize(nlevels);
  m_coar_to_coar_redist[phase::solid].resize(nlevels);
  m_fine_to_coar_redist[phase::solid].resize(nlevels);


  m_dx[0] = (m_physdom->get_prob_hi()[0] - m_physdom->get_prob_lo()[0])/m_num_cells[0];
  m_domains[0] = ProblemDomain(IntVect::Zero, m_num_cells - IntVect::Unit);

  for (int lvl = 1; lvl <= m_max_amr_depth; lvl++){
    m_dx[lvl]      = m_dx[lvl-1]/m_ref_ratio;
    m_domains[lvl] = m_domains[lvl-1];
    m_domains[lvl].refine(m_ref_ratio);
  }
}

void amr_mesh::regrid(const Vector<IntVectSet>& a_tags){
  CH_TIME("amr_mesh::regrid");
  if(m_verbosity > 3){
    pout() << "amr_mesh::regrid" << endl;
  }

  if(a_tags.size() > 0){ // Not regridding if I don't get tags

    Vector<IntVectSet> tags = a_tags; // build_grids destroys tags, so copy them
    this->build_grids(tags);  

    this->define_eblevelgrid(); // Define EBLevelGrid objects on both phases
    this->define_eb_coar_ave(); // Define EBCoarseAverage on both phases
    this->define_eb_quad_cfi(); // Define nwoebquadcfinterp on both phases. This crashes for ref_rat = 4
    this->define_flux_reg();    // Define flux register (phase::gas only)
    this->define_redist_oper(); // Define redistribution (phase::gas only)
    this->define_irreg_sten();  // Define irregular stencils
  }
}

void amr_mesh::build_grids(Vector<IntVectSet>& a_tags){
  CH_TIME("amr_mesh::build_grids");
  if(m_verbosity > 5){
    pout() << "amr_mesh::build_grids" << endl;
  }

  const int base      = 0;                       // I don't want this level to change
  const int top_level = a_tags.size() - 1;       // top_level is the finest level where we have tags
  Vector<Vector<Box> > new_boxes(1 + top_level); // New boxes to be load balance
  Vector<Vector<Box> > old_boxes(1 + top_level); // Old grids. 

  for (int lvl = 0; lvl <= top_level; lvl++){
    old_boxes[lvl].push_back(m_domains[lvl].domainBox()); // Create old grids from scratch 
  }

  // Berger-Rigoutsos grid generation
  BRMeshRefine mesh_refine(m_domains[0], m_ref_ratios, m_fill_ratio, m_blocking_factor, m_buffer_size, m_max_box_size);
  int new_finest_level = mesh_refine.regrid(new_boxes, a_tags, base, top_level, old_boxes);

  m_finest_level = Min(new_finest_level, m_max_amr_depth); // Don't exceed m_max_amr_depth


  for (int lvl = 0; lvl <= m_finest_level; lvl++){   // Generate DisjointBoxLayouts on each level
    Vector<int> proc_assign;
    this->load_balance(proc_assign, new_boxes[lvl], lvl);             // Load balance, assign boxes to grids
    m_grids[lvl].define(new_boxes[lvl], proc_assign, m_domains[lvl]); // Define grids
  }
}

void amr_mesh::load_balance(Vector<int>& a_proc_assign, Vector<Box>& a_boxes, const int a_lvl){
  CH_TIME("amr_mesh::load_balance");
  if(m_verbosity > 5){
    pout() << "amr_mesh::load_balance" << endl;
  }

  EBEllipticLoadBalance(a_proc_assign, a_boxes, m_domains[a_lvl], false, m_mfis->get_ebis(phase::gas)); // Loads for each box
}

void amr_mesh::compute_gradient(EBAMRCellData& a_gradient, EBAMRCellData& a_phi){
  CH_TIME("amr_mesh::compute_gradient");
  if(m_verbosity > 5){
    pout() << "amr_mesh::compute_gradient" << endl;
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    CH_assert(a_phi[lvl]->nComp()      == 1);
    CH_assert(a_gradient[lvl]->nComp() == SpaceDim);
    
    const DisjointBoxLayout& dbl = m_grids[lvl]; // Doing this since I assume everything is defined over m_grids
    
    LayoutData<IntVectSet> cfivs;
    EBArith::defineCFIVS(cfivs, dbl, m_domains[lvl]);

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& grad        = (*a_gradient[lvl])[dit()];
      const EBCellFAB& phi   = (*a_phi[lvl])[dit()];
      const EBISBox& ebisbox = phi.getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(dbl.get(dit()));

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();

	for (int dir = 0; dir < SpaceDim; dir++){
	  
	  grad(vof, dir) = 0.;

	  VoFStencil sten;
	  EBArith::getFirstDerivStencil(sten, vof, ebisbox, dir, m_dx[lvl], &cfivs[dit()], 0);
	  for (int i = 0; i < sten.size(); i++){
	    const VolIndex& ivof = sten.vof(i);
	    const Real& iweight  = sten.weight(i);
	    
	    grad(vof, dir) += phi(ivof, 0)*iweight;
	  }
	}
      }
    }
  }
}

void amr_mesh::compute_gradient(MFAMRCellData& a_gradient, MFAMRCellData& a_phi){
  CH_TIME("amr_mesh::compute_gradient(mf)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::compute_gradient(mf)" << endl;
  }

  for (int iphase = 0; iphase < m_mfis->num_phases(); iphase++){
    EBAMRCellData alias_grad(1 + m_finest_level);
    EBAMRCellData alias_phi(1 + m_finest_level);

    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      alias_grad[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      alias_phi[lvl]  = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
      mfalias::aliasMF(*alias_grad[lvl], iphase, *a_gradient[lvl]);
      mfalias::aliasMF(*alias_phi[lvl],  iphase, *a_phi[lvl]);
    }

    this->compute_gradient(alias_grad, alias_phi);
  }
}

void amr_mesh::define_eblevelgrid(){
  CH_TIME("amr_mesh::define_eblevelgrid");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eblevelgrid" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(!ebis_sol.isNull()){
      m_eblg[phase::gas][lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_ebghost, ebis_gas));
      m_eblg[phase::gas][lvl]->setMaxCoarseningRatio(m_ref_ratios[lvl], ebis_gas);
      m_ebisl[phase::gas][lvl] = m_eblg[phase::gas][lvl]->getEBISL();
    }
    if(!ebis_sol.isNull()){
      m_eblg[phase::solid][lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_ebghost, ebis_sol));
      m_eblg[phase::solid][lvl]->setMaxCoarseningRatio(m_ref_ratios[lvl], ebis_sol);
      m_ebisl[phase::solid][lvl] = m_eblg[phase::solid][lvl]->getEBISL();
    }
  }
}

void amr_mesh::define_eb_coar_ave(){
  CH_TIME("amr_mesh::define_eb_coar_ave");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eb_coar_ave" << endl;
  }

  const int comps = SpaceDim;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	m_coarave[phase::gas][lvl] = RefCountedPtr<EBCoarseAverage> (new EBCoarseAverage(m_grids[lvl],
											 m_grids[lvl-1],
											 m_ebisl[phase::gas][lvl],
											 m_ebisl[phase::gas][lvl-1],
											 m_domains[lvl-1],
											 m_ref_ratios[lvl-1],
											 comps,
											 ebis_gas));
      }
      if(!ebis_sol.isNull()){
	m_coarave[phase::solid][lvl] = RefCountedPtr<EBCoarseAverage> (new EBCoarseAverage(m_grids[lvl],
											   m_grids[lvl-1],
											   m_ebisl[phase::solid][lvl],
											   m_ebisl[phase::solid][lvl-1],
											   m_domains[lvl-1],
											   m_ref_ratios[lvl-1],
											   comps,
											   ebis_sol));
      }
    }
  }
}

void amr_mesh::define_eb_quad_cfi(){
  CH_TIME("amr_mesh::define_eb_quad_cfi");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eb_quad_cfi" << endl;
  }

  const int comps = SpaceDim;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	const LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::gas][lvl]->getCFIVS());
	m_quadcfi[phase::gas][lvl] = RefCountedPtr<nwoebquadcfinterp> (new nwoebquadcfinterp(m_grids[lvl],
											     m_grids[lvl-1],
											     m_ebisl[phase::gas][lvl],
											     m_ebisl[phase::gas][lvl-1],
											     m_domains[lvl-1],
											     m_ref_ratios[lvl-1],
											     comps,
											     m_dx[lvl],
											     m_num_ghost,
											     cfivs,
											     ebis_gas));

	m_old_quadcfi[phase::gas][lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(m_grids[lvl],
											   m_grids[lvl-1],
											   m_ebisl[phase::gas][lvl],
											   m_ebisl[phase::gas][lvl-1],
											   m_domains[lvl-1],
											   m_ref_ratios[lvl-1],
											   1,
											   cfivs,
											   ebis_gas));
											   
      }
      if(!ebis_sol.isNull()){
	LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::solid][lvl]->getCFIVS());
	m_quadcfi[phase::solid][lvl] = RefCountedPtr<nwoebquadcfinterp> (new nwoebquadcfinterp(m_grids[lvl],
											       m_grids[lvl-1],
											       m_ebisl[phase::solid][lvl],
											       m_ebisl[phase::solid][lvl-1],
											       m_domains[lvl-1],
											       m_ref_ratios[lvl-1],
											       comps,
											       m_dx[lvl],
											       m_num_ghost,
											       cfivs,
											       ebis_sol));
	
	m_old_quadcfi[phase::solid][lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(m_grids[lvl],
											   m_grids[lvl-1],
											   m_ebisl[phase::solid][lvl],
											   m_ebisl[phase::solid][lvl-1],
											   m_domains[lvl-1],
											   m_ref_ratios[lvl-1],
											   1,
											   cfivs,
											   ebis_gas));
      }
    }
  }
}

void amr_mesh::define_flux_reg(){
  CH_TIME("amr_mesh::define_flux_reg");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_flux_reg" << endl;
  }

  const int comps = 1;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_fine = lvl < m_finest_level;

    if(!ebis_gas.isNull()){
      if(has_fine){
	if(!ebis_gas.isNull()){
	  m_flux_reg[phase::solid][lvl] = RefCountedPtr<EBFastFR> (new EBFastFR(*m_eblg[phase::gas][lvl+1],
										*m_eblg[phase::gas][lvl],
										m_ref_ratios[lvl],
										comps,
										!m_ebcf));
	}
      }
    }
  }
}

void amr_mesh::define_redist_oper(){
  CH_TIME("amr_mesh::define_redist_oper");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_redist_oper" << endl;
  }

  const int comps = 1;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl > m_finest_level;

    if(!ebis_gas.isNull()){
      m_level_redist[phase::gas][lvl] = RefCountedPtr<EBLevelRedist> (new EBLevelRedist(m_grids[lvl],
											m_ebisl[phase::gas][lvl],
											m_domains[lvl],
											comps,
											m_redist_rad));

    
      if(m_ebcf){
	if(has_coar){
	  m_fine_to_coar_redist[phase::gas][lvl] = RefCountedPtr<EBFineToCoarRedist> (new EBFineToCoarRedist());
	  m_fine_to_coar_redist[phase::gas][lvl]->define(m_grids[lvl],
							 m_grids[lvl-1],
							 m_ebisl[phase::gas][lvl],
							 m_ebisl[phase::gas][lvl-1],
							 m_domains[lvl-1].domainBox(),
							 m_ref_ratios[lvl-1],
							 comps,
							 m_redist_rad,
							 ebis_gas);

	  // Set register to zero
	  m_fine_to_coar_redist[phase::gas][lvl]->setToZero();
	}

	if(has_fine){

	  m_coar_to_fine_redist[phase::gas][lvl] = RefCountedPtr<EBCoarToFineRedist> (new EBCoarToFineRedist());
	  m_coar_to_fine_redist[phase::gas][lvl]->define(m_grids[lvl+1],
							 m_grids[lvl],
							 m_ebisl[phase::gas][lvl],
							 m_domains[lvl].domainBox(),
							 m_ref_ratios[lvl],
							 comps,
							 m_redist_rad,
							 ebis_gas);

	  // Coarse to coarse redistribution
	  m_coar_to_coar_redist[phase::gas][lvl] = RefCountedPtr<EBCoarToCoarRedist> (new EBCoarToCoarRedist());
	  m_coar_to_coar_redist[phase::gas][lvl]->define(*m_eblg[phase::gas][lvl+1],
							 *m_eblg[phase::gas][lvl],
							 m_ref_ratios[lvl],
							 comps,
							 m_redist_rad);


	  // Set registers to zero
	  m_coar_to_fine_redist[phase::gas][lvl]->setToZero();
	  m_coar_to_coar_redist[phase::gas][lvl]->setToZero();
	}
      }
    }
  }
}

void amr_mesh::define_irreg_sten(){
  CH_TIME("amr_mesh::define_irreg_sten");
  if(m_verbosity > 3){
    pout() << "amr_mesh::define_irreg_sten" << endl;
  }

  const int order = m_irreg_sten_order;
  const int rad   = m_irreg_sten_radius;


  m_centroid_interp[phase::gas] = RefCountedPtr<irreg_amr_stencil<centroid_interp> >
    (new irreg_amr_stencil<centroid_interp>(m_grids, m_ebisl[phase::gas], m_domains, m_dx, m_finest_level, order, rad));
  m_eb_centroid_interp[phase::gas] = RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> >
    (new irreg_amr_stencil<eb_centroid_interp>(m_grids, m_ebisl[phase::gas], m_domains, m_dx, m_finest_level, order, rad));

  
  m_centroid_interp[phase::solid] = RefCountedPtr<irreg_amr_stencil<centroid_interp> >
    (new irreg_amr_stencil<centroid_interp>(m_grids, m_ebisl[phase::solid], m_domains, m_dx, m_finest_level, order, rad));
  m_eb_centroid_interp[phase::solid] = RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> >
    (new irreg_amr_stencil<eb_centroid_interp>(m_grids, m_ebisl[phase::solid], m_domains, m_dx, m_finest_level, order, rad));
}

void amr_mesh::average_down(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(ebcell)" << endl;
  }

  //
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

    a_data[lvl]->exchange(interv);
  }
}

void amr_mesh::average_down(MFAMRCellData& a_data){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(mf)" << endl;
  }
  
  EBAMRCellData alias_g(1 + m_finest_level);
  EBAMRCellData alias_s(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
  }

  this->average_down(alias_g, phase::gas);
  this->average_down(alias_s, phase::solid);
}

void amr_mesh::average_down(EBAMRFluxData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(face)" << endl;
  }

  //
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

    a_data[lvl]->exchange(interv);
  }
}

void amr_mesh::average_down(EBAMRIVData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(iv)" << endl;
  }

  //
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

    a_data[lvl]->exchange(interv);
  }
}

void amr_mesh::interp_ghost(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    m_quadcfi[a_phase][lvl]->coarseFineInterp(*a_data[lvl], *a_data[lvl-1], 0, 0, ncomps);

    a_data[lvl]->exchange(interv);
  }
}

void amr_mesh::interp_ghost(MFAMRCellData& a_data){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(mf)" << endl;
  }
  
  EBAMRCellData alias_g(1 + m_finest_level);
  EBAMRCellData alias_s(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
  }

  this->interp_ghost(alias_g, phase::gas);
  this->interp_ghost(alias_s, phase::solid);
}

void amr_mesh::set_verbosity(const int a_verbosity){
  CH_TIME("amr_mesh::set_verbosity");

  m_verbosity = a_verbosity;
}

void amr_mesh::set_coarsest_num_cells(const IntVect a_num_cells){
  m_num_cells = a_num_cells;
}

void amr_mesh::set_max_amr_depth(const int a_max_amr_depth){
  m_max_amr_depth = a_max_amr_depth;
}

void amr_mesh::set_ebcf(const bool a_ebcf){
  m_ebcf = a_ebcf;
}

void amr_mesh::set_refinement_ratio(const int a_refinement_ratio){
  CH_TIME("amr_mesh::set_refinement_ratio");
  m_ref_ratio = a_refinement_ratio;
}

void amr_mesh::set_fill_ratio(const Real a_fill_ratio){
  m_fill_ratio = a_fill_ratio;
}

void amr_mesh::set_max_box_size(const int a_max_box_size){
  CH_TIME("amr_mesh::set_max_box_size");
  m_max_box_size = a_max_box_size;
}

void amr_mesh::set_buffer_size(const int a_buffer_size){
  m_buffer_size = a_buffer_size;
}

void amr_mesh::set_blocking_factor(const int a_blocking_factor){
  CH_TIME("amr_mesh::set_blocking_factor");
  m_blocking_factor = a_blocking_factor;
}

void amr_mesh::set_eb_ghost(const int a_ebghost){
  CH_TIME("amr_mesh:.set_eb_ghost");
  m_ebghost = a_ebghost;
}

void amr_mesh::set_num_ghost(const int a_num_ghost){
  CH_TIME("amr_mesh::set_num_ghost");
  m_num_ghost = a_num_ghost;
}

void amr_mesh::set_redist_rad(const int a_redist_rad){
  CH_TIME("amr_mesh::set_redist_rads");
  m_redist_rad = a_redist_rad;
}

void amr_mesh::set_irreg_sten_order(const int a_irreg_sten_order){
  CH_TIME("amr_mesh::irreg_sten_order");
  m_irreg_sten_order = a_irreg_sten_order;
}

void amr_mesh::set_irreg_sten_radius(const int a_irreg_sten_radius){
  CH_TIME("amr_mesh::irreg_sten_radius");
  m_irreg_sten_radius = a_irreg_sten_radius;
}

void amr_mesh::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  m_physdom = a_physdom;
}

void amr_mesh::sanity_check(){
  CH_TIME("amr_mesh::sanity_check");
  if(m_verbosity > 1){
    pout() << "amr_mesh::sanity_check" << endl;
  }

  CH_assert(m_max_amr_depth >= 0);
  CH_assert(m_ref_ratio == 2 || m_ref_ratio == 4);
  CH_assert(m_blocking_factor >= 4 && m_blocking_factor % m_ref_ratio == 0);
  CH_assert(m_max_box_size >= 8 && m_max_box_size % m_blocking_factor == 0);
  CH_assert(m_fill_ratio > 0. && m_fill_ratio <= 1.0);
  CH_assert(m_buffer_size > 0);
  CH_assert(m_irreg_sten_order == 1 ||  m_irreg_sten_order == 2);
  CH_assert(m_irreg_sten_radius == 1 || m_irreg_sten_radius == 2);

  // Make sure that the isotropic constraint isn't violated
  const RealVect realbox = m_physdom->get_prob_hi() - m_physdom->get_prob_lo();
  for (int dir = 0; dir < SpaceDim - 1; dir++){
    CH_assert(realbox[dir]/m_num_cells[dir] == realbox[dir+1]/m_num_cells[dir+1]);
  }
}

int amr_mesh::get_finest_level(){
  return m_finest_level;
}

int amr_mesh::get_max_amr_depth(){
  return m_max_amr_depth;
}

int amr_mesh::get_refinement_ratio(){
  return m_ref_ratio;
}

int amr_mesh::get_num_ghost(){
  return m_num_ghost;
}

int amr_mesh::get_blocking_factor(){
  return m_blocking_factor;
}

int amr_mesh::get_max_box_size(){
  return m_max_box_size;
}

ProblemDomain amr_mesh::get_finest_domain(){
  return m_domains[m_max_amr_depth];
}

Real amr_mesh::get_finest_dx(){
  return m_dx[m_max_amr_depth];
}

Vector<Real>& amr_mesh::get_dx(){
  return m_dx;
}

Vector<int>& amr_mesh::get_ref_rat(){
  return m_ref_ratios;
}

Vector<DisjointBoxLayout>& amr_mesh::get_grids(){
  return m_grids;
}

Vector<ProblemDomain>& amr_mesh::get_domains(){
  return m_domains;
}

Vector<EBISLayout>& amr_mesh::get_ebisl(phase::which_phase a_phase){
  return m_ebisl[a_phase];
}

Vector<RefCountedPtr<EBLevelGrid> >& amr_mesh::get_eblg(phase::which_phase a_phase){
  return m_eblg[a_phase];
}

Vector<RefCountedPtr<EBCoarseAverage> >& amr_mesh::get_coarave(phase::which_phase a_phase){
  return m_coarave[a_phase];
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& amr_mesh::get_quadcfi(phase::which_phase a_phase){
  return m_quadcfi[a_phase];
}

Vector<RefCountedPtr<EBQuadCFInterp> >& amr_mesh::get_old_quadcfi(phase::which_phase a_phase){
  return m_old_quadcfi[a_phase];
}

Vector<RefCountedPtr<EBFastFR> >&  amr_mesh::get_flux_reg(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase. 
  return m_flux_reg[a_phase];
}

Vector<RefCountedPtr<EBLevelRedist> >& amr_mesh::get_level_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase. 
  return m_level_redist[a_phase];
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  amr_mesh::get_coar_to_fine_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase.

  return m_coar_to_fine_redist[a_phase];;
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  amr_mesh::get_coar_to_coar_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase.

  return m_coar_to_coar_redist[a_phase];;
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  amr_mesh::get_fine_to_coar_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase.

  return m_fine_to_coar_redist[a_phase];;
}

irreg_amr_stencil<centroid_interp>& amr_mesh::get_centroid_interp_stencils(phase::which_phase a_phase){
  return *m_centroid_interp[a_phase];
}

irreg_amr_stencil<eb_centroid_interp>& amr_mesh::get_eb_centroid_interp_stencils(phase::which_phase a_phase){
  return *m_eb_centroid_interp[a_phase];
}
