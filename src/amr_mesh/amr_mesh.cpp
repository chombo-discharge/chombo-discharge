/*!
  @file amr_mesh.cpp
  @brief Implementation of amr_mesh.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "amr_mesh.H"
#include "mfalias.H"
#include "load_balance.H"

#include <BRMeshRefine.H>
#include <EBEllipticLoadBalance.H>
#include <EBLevelDataOps.H>
#include <MFLevelDataOps.H>
#include <EBArith.H>
#include <ParmParse.H>

#define advect_level_on 1

amr_mesh::amr_mesh(){

  this->set_verbosity(10);
  this->set_coarsest_num_cells(128*IntVect::Unit);
  this->set_max_amr_depth(0);
  this->set_max_simulation_depth(0);
  this->set_refinement_ratio(2);
#if CH_SPACEDIM == 2
  this->set_blocking_factor(8);
#elif CH_SPACEDIM==3
  this->set_blocking_factor(16);
#endif
  this->set_max_box_size(32);
  this->set_max_ebis_box_size(m_max_box_size); 
  this->set_buffer_size(1);
  this->set_irreg_growth(2);
  this->set_ebcf(false);
  this->set_fill_ratio(1.0);
  this->set_redist_rad(1);
  this->set_num_ghost(3);
  this->set_eb_ghost(4);
  this->set_irreg_sten_order(1);
  this->set_irreg_sten_radius(1);
  this->set_balance(load_balance::knapsack);
  this->set_irreg_sten_type(stencil_type::linear);
  this->set_ghost_interpolation(ghost_interpolation::pwl);
}

amr_mesh::~amr_mesh(){
  
}

void amr_mesh::alias(EBAMRCellData&           a_data,
		     const phase::which_phase a_phase,
		     const MFAMRCellData&     a_mfdata,
		     const int                a_finest_level){
  CH_TIME("amr_mesh::alias(hardcap)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias(hardcap)" << endl;
  }

  for (int lvl = 0; lvl <= a_finest_level; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void amr_mesh::alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata){
  CH_TIME("amr_mesh::alias");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias" << endl;
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void amr_mesh::allocate_ptr(EBAMRCellData& a_data){
  CH_TIME("amr_mesh::allocate_ptr");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate_ptr" << endl;
  }

  a_data.resize(1 + m_finest_level);
  
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB> >(new LevelData<EBCellFAB>());
  }
}

void amr_mesh::allocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
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

void amr_mesh::allocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
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

void amr_mesh::allocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
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
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){ // Single phase
      ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    }
    
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
  const int ignored = a_ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){
      ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    }
    
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
  const int ignored = a_ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){
      ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    }
    
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

  MayDay::Abort("amr_mesh::allocate_interface - this is not the way to do it. Please use MFLevelGrid");
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
  if(m_verbosity > 5){
    pout() << "amr_mesh::set_mfis" << endl;
  }

  m_mfis = a_mfis;
}

void amr_mesh::set_balance(load_balance::which_balance a_load){
  CH_TIME("amr_mesh::set_balance");
  if(m_verbosity > 5){
    pout() << "amr_mesh::set_balance" << endl;
  }
  
  m_which_balance = a_load;

  std::string balance;
  ParmParse pp("amr");
  pp.query("load_balance", balance);

  if(balance == "knapsack"){
    m_which_balance = load_balance::knapsack;
  }
  else if(balance == "elliptic"){
    m_which_balance = load_balance::elliptic;
  }
}

void amr_mesh::set_ghost_interpolation(const ghost_interpolation::which_type a_interp){
  CH_TIME("amr_mesh::set_ghost_interpolation");
  if(m_verbosity > 5){
    pout() << "amr_mesh::set_ghost_interpolation" << endl;
  }

  m_interp_type = a_interp;

  std::string interp_type;
  ParmParse pp("amr");
  pp.query("ghost_interp", interp_type);
  if(interp_type == "pwl"){
    m_interp_type = ghost_interpolation::pwl;
  }
  else if(interp_type == "quad"){
    m_interp_type = ghost_interpolation::quad;
  }
}

void amr_mesh::build_domains(){
  CH_TIME("amr_mesh::build_domains");
  if(m_verbosity > 5){
    pout() << "amr_mesh::build_domains" << endl;
  }

  const int nlevels = 1 + m_max_amr_depth;
  
  m_domains.resize(nlevels);
  m_dx.resize(nlevels);
  m_grids.resize(nlevels);
  m_mflg.resize(nlevels);

  m_eblg.resize(phase::num_phases);
  m_ebisl.resize(phase::num_phases);
  m_coarave.resize(phase::num_phases);
  m_quadcfi.resize(phase::num_phases);
  m_flux_reg.resize(phase::num_phases);
  m_old_quadcfi.resize(phase::num_phases);
  m_level_redist.resize(phase::num_phases);
  m_pwl_fillpatch.resize(phase::num_phases);
  m_pwl_interp.resize(phase::num_phases);
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
  m_pwl_fillpatch[phase::gas].resize(nlevels);
  m_pwl_interp[phase::gas].resize(nlevels);
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
  m_pwl_fillpatch[phase::solid].resize(nlevels);
  m_pwl_interp[phase::solid].resize(nlevels);
  m_coar_to_fine_redist[phase::solid].resize(nlevels);
  m_coar_to_coar_redist[phase::solid].resize(nlevels);
  m_fine_to_coar_redist[phase::solid].resize(nlevels);

  m_dx[0] = (m_physdom->get_prob_hi()[0] - m_physdom->get_prob_lo()[0])/m_num_cells[0];
  m_domains[0] = ProblemDomain(IntVect::Zero, m_num_cells - IntVect::Unit);

  for (int lvl = 1; lvl <= m_max_amr_depth; lvl++){
    m_dx[lvl]      = m_dx[lvl-1]/m_ref_ratios[lvl-1];
    m_domains[lvl] = m_domains[lvl-1];
    m_domains[lvl].refine(m_ref_ratios[lvl-1]);
  }
}

void amr_mesh::regrid(const Vector<IntVectSet>& a_tags, const int a_hardcap){
  CH_TIME("amr_mesh::regrid");
  if(m_verbosity > 1){
    pout() << "amr_mesh::regrid" << endl;
  }

  Vector<IntVectSet> tags = a_tags; // build_grids destroys tags, so copy them
  this->build_grids(tags, a_hardcap);
  this->define_eblevelgrid();  // Define EBLevelGrid objects on both phases
  this->define_mflevelgrid();  // Define MFLevelGrid
  this->define_eb_coar_ave();  // Define ebcoarseaverage on both phases
  this->define_eb_quad_cfi();  // Define nwoebquadcfinterp on both phases. 
  this->define_fillpatch();    // Define operator for piecewise linear interpolation of ghost cells
  this->define_ebpwl_interp(); // Define interpolator for piecewise interpolation of interior points
  this->define_flux_reg();     // Define flux register (phase::gas only)
  this->define_redist_oper();  // Define redistribution (phase::gas only)
  this->define_irreg_sten();   // Define irregular stencils
}

void amr_mesh::build_grids(Vector<IntVectSet>& a_tags, const int a_hardcap){
  CH_TIME("amr_mesh::build_grids");
  if(m_verbosity > 2){
    pout() << "amr_mesh::build_grids" << endl;
  }

  const int base      = 0;                        // I don't want this level to change
  const int top_level = a_tags.size() - 1;        // top_level is the finest level where we have tags
  Vector<Vector<Box> > new_boxes(1 + top_level);  // New boxes to be load balance
  Vector<Vector<Box> > old_boxes(1 + top_level);  // Old grids.

  const int hardcap = (a_hardcap == -1) ? m_max_amr_depth : a_hardcap;

  if(m_max_amr_depth > 0 && hardcap > 0){
    for (int lvl = 0; lvl <= top_level; lvl++){
      domainSplit(m_domains[lvl], old_boxes[lvl], m_max_box_size, m_blocking_factor);
    }
    
    // Berger-Rigoutsos grid generation
    BRMeshRefine mesh_refine(m_domains[0], m_ref_ratios, m_fill_ratio, m_blocking_factor, m_buffer_size, m_max_box_size);
    int new_finest_level = mesh_refine.regrid(new_boxes, a_tags, 0, top_level, old_boxes);
    m_finest_level = Min(new_finest_level, m_max_amr_depth); // Don't exceed m_max_amr_depth
    m_finest_level = Min(m_finest_level,   m_max_sim_depth); // Don't exceed maximum simulation depth
    m_finest_level = Min(m_finest_level,   hardcap);         // Don't exceed hardcap
  }
  else{
    new_boxes.resize(1);
    domainSplit(m_domains[0], new_boxes[0], m_max_box_size, m_blocking_factor);
    
    m_finest_level = 0;
  }

  // Do Morton ordering
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    mortonOrdering((Vector<Box>&)new_boxes[lvl]);
  }

  // Load balance boxes
  Vector<Vector<int> > proc_assign(1 + m_finest_level);
  this->loadbalance(proc_assign, new_boxes);

  // Define grids
  m_grids.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_grids[lvl] = DisjointBoxLayout();
    m_grids[lvl].define(new_boxes[lvl], proc_assign[lvl], m_domains[lvl]);
    m_grids[lvl].close();
  }
}

void amr_mesh::loadbalance(Vector<Vector<int> >& a_procs, Vector<Vector<Box> >& a_boxes){
  CH_TIME("amr_mesh::loadbalance");
  if(m_verbosity > 5){
    pout() << "amr_mesh::loadbalance" << endl;
  }

  // Level-by-level load balancing. This might change in the future. 
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(m_which_balance == load_balance::knapsack){
      load_balance::balance_knapsack(a_procs[lvl], a_boxes[lvl]);
    }
    else if(m_which_balance == load_balance::elliptic){
      load_balance::balance_elliptic(a_procs[lvl], a_boxes[lvl], m_mfis->get_ebis(phase::gas), m_domains[lvl], false);
    }
    else if(m_which_balance == load_balance::multifluid){
      load_balance::balance_multifluid(a_procs[lvl], a_boxes[lvl], m_mfis, m_domains[lvl], false);
    }
  }
}

void amr_mesh::compute_gradient(EBAMRCellData& a_gradient, const EBAMRCellData& a_phi){
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
	  EBArith::getFirstDerivStencilWidthOne(sten, vof, ebisbox, dir, m_dx[lvl], &cfivs[dit()], 0);
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

void amr_mesh::compute_gradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi){
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

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(!ebis_gas.isNull()){
      m_eblg[phase::gas][lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl],
									     m_domains[lvl],
									     m_ebghost,
									     ebis_gas));
      if(lvl > 0){
	m_eblg[phase::gas][lvl]->setMaxCoarseningRatio(m_ref_ratios[lvl-1], ebis_gas);
      }
      m_ebisl[phase::gas][lvl] = m_eblg[phase::gas][lvl]->getEBISL();
    }
    if(!ebis_sol.isNull()){
      m_eblg[phase::solid][lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl],
									       m_domains[lvl],
									       m_ebghost,
									       ebis_sol));
      if(lvl > 0){
	m_eblg[phase::solid][lvl]->setMaxCoarseningRatio(m_ref_ratios[lvl-1], ebis_sol);
      }
      m_ebisl[phase::solid][lvl] = m_eblg[phase::solid][lvl]->getEBISL();
    }
  }
}

void amr_mesh::define_mflevelgrid(){
  CH_TIME("amr_mesh::define_mflevelgrid");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_mflevelgrid" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBLevelGrid> eblgs;
    if(!ebis_gas.isNull()){
      eblgs.push_back(*m_eblg[phase::gas][lvl]);
    }
    if(!ebis_sol.isNull()){
      eblgs.push_back(*m_eblg[phase::solid][lvl]);
    }

    m_mflg[lvl] = RefCountedPtr<MFLevelGrid> (new MFLevelGrid(m_mfis, eblgs));
  }
}

void amr_mesh::define_eb_coar_ave(){
  CH_TIME("amr_mesh::define_eb_coar_ave");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eb_coar_ave" << endl;
  }

  const int comps = SpaceDim;

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
      	m_coarave[phase::gas][lvl] = RefCountedPtr<ebcoarseaverage> (new ebcoarseaverage(m_grids[lvl],
      											 m_grids[lvl-1],
      											 m_ebisl[phase::gas][lvl],
      											 m_ebisl[phase::gas][lvl-1],
      											 m_domains[lvl-1],
      											 m_ref_ratios[lvl-1],
      											 comps,
      											 &(*ebis_gas)));
      }

      if(!ebis_sol.isNull()){
      	m_coarave[phase::solid][lvl] = RefCountedPtr<ebcoarseaverage> (new ebcoarseaverage(m_grids[lvl],
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

void amr_mesh::define_fillpatch(){
  CH_TIME("amr_mesh::define_fillpatch");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_fillpatch" << endl;
  }

  const int comps     = SpaceDim;

  // Should these be input somehow?
  const int radius    = 1;
  const IntVect ghost = m_num_ghost*IntVect::Unit;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	const LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::gas][lvl]->getCFIVS());
	m_pwl_fillpatch[phase::gas][lvl] = RefCountedPtr<AggEBPWLFillPatch> (new AggEBPWLFillPatch(m_grids[lvl],
												   m_grids[lvl-1],
												   m_ebisl[phase::gas][lvl],
												   m_ebisl[phase::gas][lvl-1],
												   m_domains[lvl-1],
												   m_ref_ratios[lvl-1],
												   comps,
												   radius,
												   ghost,
												   !m_ebcf,
												   ebis_gas));
      }
      if(!ebis_sol.isNull()){
	const LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::solid][lvl]->getCFIVS());
	m_pwl_fillpatch[phase::solid][lvl] = RefCountedPtr<AggEBPWLFillPatch> (new AggEBPWLFillPatch(m_grids[lvl],
												     m_grids[lvl-1],
												     m_ebisl[phase::solid][lvl],
												     m_ebisl[phase::solid][lvl-1],
												     m_domains[lvl-1],
												     m_ref_ratios[lvl-1],
												     comps,
												     radius,
												     ghost,
												     !m_ebcf,
												     ebis_sol));
      }
    }
  }
}

void amr_mesh::define_ebpwl_interp(){
  CH_TIME("amr_mesh::define_ebpwl_interp");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_ebpwl_interp" << endl;
  }

  const int comps     = SpaceDim;

  // Should these be input somehow?
  const int radius    = 1;
  const IntVect ghost = m_num_ghost*IntVect::Unit;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	m_pwl_interp[phase::gas][lvl] = RefCountedPtr<EBPWLFineInterp> (new EBPWLFineInterp(m_grids[lvl],
											    m_grids[lvl-1],
											    m_ebisl[phase::gas][lvl],
											    m_ebisl[phase::gas][lvl-1],
											    m_domains[lvl-1],
											    m_ref_ratios[lvl-1],
											    comps,
											    ebis_gas));
      }
      if(!ebis_sol.isNull()){
	m_pwl_interp[phase::solid][lvl] = RefCountedPtr<EBPWLFineInterp> (new EBPWLFineInterp(m_grids[lvl],
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
	m_flux_reg[phase::gas][lvl] = RefCountedPtr<EBFluxRegister> (new EBFluxRegister(m_grids[lvl+1],
											m_grids[lvl],
											m_ebisl[phase::gas][lvl+1],
											m_ebisl[phase::gas][lvl],
											m_domains[lvl].domainBox(),
											m_ref_ratios[lvl],
											comps,
											ebis_gas));
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
    const bool has_fine = lvl < m_finest_level;

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
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_irreg_sten" << endl;
  }

  const int order = m_irreg_sten_order;
  const int rad   = m_irreg_sten_radius;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  if(!ebis_gas.isNull()){
    m_centroid_interp[phase::gas] = RefCountedPtr<irreg_amr_stencil<centroid_interp> >
      (new irreg_amr_stencil<centroid_interp>(m_grids,
					      m_ebisl[phase::gas],
					      m_domains,
					      m_dx,
					      m_finest_level,
					      order,
					      rad,
					      m_stencil_type));

    m_eb_centroid_interp[phase::gas] = RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> >
      (new irreg_amr_stencil<eb_centroid_interp>(m_grids,
						 m_ebisl[phase::gas],
						 m_domains,
						 m_dx,
						 m_finest_level,
						 order,
						 rad,
						 m_stencil_type));
  }

  if(!ebis_sol.isNull()){
    m_centroid_interp[phase::solid] = RefCountedPtr<irreg_amr_stencil<centroid_interp> >
      (new irreg_amr_stencil<centroid_interp>(m_grids,
					      m_ebisl[phase::solid],
					      m_domains,
					      m_dx,
					      m_finest_level,
					      order,
					      rad,
					      m_stencil_type));
    m_eb_centroid_interp[phase::solid] = RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> >
      (new irreg_amr_stencil<eb_centroid_interp>(m_grids,
						 m_ebisl[phase::solid],
						 m_domains,
						 m_dx,
						 m_finest_level,
						 order,
						 rad,
						 m_stencil_type));
  }
}

void amr_mesh::average_down(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(ebcell)" << endl;
  }
  const RefCountedPtr<EBIndexSpace>& ebis = m_mfis->get_ebis(a_phase);

  if(!ebis.isNull()){
    for (int lvl = m_finest_level; lvl > 0; lvl--){
      const int ncomps = a_data[lvl]->nComp();
      const Interval interv (0, ncomps-1);

      m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);
    }

    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      a_data[lvl]->exchange();
    }
  }
}

void amr_mesh::average_down(MFAMRFluxData& a_data){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(mfflux)" << endl;
  }

  EBAMRFluxData alias_g(1 + m_finest_level);
  EBAMRFluxData alias_s(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>());
      
    mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()){
      mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
    }
  }

  this->average_down(alias_g, phase::gas);
  if(!ebis_sol.isNull()){
    this->average_down(alias_s, phase::solid);
  }
}

void amr_mesh::average_down(MFAMRCellData& a_data){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(mfcell)" << endl;
  }
  
  EBAMRCellData alias_g(1 + m_finest_level);
  EBAMRCellData alias_s(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()){
      mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
    }
  }

  this->average_down(alias_g, phase::gas);
  if(!ebis_sol.isNull()){
    this->average_down(alias_s, phase::solid);
  }
}

void amr_mesh::average_down(EBAMRFluxData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(face)" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::average_down(EBAMRIVData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(iv)" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }

}

void amr_mesh::interp_ghost(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost(eb)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(eb)" << endl;
  }

  // for (int lvl = m_finest_level; lvl > 0; lvl--){
  //   const int ncomps = a_data[lvl]->nComp();
  //   const Interval interv(0, ncomps -1);

  //   CH_assert(a_data[lvl]->ghostVect() == m_num_ghost*IntVect::Unit);

  //   m_quadcfi[a_phase][lvl]->coarseFineInterp(*a_data[lvl], *a_data[lvl-1], 0, 0, ncomps);
  // }

  // for (int lvl = 0; lvl <= m_finest_level; lvl++){
  //   a_data[lvl]->exchange();
  // }

  if(m_interp_type == ghost_interpolation::pwl){
    this->interp_ghost_pwl(a_data, a_phase);
  }
  else if(m_interp_type == ghost_interpolation::quad){
    this->interp_ghost_quad(a_data, a_phase);
  }
  else{
    MayDay::Abort("amr_mesh::interp_ghost - unsupported interpolation type requested");
  }
}

void amr_mesh::interp_ghost(MFAMRCellData& a_data){
  CH_TIME("amr_mesh::interp_ghost(mf)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(mf)" << endl;
  }
  
  EBAMRCellData alias_g(1 + m_finest_level);
  EBAMRCellData alias_s(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);
  
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()){
      mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
    }
  }

  this->interp_ghost(alias_g, phase::gas);
  if(!ebis_sol.isNull()){
    this->interp_ghost(alias_s, phase::solid);
  }
}

void amr_mesh::interp_ghost_quad(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost(quad)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(quad)" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    CH_assert(a_data[lvl]->ghostVect() == m_num_ghost*IntVect::Unit);

    m_quadcfi[a_phase][lvl]->coarseFineInterp(*a_data[lvl], *a_data[lvl-1], 0, 0, ncomps);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::interp_ghost_pwl(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost(pwl)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(pwl)" << endl;
  }
  
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    CH_assert(a_data[lvl]->ghostVect() == m_num_ghost*IntVect::Unit);

    m_pwl_fillpatch[a_phase][lvl]->interpolate(*a_data[lvl], *a_data[lvl-1], *a_data[lvl-1], 0.0, 0.0, 0.0, interv);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }  
}

void amr_mesh::interpolate_to_centroids(EBAMRCellData& a_data, phase::which_phase a_phase){
  m_centroid_interp[a_phase]->apply(a_data);
}

void amr_mesh::set_verbosity(const int a_verbosity){
  CH_TIME("amr_mesh::set_verbosity");

  m_verbosity = a_verbosity;

  ParmParse pp("amr");
  pp.query("verbosity", m_verbosity);
}

void amr_mesh::set_coarsest_num_cells(const IntVect a_num_cells){
  m_num_cells = a_num_cells;

  ParmParse pp("amr");
  if(pp.contains("coarsest_domain")){
    Vector<int> cells;
    cells.resize(pp.countval("coarsest_domain"));
    CH_assert(cells.size() >= SpaceDim);
    pp.getarr("coarsest_domain", cells, 0, SpaceDim);

    m_num_cells = IntVect(D_DECL(cells[0], cells[1], cells[2]));
  }
}

void amr_mesh::set_max_amr_depth(const int a_max_amr_depth){
  m_max_amr_depth = a_max_amr_depth;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("max_amr_depth")){
      int depth;
      pp.get("max_amr_depth", depth);
      if(depth >= 0){
	m_max_amr_depth = depth;
      }
    }
  }
}

void amr_mesh::set_max_simulation_depth(const int a_max_sim_depth){
  m_max_sim_depth = a_max_sim_depth >= 0 ? a_max_sim_depth : 0;


  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("max_sim_depth")){
      int depth;
      pp.get("max_sim_depth", depth);
      if(depth >= 0){
	m_max_sim_depth = depth;
      }
    }
  }
}

void amr_mesh::set_ebcf(const bool a_ebcf){
  m_ebcf = a_ebcf;


  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("ebcf")){
      std::string str;
      pp.get("ebcf", str);
      if(str == "true"){
	m_ebcf = true;
      }
      else if(str == "false"){
	m_ebcf = false;
      }
    }
  }
}

void amr_mesh::set_refinement_ratio(const int a_refinement_ratio){
  CH_TIME("amr_mesh::set_refinement_ratio");
  m_ref_ratios.resize(1 + m_max_amr_depth);
  for (int lvl = 0; lvl <= m_max_amr_depth; lvl++){
    m_ref_ratios[lvl] = a_refinement_ratio;
  }

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("ref_rat")){
      Vector<int> ratios;
      ratios.resize(pp.countval("ref_rat"));
      pp.getarr("ref_rat", ratios, 0, ratios.size());
      
      m_ref_ratios = ratios;
    }
  }
}

void amr_mesh::set_refinement_ratios(const Vector<int> a_ref_ratios){
  m_ref_ratios = a_ref_ratios;
}

void amr_mesh::set_fill_ratio(const Real a_fill_ratio){
  m_fill_ratio = a_fill_ratio;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("fill_ratio")){
      Real fill_ratio;
      pp.get("fill_ratio", fill_ratio);
      if(fill_ratio > 0.0 && fill_ratio <= 1.0){
	m_fill_ratio = fill_ratio;
      }
    }
  }
}

void amr_mesh::set_finest_level(const int a_finest_level){
  m_finest_level = a_finest_level;
  m_finest_level = Min(m_finest_level, m_max_amr_depth); // Don't exceed m_max_amr_depth
  m_finest_level = Min(m_finest_level, m_max_sim_depth); // Don't exceed maximum simulation depth
}

void amr_mesh::set_grids(Vector<Vector<Box> >& a_boxes){

  // Do Morton ordering
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    mortonOrdering((Vector<Box>&)a_boxes[lvl]);
  }

  // Load balance boxes
  Vector<Vector<int> > proc_assign(1 + m_finest_level);
  this->loadbalance(proc_assign, a_boxes);

  // Define grids
  m_grids.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_grids[lvl] = DisjointBoxLayout();
    m_grids[lvl].define(a_boxes[lvl], proc_assign[lvl], m_domains[lvl]);
    m_grids[lvl].close();
  }

  this->define_eblevelgrid();  // Define EBLevelGrid objects on both phases
  this->define_mflevelgrid();  // Define MFLevelGrid
  this->define_eb_coar_ave();  // Define ebcoarseaverage on both phases
  this->define_eb_quad_cfi();  // Define nwoebquadcfinterp on both phases.
  this->define_fillpatch();    // Define operator for piecewise linear interpolation of ghost cells
  this->define_ebpwl_interp(); // Define interpolator for piecewise interpolation of interior points
  this->define_flux_reg();     // Define flux register (phase::gas only)
  this->define_redist_oper();  // Define redistribution (phase::gas only)
  this->define_irreg_sten();   // Define irregular stencils
}

void amr_mesh::set_max_box_size(const int a_max_box_size){
  CH_TIME("amr_mesh::set_max_box_size");
  m_max_box_size = a_max_box_size;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("max_box_size")){
      int box_size;
      pp.get("max_box_size", box_size);
      if(box_size >= 8 && box_size % 2 == 0){
	m_max_box_size = box_size;
      }
    }
  }
}

void amr_mesh::set_max_ebis_box_size(const int a_max_box_size){
  CH_TIME("amr_mesh::set_max_box_size");
  
  m_max_ebis_box_size = a_max_box_size;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("max_ebis_box")){
      int box_size;
      pp.get("max_ebis_box", box_size);
      if(box_size >= 8 && box_size % 2 == 0){
	m_max_ebis_box_size = box_size;
      }
    }
  }
}

void amr_mesh::set_buffer_size(const int a_buffer_size){
  m_buffer_size = a_buffer_size;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("buffer_size")){
      int buffer;
      pp.get("buffer_size", buffer);
      if(buffer > 0){
	m_buffer_size = buffer;
      }
    }
  }
}

void amr_mesh::set_irreg_growth(const int a_irreg_growth){
  m_irreg_growth = a_irreg_growth;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("irreg_growth")){
      int buffer;
      pp.get("irreg_growth", buffer);
      if(buffer > 0){
	m_irreg_growth = buffer;
      }
    }
  }
}

void amr_mesh::set_blocking_factor(const int a_blocking_factor){
  CH_TIME("amr_mesh::set_blocking_factor");
  m_blocking_factor = a_blocking_factor;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("blocking_factor")){
      int blocking;
      pp.get("blocking_factor", blocking);
      if(blocking >= 4 && blocking % 2 == 0){
	m_blocking_factor = blocking;
      }
    }
  }
}

void amr_mesh::set_eb_ghost(const int a_ebghost){
  CH_TIME("amr_mesh:.set_eb_ghost");
  m_ebghost = a_ebghost;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("eb_ghost")){
      int ebghost;
      pp.get("eb_ghost", ebghost);
      if(ebghost >= 2){
	m_ebghost = ebghost;
      }
    }
  }
}

void amr_mesh::set_num_ghost(const int a_num_ghost){
  CH_TIME("amr_mesh::set_num_ghost");
  m_num_ghost = a_num_ghost;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("num_ghost")){
      int ghost;
      pp.get("num_ghost", ghost);
      if(ghost >= 2){
	m_num_ghost = a_num_ghost;
      }
    }
  }
}

void amr_mesh::set_redist_rad(const int a_redist_rad){
  CH_TIME("amr_mesh::set_redist_rads");
  m_redist_rad = a_redist_rad;

  { // Get from input script
    ParmParse pp("amr");
    if(pp.contains("redist_radius")){
      int rad = 1;
      pp.get("redist_radius", rad);
      if(rad == 1 || rad == 2){
	m_redist_rad = a_redist_rad;
      }
    }
  }
}

void amr_mesh::set_irreg_sten_type(const stencil_type::which_type a_type){
  CH_TIME("amr_mesh::set_irreg_sten_type");
  m_stencil_type = a_type;

  
  std::string str = "taylor";
  ParmParse pp("amr");
  pp.query("stencil_type", str);
  if(str == "linear"){
    m_stencil_type = stencil_type::linear;
  }
  else if(str == "taylor"){
    m_stencil_type = stencil_type::taylor;
  }
  else if(str == "lsq"){
    m_stencil_type = stencil_type::lsq;
  }
}

void amr_mesh::set_irreg_sten_order(const int a_irreg_sten_order){
  CH_TIME("amr_mesh::irreg_sten_order");
  m_irreg_sten_order = a_irreg_sten_order;

  ParmParse pp("amr");
  int order = 1;
  pp.query("stencil_order", order);
  if(order == 1 || order == 2){
    m_irreg_sten_order = order;
  }
}

void amr_mesh::set_irreg_sten_radius(const int a_irreg_sten_radius){
  CH_TIME("amr_mesh::irreg_sten_radius");
  m_irreg_sten_radius = a_irreg_sten_radius;

  int radius = 1;
  ParmParse pp("amr");
  pp.query("stencil_radius", radius);
  if(radius == 1 || radius == 2){
    m_irreg_sten_radius = radius;
  }
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
  for (int lvl = 0; lvl < m_ref_ratios.size(); lvl++){
    CH_assert(m_ref_ratios[lvl] == 2 || m_ref_ratios[lvl] == 4);
    CH_assert(m_blocking_factor >= 4 && m_blocking_factor % m_ref_ratios[lvl] == 0);
  }
  CH_assert(m_max_box_size >= 8 && m_max_box_size % m_blocking_factor == 0);
  CH_assert(m_fill_ratio > 0. && m_fill_ratio <= 1.0);
  CH_assert(m_buffer_size > 0);
  CH_assert(m_irreg_sten_order == 1 ||  m_irreg_sten_order == 2);
  CH_assert(m_irreg_sten_radius == 1 || m_irreg_sten_radius == 2);

  // Make sure that the isotropic constraint isn't violated
  CH_assert(!m_physdom.isNull());
  const RealVect realbox = m_physdom->get_prob_hi() - m_physdom->get_prob_lo();
  for (int dir = 0; dir < SpaceDim - 1; dir++){
    CH_assert(realbox[dir]/m_num_cells[dir] == realbox[dir+1]/m_num_cells[dir+1]);
  }
}

bool amr_mesh::get_ebcf(){
  return m_ebcf;
}

int amr_mesh::get_finest_level(){
  return m_finest_level;
}

int amr_mesh::get_irreg_growth(){
  return m_irreg_growth;
}

int amr_mesh::get_max_amr_depth(){
  return m_max_amr_depth;
}

int amr_mesh::get_max_sim_depth(){
  return m_max_sim_depth;
}

int amr_mesh::get_num_ghost(){
  return m_num_ghost;
}

int amr_mesh::get_redist_rad(){
  return m_redist_rad;
}

int amr_mesh::get_blocking_factor(){
  return m_blocking_factor;
}

int amr_mesh::get_max_box_size(){
  return m_max_box_size;
}

int amr_mesh::get_max_ebis_box_size(){
  return m_max_ebis_box_size;
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

Vector<RefCountedPtr<MFLevelGrid> >& amr_mesh::get_mflg(){
  return m_mflg;
}

Vector<RefCountedPtr<ebcoarseaverage> >& amr_mesh::get_coarave(phase::which_phase a_phase){
  return m_coarave[a_phase];
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& amr_mesh::get_quadcfi(phase::which_phase a_phase){
  return m_quadcfi[a_phase];
}

Vector<RefCountedPtr<EBQuadCFInterp> >& amr_mesh::get_old_quadcfi(phase::which_phase a_phase){
  return m_old_quadcfi[a_phase];
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& amr_mesh::get_fillpatch(phase::which_phase a_phase){
  return m_pwl_fillpatch[a_phase];
}

Vector<RefCountedPtr<EBPWLFineInterp> >& amr_mesh::get_eb_pwl_interp(phase::which_phase a_phase){
  return m_pwl_interp[a_phase];
}

Vector<RefCountedPtr<EBFluxRegister> >&  amr_mesh::get_flux_reg(phase::which_phase a_phase){
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
