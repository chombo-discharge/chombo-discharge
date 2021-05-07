/*!
  @file   eddington_sp1.cpp
  @brief  Implementation of eddington_sp1.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Also plot kappa on output
*/

#include "eddington_sp1.H"
#include "data_ops.H"
#include "units.H"
#include "conductivitydomainbc_wrapper.H"

#include <ParmParse.H>
#include <EBAMRIO.H>
#include <BRMeshRefine.H>

#include <chrono>
#include <time.h>

#define eddington_sp1_feature 1 // Comment Feb. 14 2018: I think we can keep this - it appears to produce the correct physics.

namespace ChomboDischarge {

  eddington_sp1::eddington_sp1() : rte_solver() {
    m_name = "eddington_sp1";
    m_class_name = "eddington_sp1";

    m_verbosity  = -1;
    m_needs_setup = true;
    m_has_mg_stuff = false;
  }

  eddington_sp1::~eddington_sp1(){
  }

  int eddington_sp1::query_ghost() const{
    return 3;
  }

  void eddington_sp1::parse_options(){
    CH_TIME("eddington_sp1::parse_options");
    if(m_verbosity > 5){
      pout() << m_name + "::parse_options" << endl;
    }
  
    parse_domain_bc();    // Parses domain BC options
    parse_stationary();   // Parse stationary solver
    parse_plot_vars();    // Parses plot variables
    parse_gmg_settings(); // Parses solver parameters for geometric multigrid
  }

  void eddington_sp1::parse_runtime_options(){
    CH_TIME("eddington_sp1::parse_runtime_options");
    if(m_verbosity > 5){
      pout() << m_name + "::parse_runtime_options" << endl;
    }
  
    parse_stationary();   // Parse stationary solver
    parse_plot_vars();    // Parses plot variables
    parse_gmg_settings(); // Parses solver parameters for geometric multigrid
  }

  void eddington_sp1::parse_domain_bc(){
    CH_TIME("eddington_sp1::parse_domain_bc");
    if(m_verbosity > 5){
      pout() << m_name + "::parse_domain_bc" << endl;
    }

    allocate_wall_bc();

    // Get BC from input script
    ParmParse pp(m_class_name.c_str());
  
    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){
	const Side::LoHiSide side = sit();
	
	std::string str_dir;
	if(dir == 0){
	  str_dir = "x";
	}
	else if(dir == 1){
	  str_dir = "y";
	}
	else if(dir == 2){
	  str_dir = "z";
	}

	if(side == Side::Lo){
	  std::string type;
	  std::string bc_string = "bc_" + str_dir + "_low";
	  if(pp.contains(bc_string.c_str())){
	    pp.get(bc_string.c_str(), type);
	    if(type == "neumann"){
	      this->set_neumann_wall_bc(dir, Side::Lo, 0.0);
	    }
	    else if(type == "robin"){
	      this->set_robin_wall_bc(dir, Side::Lo, 0.0);
	    }
	    else {
	      std::string error = "eddington_sp1::eddington_sp1 - unknown bc requested for " + bc_string;
	      MayDay::Abort(error.c_str());
	    }
	  }
	}
	else if(side == Side::Hi){
	  std::string type;
	  std::string bc_string = "bc_" + str_dir + "_high";
	  if(pp.contains(bc_string.c_str())){
	    pp.get(bc_string.c_str(), type);
	    if(type == "neumann"){
	      this->set_neumann_wall_bc(dir, Side::Hi, 0.0);
	    }
	    else if(type == "robin"){
	      this->set_robin_wall_bc(dir, Side::Hi, 0.0);
	    }
	    else {
	      std::string error = "eddington_sp1::eddington_sp1 - unknown bc requested for " + bc_string;
	      MayDay::Abort(error.c_str());
	    }
	  }
	}
      }
    }
  }

  void eddington_sp1::parse_stationary(){
    CH_TIME("eddington_sp1::parse_stationary");
    if(m_verbosity > 5){
      pout() << m_name + "::parse_stationary" << endl;
    }

    ParmParse pp(m_class_name.c_str());
    std::string str;
  
    pp.get("stationary", str);
    m_stationary = (str == "true") ? true : false;
  
    pp.get("use_tga", str);
    m_use_tga = (str == "true") ? true : false;
  }

  void eddington_sp1::parse_plot_vars(){
    CH_TIME("eddington_sp1::parse_plot_vars");
    if(m_verbosity > 5){
      pout() << m_name + "::parse_plot_vars" << endl;
    }

    m_plot_phi = false;
    m_plot_src = false;

    ParmParse pp(m_class_name.c_str());
    const int num = pp.countval("plt_vars");
    Vector<std::string> str(num);
    pp.getarr("plt_vars", str, 0, num);

    for (int i = 0; i < num; i++){
      if(     str[i] == "phi") m_plot_phi = true;
      else if(str[i] == "src") m_plot_src = true;
    }
  }

  void eddington_sp1::parse_gmg_settings(){
    ParmParse pp(m_class_name.c_str());

    std::string str;

    pp.get("gmg_verbosity",   m_gmg_verbosity);
    pp.get("gmg_coarsen",     m_gmg_coarsen);
    pp.get("gmg_pre_smooth",  m_gmg_pre_smooth);
    pp.get("gmg_post_smooth", m_gmg_post_smooth);
    pp.get("gmg_bott_smooth", m_gmg_bot_smooth);
    pp.get("gmg_max_iter",    m_gmg_max_iter);
    pp.get("gmg_min_iter",    m_gmg_min_iter);
    pp.get("gmg_tolerance",   m_gmg_eps);
    pp.get("gmg_hang",        m_gmg_hang);
    pp.get("gmg_bottom_drop", m_bottom_drop);

    // Bottom solver
    pp.get("gmg_bottom_solver", str);
    if(str == "simple"){
      m_bottomsolver = 0;
    }
    else if(str == "bicgstab"){
      m_bottomsolver = 1;
    }
    else{
      MayDay::Abort("eddington_sp1::parse_gmg_settings - unknown bottom solver requested");
    }

    // Relaxation type
    pp.get("gmg_relax_type", str);
    if(str == "gsrb"){
      m_gmg_relax_type = relax::gsrb_fast;
    }
    else if(str == "jacobi"){
      m_gmg_relax_type = relax::jacobi;
    }
    else if(str == "gauss_seidel"){
      m_gmg_relax_type = relax::gauss_seidel;
    }
    else{
      MayDay::Abort("eddington_sp1::parse_gmg_settings - unknown relaxation method requested");
    }

    // Cycle type
    pp.get("gmg_cycle", str);
    if(str == "vcycle"){
      m_gmg_type = amrmg::vcycle;
    }
    else{
      MayDay::Abort("eddington_sp1::parse_gmg_settings - unknown cycle type requested");
    }

    // No lower than 2. 
    if(m_bottom_drop < 2){
      m_bottom_drop = 2;
    }
  }

  void eddington_sp1::parse_reflection(){
    CH_TIME("eddington_sp1::parse_reflectivity");
    if(m_verbosity > 5){
      pout() << m_name + "::parse_reflectivity" << endl;
    }
  
    ParmParse pp(m_class_name.c_str());

    Real r;
    pp.get("reflectivity", r);

    m_r1 = r/(2.0);
    m_r2 = r/(3.0);
  }

  void eddington_sp1::allocate_wall_bc(){
    CH_TIME("eddington_sp1::allocate_wall_bc");
    if(m_verbosity > 5){
      pout() << "eddington_sp1::allocate_wall_bc" << endl;
    }
    m_wallbc.resize(2*SpaceDim);
    for (int i = 0; i < 2*SpaceDim; i++){
      m_wallbc[i] = RefCountedPtr<wall_bc> (NULL);
    }
  }

  void eddington_sp1::pre_regrid(const int a_base, const int a_old_finest_level){
    CH_TIME("eddington_sp1::pre_regrid");
    if(m_verbosity > 5){
      pout() << m_name + "::pre_regrid" << endl;
    }

    const int ncomp = 1;
    const int finest_level = m_amr->get_finest_level();

    m_amr->allocate(m_cache, m_realm, m_phase, ncomp);

    for (int lvl = 0; lvl <= finest_level; lvl++){
      m_state[lvl]->localCopyTo(*m_cache[lvl]);
    }
  }

  void eddington_sp1::set_reflection_coefficients(const Real a_r1, const Real a_r2){
    CH_TIME("eddington_sp1::set_reflection_coefficients");
    if(m_verbosity > 5){
      pout() << m_name + "::set_reflection_coefficients" << endl;
    }
  
    m_r1 = a_r1;
    m_r2 = a_r2;
  }

  void eddington_sp1::allocate_internals(){
    CH_TIME("eddington_sp1::allocate_internals");
    if(m_verbosity > 5){
      pout() << m_name + "::allocate_internals" << endl;
    }
  
    const int ncomp = 1;

    m_amr->allocate(m_aco,       m_realm, m_phase, ncomp);
    m_amr->allocate(m_bco,       m_realm, m_phase, ncomp);
    m_amr->allocate(m_bco_irreg, m_realm, m_phase, ncomp);
    m_amr->allocate(m_state,     m_realm, m_phase, ncomp);
    m_amr->allocate(m_source,    m_realm, m_phase, ncomp);
    m_amr->allocate(m_resid,     m_realm, m_phase, ncomp);

    data_ops::set_value(m_resid,  0.0);
    data_ops::set_value(m_state,  0.0);
    data_ops::set_value(m_source, 0.0);

    this->set_aco_and_bco();
  }

  void eddington_sp1::deallocate_internals(){
    m_amr->deallocate(m_aco);
    m_amr->deallocate(m_bco);
    m_amr->deallocate(m_bco_irreg);
    m_amr->deallocate(m_state);
    m_amr->deallocate(m_source);
    m_amr->deallocate(m_resid);
  }

  void eddington_sp1::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {
    CH_TIME("eddington_sp1::regrid");
    if(m_verbosity > 5){
      pout() << m_name + "::regrid" << endl;
    }

    const int comp  = 0;
    const int ncomp = 1;
    const Interval interv(comp, comp);

    this->allocate_internals();

    Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->get_eb_pwl_interp(m_realm, m_phase);

    // These levels have not changed
    for (int lvl = 0; lvl <= Max(0, a_lmin-1); lvl++){
      m_cache[lvl]->copyTo(*m_state[lvl]); // Base level should never change, but ownership might
    }

    // These levels have changed
    for (int lvl = Max(1,a_lmin); lvl <= a_new_finest_level; lvl++){
      interpolator[lvl]->interpolate(*m_state[lvl], *m_state[lvl-1], interv);

      if(lvl <= a_old_finest_level){
	m_cache[lvl]->copyTo(*m_state[lvl]);
      }
    }

    m_needs_setup = true;
  }

  void eddington_sp1::register_operators(){
    CH_TIME("eddington_sp1::register_operators");
    if(m_verbosity > 5){
      pout() << m_name + "::register_operators" << endl;
    }

    if(m_amr.isNull()){
      MayDay::Abort("eddington_sp1::register_operators - need to set amr_mesh!");
    }
    else{
      m_amr->register_operator(s_eb_coar_ave,     m_realm, m_phase);
      m_amr->register_operator(s_eb_fill_patch,   m_realm, m_phase);
      m_amr->register_operator(s_eb_flux_reg,     m_realm, m_phase);
      m_amr->register_operator(s_eb_quad_cfi,     m_realm, m_phase);
      m_amr->register_operator(s_eb_gradient,     m_realm, m_phase);
      m_amr->register_operator(s_eb_irreg_interp, m_realm, m_phase);
    }
  }

  bool eddington_sp1::advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi){
    CH_TIME("eddington_sp1::advance(ebamrcell, ebamrcell)");
    if(m_verbosity > 5){
      pout() << m_name + "::advance(ebamrcell, ebamrcell)" << endl;
    }

    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    if(m_needs_setup){
      this->setup_gmg();
    }

    bool converged;

    // Must have a dummy for checking initial residual
    EBAMRCellData dummy;
    EBAMRCellData source;
    m_amr->allocate(dummy,  m_realm, m_phase, ncomp);
    m_amr->allocate(source, m_realm, m_phase, ncomp);
    data_ops::set_value(dummy, 0.0);

    // Various source term manipulations. 
    data_ops::set_value(source, 0.0);
    data_ops::incr(source, a_source, 1.0);
#if eddington_sp1_feature
    data_ops::scale(source, 1./units::s_c0); // Source should be scaled by 1./c0
#endif
    if(m_stationary){ // Should kappa-scale for transient solvres
      data_ops::kappa_scale(source);
    }

    Vector<LevelData<EBCellFAB>* > phi, rhs, res, zero;
    m_amr->alias(phi,  a_state);
    m_amr->alias(rhs,  source);
    m_amr->alias(res,  m_resid);
    m_amr->alias(zero, dummy);

    if(m_stationary){
      const Real phi_resid  = m_gmg_solver->computeAMRResidual(phi,  rhs, finest_level, 0); // Incoming residual
      const Real zero_resid = m_gmg_solver->computeAMRResidual(zero, rhs, finest_level, 0); // Zero residual

      if(phi_resid > zero_resid*m_gmg_eps){ // Residual is too large
	m_gmg_solver->m_convergenceMetric = zero_resid;
	m_gmg_solver->solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);
      
	const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
	if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
	  converged = true;
	}
      }
      else{ // Solution is already good enough
	converged = true;
      }
      m_gmg_solver->revert(phi, rhs, finest_level, 0);

      data_ops::set_covered_value(a_state, 0, 0.0);
    }
    else{
      if(m_use_tga){
#if eddington_sp1_feature
	m_tgasolver->oneStep(res, phi, rhs, units::s_c0*a_dt, 0, finest_level, m_time);
#else
	m_tgasolver->oneStep(res, phi, rhs, a_dt, 0, finest_level, m_time);
#endif
      }
      else{
#if eddington_sp1_feature
	m_eulersolver->oneStep(res, phi, rhs, units::s_c0*a_dt, 0, finest_level, false);
#else
	m_eulersolver->oneStep(res, phi, rhs, a_dt, 0, finest_level, false);
#endif
      }
    
      const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
      if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
	converged = true;
      }

      // We solve onto res, copy back to state
      data_ops::copy(a_state, m_resid);
    }

    m_amr->average_down(a_state, m_realm, m_phase);
    m_amr->interp_ghost(a_state, m_realm, m_phase);

    data_ops::floor(a_state, 0.0);

    return converged;
  }

  void eddington_sp1::setup_gmg(){
    CH_TIME("eddington_sp1::setup_gmg");
    if(m_verbosity > 5){
      pout() << m_name + "::setup_gmg" << endl;
    }

    if(!m_has_mg_stuff){
      this->define_mg_levels();
      m_has_mg_stuff = true;
    }
  
    this->set_coefficients();       // Set coefficients, kappa, aco, bco
    this->setup_operator_factory(); // Set the operator factory
    this->setup_multigrid();        // Set up the AMR multigrid solver

    if(!m_stationary){
      if(m_use_tga){
	this->setup_tga();
      }
      else{
	this->setup_euler();
      }
    }

    m_needs_setup = false;
  }

  void eddington_sp1::set_coefficients(){
    CH_TIME("eddington_sp1::set_coefficients");
    if(m_verbosity > 5){
      pout() << m_name + "::set_coefficients" << endl;
    }

    const int ncomp = 1;
    const int ghost = 3;

    m_amr->allocate(m_aco,        m_realm, m_phase, ncomp, ghost);
    m_amr->allocate(m_bco,        m_realm, m_phase, ncomp, ghost);
    m_amr->allocate(m_bco_irreg,  m_realm, m_phase, ncomp, ghost);

    this->set_aco_and_bco();
  }

  void eddington_sp1::set_aco_and_bco(){
    CH_TIME("eddington_sp1::set_aco_and_bco");
    if(m_verbosity > 5){
      pout() << m_name + "::set_aco_and_bco" << endl;
    }

    // This loop fills aco with kappa and bco_irreg with 1./kappa
    if(m_rte_species->constant_kappa()){
      const Real kap = m_rte_species->get_kappa(RealVect::Zero);
      data_ops::set_value(m_aco, kap);
      data_ops::set_value(m_bco, 1./kap);
      data_ops::set_value(m_bco_irreg, 1./kap);
    }
    else{ // If kappa is not constant, we need to go through each cell to determine it
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const RealVect origin = m_amr->get_prob_lo();
	const Real dx         = m_amr->get_dx()[lvl];
    
	LevelData<EBCellFAB>& aco            = *m_aco[lvl];
	LevelData<EBFluxFAB>& bco            = *m_bco[lvl];
	LevelData<BaseIVFAB<Real> >& bco_irr = *m_bco_irreg[lvl];

	for (DataIterator dit = aco.dataIterator(); dit.ok(); ++dit){
	  const Box box = (m_amr->get_grids(m_realm)[lvl]).get(dit());
	  this->set_aco_and_bco_box(aco[dit()], bco_irr[dit()], box, origin, dx, lvl, dit());
	}
      }

      m_amr->average_down(m_aco, m_realm, m_phase);
      m_amr->interp_ghost(m_aco, m_realm, m_phase);
      data_ops::average_cell_to_face_allcomps(m_bco, m_aco, m_amr->get_domains()); // Average aco onto face
      data_ops::invert(m_bco); // Make m_bco = 1./kappa
    }

#if eddington_sp1_feature // Different scaling for the RTE
    data_ops::scale(m_aco,       1.0);       // aco = c*kappa
    data_ops::scale(m_bco,       1.0/(3.0)); // bco = c/(3*kappa)
    data_ops::scale(m_bco_irreg, 1.0/(3.0)); // bco = c/(3*kappa)
#else // Original code before different scaling
    data_ops::scale(m_aco,       units::s_c0);       // aco = c*kappa
    data_ops::scale(m_bco,       units::s_c0/(3.0)); // bco = c/(3*kappa)
    data_ops::scale(m_bco_irreg, units::s_c0/(3.0)); // bco = c/(3*kappa)
#endif
  }

  void eddington_sp1::set_aco_and_bco_box(EBCellFAB&       a_aco,
					  BaseIVFAB<Real>& a_bco,
					  const Box        a_box,
					  const RealVect   a_origin,
					  const Real       a_dx,
					  const int        a_lvl,
					  const DataIndex& a_dit){
    CH_TIME("eddington_sp1::set_aco_and_bco_box");
    if(m_verbosity > 10){
      pout() << m_name + "::set_aco_and_bco_box" << endl;
    }
  
    const int comp  = 0;
    const int ncomp = 1;
  
    const EBISBox& ebisbox = a_aco.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    // Regular aco
    BaseFab<Real>& aco_fab = a_aco.getSingleValuedFAB();
    for (BoxIterator bit(a_box); bit.ok(); ++bit){
      const IntVect iv = bit();

      const RealVect pos = a_origin + iv*a_dx*RealVect::Unit;
      aco_fab(iv, comp) = m_rte_species->get_kappa(pos);
    }


    // Irregular stuff
    VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_lvl])[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, a_origin);
      const Real tmp = m_rte_species->get_kappa(pos);
      a_aco(vof, comp) = tmp;
      a_bco(vof, comp) = 1./tmp;
    }
  }

  void eddington_sp1::define_mg_levels(){
    CH_TIME("eddington_sp1::define_mg_levels");
    if(m_verbosity > 5){
      pout() << m_name + "::define_mg_levels" << endl;
    }

    const int coar_ref = 2;

    Vector<ProblemDomain> m_mg_domains(0);
    Vector<DisjointBoxLayout> m_mg_grids(0);
  
    m_mg_domains.resize(0);
    m_mg_grids.resize(0);
    m_mg_levelgrids.resize(0);

    // Get some stuff from amr_mesh on how to decompose the levels
    const Vector<ProblemDomain>& domains = m_amr->get_domains();
    const int max_box_size               = m_amr->get_max_box_size();
    const int blocking_factor            = m_amr->get_blocking_factor();
    const int num_ebghost                = m_amr->get_eb_ghost();

    int num_coar       = 0;
    bool has_coar      = true;
    ProblemDomain fine = domains[0];

    // Coarsen problem domains and create grids
    while(num_coar < m_gmg_coarsen || !has_coar){

      // Check if we can coarsen
      const ProblemDomain coar = fine.coarsen(coar_ref);
      const Box coar_box       = coar.domainBox();
      for (int dir = 0; dir < SpaceDim; dir++){
	if(coar_box.size()[dir] < max_box_size || coar_box.size()[dir]%max_box_size != 0){
	  has_coar = false;
	}
      }

      if(has_coar){
	// Split the domain into pieces, then order and load balance them
	Vector<Box> boxes;
	Vector<int> proc_assign;
	domainSplit(coar, boxes, max_box_size, blocking_factor);
	mortonOrdering(boxes);
	load_balance::make_balance(proc_assign, boxes);

	// Add problem domain and grid
	m_mg_domains.push_back(coar);
	m_mg_grids.push_back(DisjointBoxLayout(boxes, proc_assign, coar));

	// Define the EBLevelGrids
	const int idx = m_mg_grids.size() - 1; // Last element added
	m_mg_levelgrids.push_back(EBLevelGrid(m_mg_grids[idx],
					      m_mg_domains[idx],
					      num_ebghost,
					      m_ebis));

	// Next iterate
	fine = coar;
	num_coar++;
      }
      else{
	break;
      }
    }
  }

  void eddington_sp1::setup_operator_factory(){
    CH_TIME("eddington_sp1::setup_operator_factory");
    if(m_verbosity > 5){
      pout() << m_name + "::setup_operator_factory" << endl;
    }

    const int finest_level                 = m_amr->get_finest_level();
    const int ghost                        = m_amr->get_num_ghost();
    const Vector<DisjointBoxLayout>& grids = m_amr->get_grids(m_realm);
    const Vector<int>& refinement_ratios   = m_amr->get_ref_rat();
    const Vector<ProblemDomain>& domains   = m_amr->get_domains();
    const Vector<Real>& dx                 = m_amr->get_dx();
    const RealVect& origin                 = m_amr->get_prob_lo();
    const Vector<EBISLayout>& ebisl        = m_amr->get_ebisl(m_realm, m_phase);
  
    const Vector<RefCountedPtr<EBQuadCFInterp> >& quadcfi  = m_amr->get_old_quadcfi(m_realm, m_phase);
    const Vector<RefCountedPtr<EBFluxRegister> >& fastFR   = m_amr->get_flux_reg(m_realm, m_phase);

    Vector<EBLevelGrid> levelgrids;

    for (int lvl = 0; lvl <= finest_level; lvl++){ 
      levelgrids.push_back(*(m_amr->get_eblg(m_realm, m_phase)[lvl])); // amr_mesh uses RefCounted levelgrids. EBConductivityOp does not. 
    }

#if 0
    Vector<EBLevelGrid> mg_levelgrids;
    Vector<RefCountedPtr<EBLevelGrid> >& mg_eblg = m_amr->get_mg_eblg(m_realm, m_phase);
    for (int lvl = 0; lvl < mg_eblg.size(); lvl++){
      mg_levelgrids.push_back(*mg_eblg[lvl]);
    }
#endif

    // Appropriate coefficients. 
    const Real alpha =  1.0;
    const Real beta  = -1.0;

    // Appropriate coefficients for this type of Robin BC
    m_robinco = RefCountedPtr<larsen_coefs> (new larsen_coefs(m_rte_species, m_r1, m_r2));

    // Domain BC
#if 0
    m_domfact = RefCountedPtr<robinconductivitydomainbcfactory> (new robinconductivitydomainbcfactory());
    m_domfact->set_coefs(m_robinco);
#else
    RefCountedPtr<BaseDomainBCFactory> domfact = RefCountedPtr<BaseDomainBCFactory>(NULL);
    conductivitydomainbc_wrapper_factory* bcfact = new conductivitydomainbc_wrapper_factory();
    Vector<RefCountedPtr<robin_coef> > coefs(2*SpaceDim, m_robinco);
    bcfact->set_wallbc(m_wallbc);
    bcfact->set_robin_coefs(coefs);
    domfact = RefCountedPtr<BaseDomainBCFactory> (bcfact);
#endif

    // EBBC
    m_ebfact  = RefCountedPtr<robinconductivityebbcfactory> (new robinconductivityebbcfactory(origin));
    m_ebfact->set_coefs(m_robinco);
    m_ebfact->set_type(stencil_type::lsq);

    //  Make relaxation type into int code
    int relax_type = 0;
    if(m_gmg_relax_type == relax::jacobi){
      relax_type = 0;
    }
    else if(m_gmg_relax_type == relax::gauss_seidel){
      relax_type = 1;
    }
    else if(m_gmg_relax_type == relax::gsrb_fast){
      relax_type = 2;
    }

    // Create operator factory.
    m_opfact = RefCountedPtr<ebconductivityopfactory> (new ebconductivityopfactory(levelgrids,
										   quadcfi,
										   fastFR,
										   alpha,
										   beta,
										   m_aco.get_data(),
										   m_bco.get_data(),
										   m_bco_irreg.get_data(),
										   dx[0],
										   refinement_ratios,
										   domfact,
										   m_ebfact,
										   ghost*IntVect::Unit,
										   ghost*IntVect::Unit,
										   relax_type,
										   m_bottom_drop,
										   -1,
										   m_mg_levelgrids));
  }

  void eddington_sp1::setup_multigrid(){
    CH_TIME("eddington_sp1::setup_multigrid");
    if(m_verbosity > 5){
      pout() << m_name + "::setup_multigrid" << endl;
    }

    const int finest_level       = m_amr->get_finest_level();
    const ProblemDomain coar_dom = m_amr->get_domains()[0];

    // Select bottom solver
    LinearSolver<LevelData<EBCellFAB> >* botsolver = NULL;
    if(m_bottomsolver == 0){
      m_simple_solver.setNumSmooths(m_numsmooth);
      botsolver = &m_simple_solver;
    }
    else{
      botsolver = &m_bicgstab;
    }

    // Make m_gmg_type into an int for multigrid
    int gmg_type;
    if(m_gmg_type == amrmg::full){
      gmg_type = 0;
    }
    else if(m_gmg_type == amrmg::vcycle){
      gmg_type = 1;
    }
    else if(m_gmg_type == amrmg::fcycle){
      gmg_type = 2;
    }

    m_gmg_solver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
    m_gmg_solver->define(coar_dom, *m_opfact, botsolver, 1 + finest_level);
    m_gmg_solver->setSolverParameters(m_gmg_pre_smooth,
				      m_gmg_post_smooth,
				      m_gmg_bot_smooth,
				      gmg_type,
				      m_gmg_max_iter,
				      m_gmg_eps,
				      m_gmg_hang,
				      1.E-99); // Residue set through other means
    m_gmg_solver->m_imin    = m_gmg_min_iter;
    m_gmg_solver->m_verbosity = m_gmg_verbosity;

    // Dummies for init
    const int ncomp = 1;
    EBAMRCellData dummy1, dummy2;
    m_amr->allocate(dummy1, m_realm, m_phase, ncomp);
    m_amr->allocate(dummy2, m_realm, m_phase, ncomp);
    data_ops::set_value(dummy1, 0.0);
    data_ops::set_value(dummy2, 0.0);

    // Aliasing
    Vector<LevelData<EBCellFAB>* > phi, rhs;
    m_amr->alias(phi, dummy1);
    m_amr->alias(rhs, dummy2);

    // Init solver
    m_gmg_solver->init(phi, rhs, finest_level, 0);
  }

  void eddington_sp1::setup_tga(){
    CH_TIME("eddington_sp1::setup_tga");
    if(m_verbosity > 5){
      pout() << m_name + "::setup_tga" << endl;
    }
  
    const int finest_level       = m_amr->get_finest_level();
    const ProblemDomain coar_dom = m_amr->get_domains()[0];
    const Vector<int> ref_rat    = m_amr->get_ref_rat();

    m_tgasolver = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
      (new AMRTGA<LevelData<EBCellFAB> > (m_gmg_solver, *m_opfact, coar_dom, ref_rat, 1 + finest_level, m_gmg_solver->m_verbosity));

    // Must init gmg for TGA
    Vector<LevelData<EBCellFAB>* > phi, rhs;
    m_amr->alias(phi, m_state);
    m_amr->alias(rhs, m_source);
    m_gmg_solver->init(phi, rhs, finest_level, 0);
  }

  void eddington_sp1::setup_euler(){
    CH_TIME("eddington_sp1::setup_euler");
    if(m_verbosity > 5){
      pout() << m_name + "::setup_euler" << endl;
    }

    const int finest_level       = m_amr->get_finest_level();
    const ProblemDomain coar_dom = m_amr->get_domains()[0];
    const Vector<int> ref_rat    = m_amr->get_ref_rat();

    m_eulersolver = RefCountedPtr<EBBackwardEuler> 
      (new EBBackwardEuler (m_gmg_solver, *m_opfact, coar_dom, ref_rat, 1 + finest_level, m_gmg_solver->m_verbosity));

    // Note: If this crashes, try to init gmg first
  }

  void eddington_sp1::compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state){
    CH_TIME("eddington_sp1::compute_boundary_flux");
    if(m_verbosity > 5){
      pout() << m_name + "::compute_boundary_flux" << endl;
    }

    const int finest_level = m_amr->get_finest_level();
  
    irreg_amr_stencil<eb_centroid_interp>& sten = m_amr->get_eb_centroid_interp_stencils(m_realm, m_phase);
    for(int lvl = 0; lvl <= finest_level; lvl++){
      sten.apply(*a_ebflux[lvl], *a_state[lvl], lvl, true);
    }

    m_amr->average_down(a_ebflux, m_realm, m_phase);

    data_ops::scale(a_ebflux, 0.5*units::s_c0);
  }

  void eddington_sp1::compute_domain_flux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_data){
    CH_TIME("eddington_sp1::compute_domain_flux");
    if(m_verbosity > 5){
      pout() << m_name + "::compute_domain_flux" << endl;
    }


    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const int ncomp = a_data[lvl]->nComp();
      
      const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const EBCellFAB& data         = (*a_data[lvl])[dit()];
	const EBISBox& ebisbox        = ebisl[dit()];
	const BaseFab<Real>& data_fab = data.getSingleValuedFAB();
      
	for (int dir = 0; dir < SpaceDim; dir++){
	  for (SideIterator sit; sit.ok(); ++sit){
	    BaseIFFAB<Real>& extrap = (*a_domainflux[lvl])[dit()](dir, sit());

	    const IntVectSet& ivs  = extrap.getIVS();
	    const EBGraph& ebgraph = extrap.getEBGraph();

	    // Extrapolate to the boundary. Use face-centered stuff for all faces (also multivalued ones)
	    const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
	    for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	      const FaceIndex& face = faceit();

	      const int sgn = sign(sit()); // Lo = -1, Hi = 1
	    
	      const VolIndex& vof = face.getVoF(flip(sit()));
	      const IntVect iv0   = vof.gridIndex();
	      const IntVect iv1   = iv0 - sgn*BASISV(dir);

	      if(ebisbox.isCovered(iv0)){ // Just provide some bogus data because the face
		for (int comp = 0; comp < ncomp; comp++){
		  extrap(face, comp) = 0.0;
		}
	      }
	      else{
		if(!ebisbox.isCovered(iv1)){ // linear extrapolation
		  for (int comp = 0; comp < ncomp; comp++){
		    extrap(face, comp) = 1.5*data_fab(iv0, comp) - 0.5*data_fab(iv1, comp); // Should be ok
		  }
		}
		else{ // Not enough cells available, use cell-centered only
		  for (int comp = 0; comp < ncomp; comp++){
		    extrap(face, comp) = data_fab(iv0, comp);
		  }
		}
	      }

	      // Necessary scaling
	      for (int comp = 0; comp < ncomp; comp++){
		extrap(face, comp) = 0.5*units::s_c0*extrap(face, comp);
	      }
	    }
	  }
	}
      }
    }
  }

  void eddington_sp1::compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state){
    CH_TIME("eddington_sp1::compute_flux");
    if(m_verbosity > 5){
      pout() << m_name + "::compute_flux" << endl;
    }

    const int finest_level = m_amr->get_finest_level();

    m_amr->compute_gradient(a_flux, a_state, m_realm, m_phase); // flux = grad(phi)
    for (int lvl = 0; lvl <= finest_level; lvl++){
      data_ops::divide_scalar(*a_flux[lvl], *m_aco[lvl]);   // flux = grad(phi)/(c*kappa)
      data_ops::scale(*a_flux[lvl], -units::s_c0*units::s_c0/3.0);  // flux = -c*grad(phi)/3.
    }

    m_amr->average_down(a_flux, m_realm, m_phase);
    m_amr->interp_ghost(a_flux, m_realm, m_phase);
  }


  void eddington_sp1::compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state){
    CH_TIME("eddington_sp1::compute_density");
    if(m_verbosity > 5){
      pout() << m_name + "::compute_density" << endl;
    }

    const int finest_level = m_amr->get_finest_level();
  
    for (int lvl = 0; lvl <= finest_level; lvl++){
      const Interval interv(0,0);
      a_state[lvl]->localCopyTo(interv, *a_isotropic[lvl], interv);
    }
  }

  void eddington_sp1::write_plot_file(){
    CH_TIME("eddington_sp1::write_plot_file");
    if(m_verbosity > 5){
      pout() << m_name + "::write_plot_file" << endl;
    }

    char file_char[1000];
    sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_step, SpaceDim);

    const int ncomps = 3 + SpaceDim;
    Vector<string> names(ncomps);
    names[0] = "density";
    names[1] = "x-flux";
    names[2] = "y-flux";
    if(SpaceDim == 3){
      names[3] = "z-flux";
      names[4] = "isotropic source";
      names[5] = "residue";
    }
    else{
      names[3] = "isotropic source";
      names[4] = "residue";
    }

    // Compute the flux
    EBAMRCellData flux;
    m_amr->allocate(flux, m_realm, m_phase, SpaceDim);
    this->compute_flux(flux, m_state);

    // Allocate output storage
    EBAMRCellData output;
    m_amr->allocate(output, m_realm, m_phase, ncomps, 1);


    for (int lvl = 0; lvl < output.size(); lvl++){
      LevelData<EBCellFAB>& state  = *m_state[lvl];
      LevelData<EBCellFAB>& source = *m_source[lvl];
      LevelData<EBCellFAB>& flx    = *flux[lvl];
      LevelData<EBCellFAB>& res    = *m_resid[lvl];

      state.localCopyTo(Interval(0,0),          *output[lvl], Interval(0,0));
      flx.localCopyTo(Interval(0,SpaceDim - 1), *output[lvl], Interval(1, SpaceDim));
      source.localCopyTo(Interval(0,0),         *output[lvl], Interval(1 + SpaceDim, 1 + SpaceDim));
      res.localCopyTo(Interval(0,0),            *output[lvl], Interval(2+SpaceDim, 2+SpaceDim));
    }

    // Transform to centroid-centered
    irreg_amr_stencil<centroid_interp>& sten = m_amr->get_centroid_interp_stencils(m_realm, phase::gas);
    sten.apply(output);

    // Alias this stuff
    Vector<LevelData<EBCellFAB>* > output_ptr;
    m_amr->alias(output_ptr, output);

    Vector<Real> covered_values(ncomps, 0.0);
    string fname(file_char);
    writeEBHDF5(fname,
		m_amr->get_grids(m_realm),
		output_ptr,
		names,
		m_amr->get_domains()[0].domainBox(),
		m_amr->get_dx()[0],
		m_dt,
		m_time,
		m_amr->get_ref_rat(),
		m_amr->get_finest_level() + 1,
		true,
		covered_values);
  }

  void eddington_sp1::write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const {
    CH_TIME("eddington_sp1::write_checkpoint_level");
    if(m_verbosity > 5){
      pout() << m_name + "::write_checkpoint_level" << endl;
    }

    // Write state vector
    write(a_handle, *m_state[a_level], m_name);
  }

  void eddington_sp1::read_checkpoint_level(HDF5Handle& a_handle, const int a_level){
    CH_TIME("eddington_sp1::read_checkpoint_level");
    if(m_verbosity > 5){
      pout() << m_name + "::read_checkpoint_level" << endl;
    }

    read<EBCellFAB>(a_handle, *m_state[a_level], m_name, m_amr->get_grids(m_realm)[a_level], Interval(0,0), false);
  }

  void eddington_sp1::set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
    CH_TIME("eddington_sp1::set_neumann_wall_bc");
    if(m_verbosity > 5){
      pout() << "eddington_sp1::set_neumann_wall_bc" << endl;
    }

    const int idx = wall_bc::map_bc(a_dir, a_side);
    m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, wallbc::neumann));
    m_wallbc[idx]->set_value(a_value);
  }

  void eddington_sp1::set_robin_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
    CH_TIME("eddington_sp1::set_robin_wall_bc");
    if(m_verbosity > 5){
      pout() << "eddington_sp1::set_robin_wall_bc" << endl;
    }

    const int idx = wall_bc::map_bc(a_dir, a_side);
    m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, wallbc::robin));
    m_wallbc[idx]->set_value(a_value);
  }
}
