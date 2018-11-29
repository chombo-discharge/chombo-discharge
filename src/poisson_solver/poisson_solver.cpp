/*!
  @file   poisson_solver.cpp
  @brief  Implementation of poisson_solver.H
  @author Robert Marskar
  @date   Nov. 2017
  @todo   set_covered_potential is buggy and breaks when we don't have dielectrics. Fix this. 
*/

#include "poisson_solver.H"
#include "MFAliasFactory.H"
#include "mfalias.H"
#include "data_ops.H"
#include "units.H"

#include <EBArith.H>
#include <iostream>
#include <ParmParse.H>
#include <MFAMRIO.H>
#include <EBAMRIO.H>

Real poisson_solver::s_constant_one(const RealVect a_pos){
  return 1.0;
}

poisson_solver::poisson_solver(){
  this->set_verbosity(-1);
  
  CH_TIME("poisson_solver::poisson_solver");
  if(m_verbosity > 5){
    pout() << "poisson_solver::poisson_solver" << endl;
  }

  this->allocate_wall_bc();
  this->set_time(0, 0., 0.);

  m_autotune = false;
  { // Check if we should use auto-tuning of the solver
    std::string str;
    ParmParse pp("poisson_solver");
    if(pp.contains("auto_tune")){
      pp.get("auto_tune", str);
      if(str == "true"){
	m_autotune = true;
      }
    }
  }

  if(SpaceDim == 2){
    this->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
    this->set_neumann_wall_bc(0,   Side::Hi, 0.0);
    this->set_dirichlet_wall_bc(1, Side::Lo, potential::ground);
    this->set_dirichlet_wall_bc(1, Side::Hi, potential::live);
  }
  else if(SpaceDim == 3){
    this->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
    this->set_neumann_wall_bc(0,   Side::Hi, 0.0);
    this->set_neumann_wall_bc(1,   Side::Lo, 0.0);                  
    this->set_neumann_wall_bc(1,   Side::Hi, 0.0);
    this->set_dirichlet_wall_bc(2, Side::Lo, potential::ground);
    this->set_dirichlet_wall_bc(2, Side::Hi, potential::live);
  }

  { // Get parameters from input script
    ParmParse pp("poisson_solver");

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
	    if(type == "dirichlet_ground"){
	      this->set_dirichlet_wall_bc(dir, Side::Lo, potential::ground);
	    }
	    else if(type == "dirichlet_live"){
	      this->set_dirichlet_wall_bc(dir, Side::Lo, potential::live);
	    }
	    else if(type == "neumann"){
	      this->set_neumann_wall_bc(dir, Side::Lo, 0.0);
	    }
	    else if(type == "robin"){
	      this->set_robin_wall_bc(dir, Side::Lo, 0.0);
	    }
	    else {
	      std::string error = "poisson_solver::poisson_solver - unknown bc requested for " + bc_string;
	      MayDay::Abort(error.c_str());
	    }
	  }
	}
	else if(side == Side::Hi){
	  std::string type;
	  std::string bc_string = "bc_" + str_dir + "_high";
	  if(pp.contains(bc_string.c_str())){
	    pp.get(bc_string.c_str(), type);
	    if(type == "dirichlet_ground"){
	      this->set_dirichlet_wall_bc(dir, Side::Hi, potential::ground);
	    }
	    else if(type == "dirichlet_live"){
	      this->set_dirichlet_wall_bc(dir, Side::Hi, potential::live);
	    }
	    else if(type == "neumann"){
	      this->set_neumann_wall_bc(dir, Side::Hi, 0.0);
	    }
	    else if(type == "robin"){
	      this->set_robin_wall_bc(dir, Side::Hi, 0.0);
	    }
	    else {
	      std::string error = "poisson_solver::poisson_solver - unknown bc requested for " + bc_string;
	      MayDay::Abort(error.c_str());
	    }
	  }
	}
      }
    }
  }

  // Set default distribution on domain edges(faces)
  m_wall_func_x_lo = poisson_solver::s_constant_one;
  m_wall_func_x_hi = poisson_solver::s_constant_one;
  m_wall_func_y_lo = poisson_solver::s_constant_one;
  m_wall_func_y_hi = poisson_solver::s_constant_one;
#if CH_SPACEDIM==3
  m_wall_func_z_lo = poisson_solver::s_constant_one;
  m_wall_func_z_hi = poisson_solver::s_constant_one;
#endif
}

poisson_solver::~poisson_solver(){
  
}

bool poisson_solver::solve(const bool a_zerophi) {
  this->solve(m_state, m_source, a_zerophi);
}

bool poisson_solver::solve(MFAMRCellData& a_state, const bool a_zerophi){
  this->solve(a_state, m_source, a_zerophi);
}

void poisson_solver::allocate_internals(){
  CH_TIME("poisson_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << "poisson_solver::allocate_internals" << endl;
  }

  const int ncomp = 1;
  
  m_amr->allocate(m_state,  ncomp, query_ghost());
  m_amr->allocate(m_source, ncomp, query_ghost());
  m_amr->allocate(m_resid,  ncomp, query_ghost());
  m_amr->allocate(m_sigma,  phase::gas, ncomp, query_ghost());

  data_ops::set_value(m_state,  0.0);
  data_ops::set_value(m_source, 0.0);
  data_ops::set_value(m_sigma,  0.0);
  data_ops::set_value(m_resid,  0.0);
}

void poisson_solver::allocate_wall_bc(){
  CH_TIME("poisson_solver::poisson_solver(full)");
  if(m_verbosity > 5){
    pout() << "poisson_solver::poisson_solver(full)" << endl;
  }
  m_wallbc.resize(2*SpaceDim);
  for (int i = 0; i < 2*SpaceDim; i++){
    m_wallbc[i] = RefCountedPtr<wall_bc> (NULL);
  }
}

void poisson_solver::cache_state(){
  CH_TIME("poisson_solver::cache_state");
  if(m_verbosity > 5){
    pout() << "poisson_solver::cache_state" << endl;
  }

  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  
  m_amr->allocate(m_cache, ncomp);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    m_state[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void poisson_solver::compute_D(MFAMRCellData& a_D, const MFAMRCellData& a_E){
  CH_TIME("poisson_solver::compute_D");
  if(m_verbosity > 5){
    pout() << "poisson_solver::compute_D" << endl;
  }

  const Vector<dielectric>& dielectrics = m_compgeom->get_dielectrics();

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    LevelData<MFCellFAB>& D = *a_D[lvl];
    LevelData<MFCellFAB>& E = *a_E[lvl];

    LevelData<EBCellFAB> D_gas, D_solid;
    LevelData<EBCellFAB> E_gas, E_solid;

    // This is all we need for the gas phase
    mfalias::aliasMF(D_gas,   phase::gas, D);
    mfalias::aliasMF(E_gas,   phase::gas, E);
    E_gas.localCopyTo(D_gas);
    data_ops::scale(D_gas,   units::s_eps0);
    
    if(m_mfis->num_phases() > 1){
      mfalias::aliasMF(D_solid, phase::solid, D);
      mfalias::aliasMF(E_solid, phase::solid, E);
      E_solid.localCopyTo(D_solid);
      data_ops::scale(D_solid, units::s_eps0);

      // Now scale by relative epsilon
      if(dielectrics.size() > 0){
	const RealVect dx            = m_amr->get_dx()[lvl]*RealVect::Unit;
	const RealVect origin        = m_physdom->get_prob_lo();
	const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  EBCellFAB& dg = D_gas[dit()];

	  const Box box = dbl.get(dit());
	  for (VoFIterator vofit(IntVectSet(box), dg.getEBISBox().getEBGraph()); vofit.ok(); ++vofit){
	    const VolIndex& vof = vofit();
	    const RealVect& pos = EBArith::getVofLocation(vof, origin, dx);

	    Real dist = 1.E99;
	    int closest = 0;
	    for (int i = 0; i < dielectrics.size(); i++){
	      const RefCountedPtr<BaseIF> func = dielectrics[i].get_function();

	      const Real cur_dist = func->value(pos);
	
	      if(cur_dist <= dist){
		dist = cur_dist;
		closest = i;
	      }
	    }
	    const Real eps = dielectrics[closest].get_permittivity(pos);

	    for (int comp = 0; comp < dg.nComp(); comp++){
	      dg(vof, comp) *= eps;
	    }
	  }
	}
      }
    }
  }
}

Real poisson_solver::compute_U(const MFAMRCellData& a_E){
  CH_TIME("poisson_solver::compute_U");
  if(m_verbosity > 5){
    pout() << "poisson_solver::compute_U" << endl;
  }

  MFAMRCellData D, EdotD;   
  m_amr->allocate(D, SpaceDim);
  m_amr->allocate(EdotD, 1);
  this->compute_D(D, a_E);
  data_ops::dot_prod(EdotD, D, a_E);

  Real U_g = 0.0;
  Real U_s = 0.0;

  // Energy in gas phase
  EBAMRCellData data_g;
  m_amr->allocate_ptr(data_g);
  m_amr->alias(data_g, phase::gas, EdotD);
  m_amr->average_down(data_g, phase::gas);
  data_ops::norm(U_g, *data_g[0], m_amr->get_domains()[0], 1);

  if(m_mfis->num_phases() > 1){
    EBAMRCellData data_s;
    m_amr->allocate_ptr(data_s);
    m_amr->alias(data_s, phase::solid, EdotD);
    m_amr->average_down(data_s, phase::solid);
    data_ops::norm(U_s, *data_s[0], m_amr->get_domains()[0], 1);
  }

  return 0.5*(U_g + U_s);
}

Real poisson_solver::compute_capacitance(){
  CH_TIME("poisson_solver::compute_capacitance");
  if(m_verbosity > 5){
    pout() << "poisson_solver::compute_capacitance" << endl;
  }

  // TLDR; We MUST compute the energy density with the Laplace field, so no sources here...
  Real C;

  MFAMRCellData phi, source;
  EBAMRIVData sigma;

  m_amr->allocate(phi, 1);
  m_amr->allocate(source, 1);
  m_amr->allocate(sigma, phase::gas, 1);

  data_ops::set_value(phi,    0.0);
  data_ops::set_value(source, 0.0);
  data_ops::set_value(sigma,  0.0);

  solve(phi, source, sigma);

  // Solve and compute energy density
  MFAMRCellData E;
  m_amr->allocate(E, SpaceDim);
  m_amr->compute_gradient(E, phi); // -E
  const Real U = this->compute_U(E); // Energy density

  // U = 0.5*CV^2
  const Real pot = m_potential(m_time);
  if(pot == 0.0){
    MayDay::Abort("poisson_solver::compute_capacitance - error, can't compute energy density with V = 0");
  }
  C = 2.0*U/(pot*pot);

  return C;
}

void poisson_solver::deallocate_internals(){
  CH_TIME("poisson_solver::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "poisson_solver::deallocate_internals" << endl;
  }
  m_amr->deallocate(m_state);
  m_amr->deallocate(m_source);
  m_amr->deallocate(m_resid);
  m_amr->deallocate(m_sigma);
}

void poisson_solver::regrid(const int a_old_finest, const int a_new_finest){
  CH_TIME("poisson_solver::regrid");
  if(m_verbosity > 5){
    pout() << "poisson_solver::regrid" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocate_internals();
  
  for (int i = 0; i < phase::num_phases; i++){
    phase::which_phase cur_phase;    
    if(i == 0){
      cur_phase = phase::gas;
    }
    else{
      cur_phase = phase::solid;
    }

    const RefCountedPtr<EBIndexSpace>& ebis = m_mfis->get_ebis(cur_phase);

    if(!ebis.isNull()){
    
      EBAMRCellData scratch_phase;
      EBAMRCellData state_phase;

      m_amr->allocate_ptr(scratch_phase);
      m_amr->allocate_ptr(state_phase);
      m_amr->alias(state_phase,   cur_phase, m_state);
      m_amr->alias(scratch_phase, cur_phase, m_cache, a_old_finest);

      Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->get_eb_pwl_interp(cur_phase);

      scratch_phase[0]->copyTo(*state_phase[0]); // Base level should never change, but ownership might
      for (int lvl = 1; lvl <= a_new_finest; lvl++){
	interpolator[lvl]->interpolate(*state_phase[lvl], *state_phase[lvl-1], interv);
	if(lvl <= a_old_finest){
	  scratch_phase[lvl]->copyTo(*state_phase[lvl]);
	}
      }
    }
  }
}

void poisson_solver::sanity_check(){
  CH_TIME("poisson_solver::sanity_check");
  if(m_verbosity > 4){
    pout() << "poisson_solver::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_physdom.isNull());

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sideit; sideit.ok(); ++sideit){
      if(m_wallbc[wall_bc::map_bc(dir, sideit())].isNull()){
	pout() << "poisson_solver::sanity_check() - bc is null at coord = " << dir << ", side = " << sideit() << endl;
  	MayDay::Abort("poisson_solver::sanity_check() failed. Wall BC has not been set properly");
      }
    }
  }
}

void poisson_solver::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("poisson_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_computational_geometry" << endl;
  }

  m_compgeom = a_compgeom;

  this->set_mfis(m_compgeom->get_mfis());
}

void poisson_solver::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("poisson_solver::set_physical_domain");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_physical_domain" << endl;
  }

  m_physdom = a_physdom;
}

void poisson_solver::set_mfis(const RefCountedPtr<mfis>& a_mfis){
  CH_TIME("poisson_solver::set_mfis");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_mfis" << endl;
  }

  m_mfis = a_mfis;
}

void poisson_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("poisson_solver::set_amr");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_amr" << endl;
  }

  m_amr = a_amr;
}

void poisson_solver::set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const potential::ground_live a_live){
  CH_TIME("poisson_solver::set_dirichlet_wall_bc");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_dirichlet_wall_bc" << endl;
  }

  const int idx = wall_bc::map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, wallbc::dirichlet));
  m_wallbc[idx]->set_live(a_live);
}

void poisson_solver::set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
  CH_TIME("poisson_solver::set_neumann_wall_bc");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_neumann_wall_bc" << endl;
  }

  const int idx = wall_bc::map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, wallbc::neumann));
  m_wallbc[idx]->set_value(a_value);
}

void poisson_solver::set_robin_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
  CH_TIME("poisson_solver::set_robin_wall_bc");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_robin_wall_bc" << endl;
  }

  const int idx = wall_bc::map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, wallbc::robin));
  m_wallbc[idx]->set_value(a_value);
}

void poisson_solver::set_potential(Real (*a_potential)(const Real a_time)){
  m_potential = a_potential;
}

void poisson_solver::set_poisson_wall_func(const int a_dir, const Side::LoHiSide a_side, Real (*a_func)(const RealVect a_pos)){
  CH_TIME("poisson_solver::set_poisson_wall_func");
  if(m_verbosity > 4){
    pout() << "poisson_solver::set_poisson_wall_func" << endl;
  }

  if(a_dir == 0){
    if(a_side == Side::Lo){
      m_wall_func_x_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_x_hi = a_func;
    }
  }
  else if(a_dir == 1){
    if(a_side == Side::Lo){
      m_wall_func_y_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_y_hi = a_func;
    }
  }
#if CH_SPACEDIM==3
  else if(a_dir == 2){
    if(a_side == Side::Lo){
      m_wall_func_z_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_z_hi = a_func;
    }
  }
#endif
}

void poisson_solver::set_verbosity(const int a_verbosity){
  CH_TIME("poisson_solver::set_verbosity");
  m_verbosity = a_verbosity;

  if(m_verbosity > 4){
    pout() << "poisson_solver::set_verbosity" << endl;
  }
}

void poisson_solver::set_time(const int a_step, const Real a_time, const Real a_dt) {
  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void poisson_solver::set_covered_potential(EBAMRCellData& a_phi, const int a_comp, const Real a_time){
  CH_TIME("poisson_solver::set_covered_potential");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_covered_potential" << endl;
  }


  const Vector<electrode>& electrodes = m_compgeom->get_electrodes();

  if(electrodes.size() > 0){
    const int finest_level = m_amr->get_finest_level();
    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl_gas  = m_amr->get_ebisl(phase::gas)[lvl];
      const EBISLayout& ebisl_sol  = m_amr->get_ebisl(phase::solid)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	EBCellFAB& phi = (*a_phi[lvl])[dit()];
	const EBISBox& ebisbox_gas = ebisl_gas[dit()];
	const EBISBox& ebisbox_sol = ebisl_sol[dit()];

	if(!ebisbox_gas.isAllRegular() && !ebisbox_sol.isAllRegular()){
	  const Box box = dbl.get(dit());
	  FArrayBox& fbox = (*a_phi[lvl])[dit()].getFArrayBox();

	  for (BoxIterator bit(box); bit.ok(); ++bit){
	    const IntVect& iv  = bit();
	    const RealVect pos = m_physdom->get_prob_lo() + dx*RealVect(iv)*RealVect::Unit;

	    if(ebisbox_gas.isCovered(iv) && ebisbox_sol.isCovered(iv)){
	      int closest;
	      Real dist = 1.E99;
	      for (int i = 0; i < electrodes.size(); i++){
		const Real cur_dist = electrodes[i].get_function()->value(pos);
		if(Abs(cur_dist) < dist){
		  dist = Abs(cur_dist);
		  closest = i;
		}
	      }

	      if(electrodes[closest].is_live()){
		fbox(iv, a_comp) = electrodes[closest].get_fraction()*m_potential(a_time);
	      }
	      else{
		fbox(iv, a_comp) = 0.0;
	      }
	    }
	  }
	}
      }
    }
  }
}

#ifdef CH_USE_HDF5
void poisson_solver::write_plot_file(){
  CH_TIME("poisson_solver::write_plot_file");
  if(m_verbosity > 5){
    pout() << "poisson_solver::write_plot_file" << endl;
  }

  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", "poisson_solver", m_step, SpaceDim);

  const int ncomps = 3 + SpaceDim;
  Vector<string> names(ncomps);
  names[0] = "potential";
  names[1] = "space charge density";
  names[2] = "residue";
  names[3] = "x-Electric field";
  names[4] = "y-Electric field";
  if(SpaceDim == 3){
    names[5] = "z-Electric field";
  }

  // Compute the electric field. 
  MFAMRCellData E;
  m_amr->allocate(E, SpaceDim, 3);
  m_amr->compute_gradient(E, m_state);
  m_amr->interp_ghost(E);
  m_amr->average_down(E);
  data_ops::scale(E, -1.0);

  Vector<RefCountedPtr<LevelData<EBCellFAB> > > output;
  m_amr->allocate(output, phase::gas, ncomps, 1);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);



  for (int lvl = 0; lvl < output.size(); lvl++){
    LevelData<EBCellFAB> state_gas,  source_gas, E_gas, resid_gas;
    LevelData<EBCellFAB> state_sol,  source_sol, E_sol, resid_sol;

    mfalias::aliasMF(state_gas,  phase::gas,   *m_state[lvl]);
    mfalias::aliasMF(source_gas, phase::gas,   *m_source[lvl]);
    mfalias::aliasMF(E_gas,      phase::gas,   *E[lvl]);
    mfalias::aliasMF(resid_gas,  phase::gas,   *m_resid[lvl]);

    if(!ebis_sol.isNull()){
      mfalias::aliasMF(state_sol,  phase::solid, *m_state[lvl]);
      mfalias::aliasMF(source_sol, phase::solid, *m_source[lvl]);
      mfalias::aliasMF(E_sol,      phase::solid, *E[lvl]);
      mfalias::aliasMF(resid_sol,  phase::solid, *m_resid[lvl]);
    }



    // Copy all covered cells from the other phase
#if 0 // This code segment is really shit - it changes the potential at irregular cells which it really shouldn't do...
    // Fix this if you really need this functionality. 
    if(!ebis_sol.isNull()){
      for (DataIterator dit = state_sol.dataIterator(); dit.ok(); ++dit){
	const Box box = state_sol.disjointBoxLayout().get(dit());
	const IntVectSet ivs(box);
	const EBISBox& ebisb_gas = state_gas[dit()].getEBISBox();
	const EBISBox& ebisb_sol = state_sol[dit()].getEBISBox();

	FArrayBox& data_gas = state_gas[dit()].getFArrayBox();
	FArrayBox& src_gas  = source_gas[dit()].getFArrayBox();
	FArrayBox& e_gas    = E_gas[dit()].getFArrayBox();
	FArrayBox& res_gas  = resid_gas[dit()].getFArrayBox();

	FArrayBox& data_sol = state_sol[dit()].getFArrayBox();
	FArrayBox& src_sol  = source_sol[dit()].getFArrayBox();
	FArrayBox& e_sol    = E_sol[dit()].getFArrayBox();
	FArrayBox& res_sol  = resid_sol[dit()].getFArrayBox();

	for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit){
	  const IntVect iv = ivsit();
	  if(ebisb_gas.isCovered(iv) && !ebisb_sol.isCovered(iv)){ // Regular cells from phase 2
	    data_gas(iv, 0) = data_sol(iv,0);
	    src_gas(iv,0)   = src_sol(iv,0);
	    res_gas(iv, 0)  = res_sol(iv, 0);
	    for (int comp = 0; comp < SpaceDim; comp++){
	      e_gas(iv,comp) = e_sol(iv, comp);
	    }
	  }
	  if(ebisb_sol.isIrregular(iv) && ebisb_gas.isIrregular(iv)){ // Irregular cells
	    data_gas(iv, 0) = 0.5*(data_gas(iv,0) + data_sol(iv,0));
	    src_gas(iv, 0)  = 0.5*(src_gas(iv,0) + src_sol(iv,0));
	    res_gas(iv, 0)  = 0.5*(res_gas(iv,0) + res_sol(iv,0));
	  }
	}
      }
    }
#endif

    state_gas.localCopyTo(Interval(0,0),        *output[lvl], Interval(0,0));
    source_gas.localCopyTo(Interval(0,0),       *output[lvl], Interval(1,1));
    resid_gas.localCopyTo(Interval(0,0),        *output[lvl], Interval(2,2));
    E_gas.localCopyTo(Interval(0,SpaceDim - 1), *output[lvl], Interval(3, 2 + SpaceDim));
  }

  const irreg_amr_stencil<centroid_interp>& sten = m_amr->get_centroid_interp_stencils(phase::gas);
  sten.apply(output);

  Vector<LevelData<EBCellFAB>* > output_ptr;
  for (int lvl = 0; lvl < m_state.size(); lvl++){
    output_ptr.push_back(&(*output[lvl]));
  }

  Vector<Real> covered_values(ncomps, 0.0);
  string fname(file_char);
  writeEBHDF5(fname,
	      m_amr->get_grids(),
	      output_ptr,
	      names,
	      m_amr->get_domains()[0].domainBox(),
	      m_amr->get_dx()[0],
	      m_dt,
	      m_time,
	      m_amr->get_ref_rat(),
	      m_amr->get_finest_level() + 1,
	      false,
	      covered_values,
	      IntVect::Unit);
}
#endif

Real poisson_solver::get_time() const{
  return m_time;
}


wall_bc& poisson_solver::get_wall_bc(const int a_dir, Side::LoHiSide a_side) const{
  CH_TIME("poisson_solver::get_wall_bc");
  if(m_verbosity > 5){
    pout() << "poisson_solver::get_wall_bc" << endl;
  }
  return *m_wallbc[wall_bc::map_bc(a_dir, a_side)];
}

MFAMRCellData& poisson_solver::get_state(){
  return m_state;
}

MFAMRCellData& poisson_solver::get_source(){
  return m_source;
}

MFAMRCellData& poisson_solver::get_resid(){
  return m_resid;
}
