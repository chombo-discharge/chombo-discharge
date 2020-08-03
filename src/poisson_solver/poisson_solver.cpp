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

Real poisson_solver::s_potential_one(const Real a_time){
  return 1.0;
}

Real poisson_solver::s_constant_one(const RealVect a_pos){
  return 1.0;
}

poisson_solver::poisson_solver(){
  m_class_name = "poisson_solver";
  this->set_verbosity(-1);
}

poisson_solver::~poisson_solver(){
  
}

bool poisson_solver::solve(const bool a_zerophi) {
  return this->solve(m_state, m_source, a_zerophi);
}

bool poisson_solver::solve(MFAMRCellData& a_state, const bool a_zerophi){
  return this->solve(a_state, m_source, a_zerophi);
}

void poisson_solver::compute_E(){
  CH_TIME("poisson_solver::compute_E()");
  if(m_verbosity > 5){
    pout() << "poisson_solver::compute_E()" << endl;
  }

  this->compute_E(m_E, m_state);
}

void poisson_solver::compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential){
  CH_TIME("poisson_solver::compute_E(mfamrcell, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "poisson_solver::compute_E(mfamrcell, mfamrcell)" << endl;
  }

  m_amr->compute_gradient(a_E, a_potential);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E);
  m_amr->interp_ghost(a_E);
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
  m_amr->allocate(m_E,      SpaceDim, query_ghost());

  data_ops::set_value(m_state,  0.0);
  data_ops::set_value(m_source, 0.0);
  data_ops::set_value(m_sigma,  0.0);
  data_ops::set_value(m_resid,  0.0);
  data_ops::set_value(m_E,  0.0);
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

void poisson_solver::pre_regrid(const int a_lbase, const int a_old_finest_level){
  CH_TIME("poisson_solver::pre_regrid");
  if(m_verbosity > 5){
    pout() << "poisson_solver::pre_regrid" << endl;
  }

  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  
  m_amr->allocate(m_cache, ncomp);
  for (int lvl = 0; lvl <= a_old_finest_level; lvl++){
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
	const RealVect origin        = m_amr->get_prob_lo();
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

void poisson_solver::regrid(const int a_lmin, const int a_old_finest, const int a_new_finest){
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
      m_amr->alias(scratch_phase, cur_phase, m_cache, Min(a_old_finest, a_new_finest));

      Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->get_eb_pwl_interp(cur_phase);

      // These levels have not changed
      for (int lvl = 0; lvl <= Max(0, a_lmin-1); lvl++){
	scratch_phase[lvl]->copyTo(*state_phase[lvl]); // Base level should never change, but ownership might
      }
      for (int lvl = a_lmin; lvl <= a_new_finest; lvl++){
	interpolator[lvl]->interpolate(*state_phase[lvl], *state_phase[lvl-1], interv);
	if(lvl <= a_old_finest){
	  scratch_phase[lvl]->copyTo(*state_phase[lvl]);
	}
      }
    }
  }

  // Now recompute E
  this->compute_E();
}

void poisson_solver::sanity_check(){
  CH_TIME("poisson_solver::sanity_check");
  if(m_verbosity > 4){
    pout() << "poisson_solver::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());

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

void poisson_solver::set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const potential a_live){
  CH_TIME("poisson_solver::set_dirichlet_wall_bc");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_dirichlet_wall_bc" << endl;
  }

  const int idx = wall_bc::map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, wallbc::dirichlet));
  if(a_live == potential::live){
    m_wallbc[idx]->set_live(true);
  }
  else{
    m_wallbc[idx]->set_live(false);
  }
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
	    const RealVect pos = m_amr->get_prob_lo() + dx*RealVect(iv)*RealVect::Unit;

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

void poisson_solver::set_output_variables(){
  CH_TIME("poisson_solver::set_output_variables");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_output_variables" << endl;
  }

  m_plot_phi = false;
  m_plot_rho = false;
  m_plot_E   = false;
  m_plot_res = false;

  ParmParse pp("poisson_solver");
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi")   m_plot_phi = true;
    else if(str[i] == "rho")   m_plot_rho = true;
    else if(str[i] == "resid") m_plot_res = true;
    else if(str[i] == "E")     m_plot_E = true;
  }
}

void poisson_solver::write_plot_file(){
  CH_TIME("poisson_solver::write_plot_file");
  if(m_verbosity > 5){
    pout() << "poisson_solver::write_plot_file" << endl;
  }

  // Number of output components and their names
  const int ncomps = get_num_plotvars();
  const Vector<std::string> names = get_plotvar_names();

  // Allocate storage for output
  EBAMRCellData output;
  m_amr->allocate(output, phase::gas, ncomps);

  // Copy internal data to be plotted over to 'output'
  int icomp = 0;
  write_plot_data(output, icomp);

  // Filename
  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", "poisson_solver", m_step, SpaceDim);

  // Alias
  Vector<LevelData<EBCellFAB>* > output_ptr(1+m_amr->get_finest_level());
  m_amr->alias(output_ptr, output);

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

void poisson_solver::write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("poisson_solver::write_checkpoint_level");
  if(m_verbosity > 5){
    pout() << "poisson_solver::write_checkpoint_level" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  // Used for aliasing phases
  LevelData<EBCellFAB> state_gas;
  LevelData<EBCellFAB> state_sol;

  if(!ebis_gas.isNull()) mfalias::aliasMF(state_gas,  phase::gas,   *m_state[a_level]);
  if(!ebis_sol.isNull()) mfalias::aliasMF(state_sol,  phase::solid, *m_state[a_level]);
  
  // Write data
  if(!ebis_gas.isNull()) write(a_handle, state_gas, "poisson_g");
  if(!ebis_sol.isNull()) write(a_handle, state_sol, "poisson_s");
}

void poisson_solver::read_checkpoint_level(HDF5Handle& a_handle, const int a_level){
  CH_TIME("poisson_solver::read_checkpoint_level");
  if(m_verbosity > 5){
    pout() << "poisson_solver::read_checkpoint_level" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  // Used for aliasing phases
  LevelData<EBCellFAB> state_gas;
  LevelData<EBCellFAB> state_sol;

  if(!ebis_gas.isNull()) mfalias::aliasMF(state_gas,  phase::gas,   *m_state[a_level]);
  if(!ebis_sol.isNull()) mfalias::aliasMF(state_sol,  phase::solid, *m_state[a_level]);
  
  // Read data
  if(!ebis_gas.isNull()) read<EBCellFAB>(a_handle, state_gas, "poisson_g", m_amr->get_grids()[a_level], Interval(0,0), false);
  if(!ebis_sol.isNull()) read<EBCellFAB>(a_handle, state_sol, "poisson_s", m_amr->get_grids()[a_level], Interval(0,0), false);
}

void poisson_solver::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("poisson_solver::write_plot_level");
  if(m_verbosity > 5){
    pout() << "poisson_solver::write_plot_level" << endl;
  }

  // Add phi to output
  if(m_plot_phi) write_mfdata(a_output, a_comp, m_state,  true);
  if(m_plot_rho) {
    write_mfdata(a_output, a_comp, m_source, false);
    data_ops::set_covered_value(a_output, a_comp-1, 0.0); // We wrote m_source to a_comp and then increment by one. Undo that. 
  }
  if(m_plot_res) write_mfdata(a_output, a_comp, m_resid,  false);
  if(m_plot_E) {
    this->compute_E();
    write_mfdata(a_output, a_comp, m_E, true);
  }
}

void poisson_solver::write_mfdata(EBAMRCellData& a_output, int& a_comp, const MFAMRCellData& a_data, const bool a_interp){
  CH_TIME("poisson_solver::write_mfdata");
  if(m_verbosity > 5){
    pout() << "poisson_solver::write_mfdata" << endl;
  }

  const int ncomp = a_data[0]->nComp();

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  // Allocate some scratch data that we can use
  EBAMRCellData scratch;
  m_amr->allocate(scratch, phase::gas, ncomp);

  //
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    LevelData<EBCellFAB> data_gas;
    LevelData<EBCellFAB> data_sol;

    if(!ebis_gas.isNull()) mfalias::aliasMF(data_gas,  phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()) mfalias::aliasMF(data_sol,  phase::solid, *a_data[lvl]);

    data_gas.localCopyTo(*scratch[lvl]);

    // Copy all covered cells from the other phase
    if(!ebis_sol.isNull()){
      for (DataIterator dit = data_sol.dataIterator(); dit.ok(); ++dit){
	const Box box = data_sol.disjointBoxLayout().get(dit());
	const IntVectSet ivs(box);
	const EBISBox& ebisb_gas = data_gas[dit()].getEBISBox();
	const EBISBox& ebisb_sol = data_sol[dit()].getEBISBox();

	FArrayBox& scratch_gas    = (*scratch[lvl])[dit()].getFArrayBox();
	const FArrayBox& fab_gas = data_gas[dit()].getFArrayBox();
	const FArrayBox& fab_sol = data_sol[dit()].getFArrayBox();

	// TLDR: There are four cases
	// 1. Both are covered => inside electrode
	// 2. Gas is covered, solid is regular => inside solid phase only
	// 3. Gas is regular, solid is covered => inside gas phase only
	// 4. Gas is !(covered || regular) and solid is !(covered || regular) => on dielectric boundary
	if(!(ebisb_gas.isAllCovered() && ebisb_sol.isAllCovered())){ // Outside electrode, lets do stuff
	  if(ebisb_gas.isAllCovered() && ebisb_sol.isAllRegular()){ // Case 2, copy directly. 
	    scratch_gas += fab_sol;
	  }
	  else if(ebisb_gas.isAllRegular() && ebisb_sol.isAllCovered()) { // Case 3. Inside gas phase. Already did this. 
	  }
	  else{ // Case 4, needs special treatment. 
	    for (BoxIterator bit(box); bit.ok(); ++bit){ // Loop through all cells here
	      const IntVect iv = bit();
	      if(ebisb_gas.isCovered(iv) && !ebisb_sol.isCovered(iv)){   // Regular cells from phase 2
		for (int comp = 0; comp < ncomp; comp++){
		  scratch_gas(iv, comp) = fab_sol(iv, comp);
		}
	      }
	      else if(ebisb_sol.isIrregular(iv) && ebisb_gas.isIrregular(iv)){ // Irregular cells. Use gas side data
		for (int comp = 0; comp < ncomp; comp++){
		  scratch_gas(iv, comp) = fab_gas(iv,comp);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Average down shit and interpolate to centroids
  m_amr->average_down(scratch, phase::gas);
  m_amr->interp_ghost(scratch, phase::gas);
  if(a_interp){
    m_amr->interpolate_to_centroids(scratch, phase::gas);
  }

  const Interval src_interv(0, ncomp-1);
  const Interval dst_interv(a_comp, a_comp + ncomp - 1);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }

  a_comp += ncomp;
}

int poisson_solver::get_num_plotvars() const {
  CH_TIME("poisson_solver::get_num_plotvars");
  if(m_verbosity > 5){
    pout() << "poisson_solver::get_num_plotvars" << endl;
  }

  int num_output = 0;

  if(m_plot_phi) num_output = num_output + 1;
  if(m_plot_rho) num_output = num_output + 1;
  if(m_plot_res) num_output = num_output + 1;
  if(m_plot_E)   num_output = num_output + SpaceDim;

  return num_output;
}

Vector<std::string> poisson_solver::get_plotvar_names() const {
  CH_TIME("poisson_solver::get_plotvar_names");
  if(m_verbosity > 5){
    pout() << "poisson_solver::get_plotvar_names" << endl;
  }
  Vector<std::string> names(0);
  
  if(m_plot_phi) names.push_back("potential");
  if(m_plot_rho) names.push_back("Space charge density");
  if(m_plot_res) names.push_back("potential_residue");
  if(m_plot_E){
    names.push_back("x-Electric field"); 
    names.push_back("y-Electric field"); 
    if(SpaceDim == 3){
      names.push_back("z-Electric field");
    }
  }
  
  return names;
}

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

MFAMRCellData& poisson_solver::get_E(){
  return m_E;
}

MFAMRCellData& poisson_solver::get_source(){
  return m_source;
}

MFAMRCellData& poisson_solver::get_resid(){
  return m_resid;
}
