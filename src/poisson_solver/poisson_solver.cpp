/*!
  @file poisson_solver.cpp
  @brief Implementation of poisson_solver.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_solver.H"
#include "MFAliasFactory.H"
#include "mfalias.H"
#include "data_ops.H" 

#include <MFAMRIO.H>
#include <EBAMRIO.H>

poisson_solver::poisson_solver(){
  CH_TIME("poisson_solver::set_verbosity");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_verbosity" << endl;
  }
  
  this->set_verbosity(-1);
  this->allocate_wall_bc();
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

void poisson_solver::set_verbosity(const int a_verbosity){
  m_verbosity = a_verbosity;
}

void poisson_solver::set_time(const Real a_time) {
  m_time = a_time;
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

#ifdef CH_USE_HDF5
void poisson_solver::write_plot_file(const int a_step){
  CH_TIME("poisson_solver::write_plot_file");
  if(m_verbosity > 5){
    pout() << "poisson_solver::write_plot_file" << endl;
  }

  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", "poisson_solver", a_step, SpaceDim);

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

  irreg_amr_stencil<centroid_interp>& sten = m_amr->get_centroid_interp_stencils(phase::gas);

  for (int lvl = 0; lvl < output.size(); lvl++){
    LevelData<EBCellFAB> state_gas,  source_gas, E_gas, resid_gas;
    LevelData<EBCellFAB> state_sol,  source_sol, E_sol, resid_sol;

    mfalias::aliasMF(state_gas,  phase::gas,   *m_state[lvl]);
    mfalias::aliasMF(source_gas, phase::gas,   *m_source[lvl]);
    mfalias::aliasMF(E_gas,      phase::gas,   *E[lvl]);
    mfalias::aliasMF(resid_gas,  phase::gas,   *m_resid[lvl]);


#if 1 // Transform to centroid-centered for irregular cells. 
    sten.apply(state_gas, lvl);
    sten.apply(E_gas, lvl);
#endif

    if(!ebis_sol.isNull()){
      mfalias::aliasMF(state_sol,  phase::solid, *m_state[lvl]);
      mfalias::aliasMF(source_sol, phase::solid, *m_source[lvl]);
      mfalias::aliasMF(E_sol,      phase::solid, *E[lvl]);
      mfalias::aliasMF(resid_sol,  phase::solid, *m_resid[lvl]);
    }

    


    // Copy all covered cells from the other phase
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

    state_gas.copyTo(Interval(0,0),        *output[lvl], Interval(0,0));
    source_gas.copyTo(Interval(0,0),       *output[lvl], Interval(1,1));
    resid_gas.copyTo(Interval(0,0),        *output[lvl], Interval(2,2));
    E_gas.copyTo(Interval(0,SpaceDim - 1), *output[lvl], Interval(3, 2 + SpaceDim));
  }

#if 0 // Can only do up to SpaceDim, need to rewrite this routine if we want this to work. 
  m_amr->average_down(output, phase::gas);
  m_amr->interp_ghost(output, phase::gas);
#endif

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
	      0.0,
	      0.0,
	      m_amr->get_ref_rat(),
	      m_amr->get_finest_level() + 1,
	      false,
	      covered_values);
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
