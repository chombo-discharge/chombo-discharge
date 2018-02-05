/*!
  @file time_stepper.cpp
  @brief Implementation of time_stepper.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "time_stepper.H"
#include "poisson_multifluid_gmg.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "units.H"

time_stepper::time_stepper(){
  this->set_verbosity(10);
  this->set_cfl(0.8);
  this->set_relax_time(2.0);
  this->set_min_dt(0.0);
  this->set_max_dt(1.E99);
}

time_stepper::~time_stepper(){
}

int time_stepper::query_ghost(){
  return 3;
}

bool time_stepper::stationary_rte(){
  CH_TIME("time_stepper::stationary_rte");
  if(m_verbosity > 5){
    pout() << "time_stepper::stationary_rte" << endl;
  }

  return m_rte->is_stationary();
}

void time_stepper::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("time_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_dt" << endl;
  }

  Real dt = 1.E99;

  const Real dt_cfl = m_cfl*m_cdr->compute_cfl_dt();
  if(dt_cfl < dt){
    dt = dt_cfl;
    a_timecode = time_code::cfl;
  }
  pout() << dt << endl;
  
  const Real dt_dif = m_cfl*m_cdr->compute_diffusive_dt();
  if(dt_dif < dt){
    dt = dt_dif;
    a_timecode = time_code::diffusion;
  }
  pout() << dt << endl;

  const Real dt_relax = m_relax_time*this->compute_relaxation_time();
  if(dt_relax < dt){
    pout() << dt_relax << endl;
    dt = dt_relax;
    a_timecode = time_code::relaxation_time;
  }
  pout() << dt << endl;

  const Real dt_restrict = this->restrict_dt();
  if(dt_restrict < dt){
    dt = dt_restrict;
    a_timecode = time_code::restricted;
  }
  pout() << dt << endl;

  if(dt < m_min_dt){
    dt = m_min_dt;
    a_timecode = time_code::hardcap;
  }
  pout() << dt << endl;

  if(dt > m_max_dt){
    dt = m_max_dt;
    a_timecode = time_code::hardcap;
  }
  pout() << dt << endl;

  CH_assert(dt > 0.0);
  a_dt = dt;
}

void time_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential){
  CH_TIME("time_stepper::compute_E(ebamrcell, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_E(ebamrcell, phase, mfamrcell)" << endl;
  }

  EBAMRCellData pot_gas;
  m_amr->allocate_ptr(pot_gas);
  m_amr->alias(pot_gas, a_phase, a_potential);

  m_amr->compute_gradient(a_E, pot_gas);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E, a_phase);
  m_amr->interp_ghost(a_E, a_phase);
}

void time_stepper::compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential){
  CH_TIME("time_stepper::compute_E(mfamrcell, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_E(mfamrcell, mfamrcell)" << endl;
  }

  m_amr->compute_gradient(a_E, a_potential);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E);
  m_amr->interp_ghost(a_E);
}

void time_stepper::compute_J(EBAMRCellData& a_J){
  CH_TIME("time_stepper::compute_J");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_J" << endl;
  }

  const int density_comp = 0;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_cdr->get_phase())[lvl];

    data_ops::set_value(*a_J[lvl], 0.0);

    for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver> solver = solver_it();
      RefCountedPtr<species> spec      = solver_it.get_species();

      const int q                      = spec->get_charge();
      const EBAMRCellData& density     = solver->get_state();
      const EBAMRCellData& velo        = solver->get_velo_cell();

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box& box         = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet ivs(box);

	EBCellFAB& J       = (*a_J[lvl])[dit()];
	const EBCellFAB& n = (*density[lvl])[dit()];
	const EBCellFAB& v = (*velo[lvl])[dit()];

	for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();

	  for (int comp = 0; comp < SpaceDim; comp++){
	    J(vof, comp) += q*n(vof,density_comp)*v(vof, comp);
	  }
	}
      }
    }

    data_ops::scale(*a_J[lvl], units::s_Qe);
  }

  m_amr->average_down(a_J, m_cdr->get_phase());
  m_amr->interp_ghost(a_J, m_cdr->get_phase());

#if 1 // override
  data_ops::set_value(a_J, 1.0);
#endif
}

void time_stepper::compute_photon_source_terms(Vector<EBAMRCellData*>        a_source,
					       const Vector<EBAMRCellData*>& a_cdr_states,
					       const EBAMRCellData&          a_E,
					       const centering::which_center a_centering){
  CH_TIME("time_stepper::compute_photon_source_terms(full)");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_photon_source_terms(full)" << endl;
  }

  const phase::which_phase rte_phase = m_rte->get_phase();
  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  const irreg_amr_stencil<centroid_interp>& interp_stencils = m_amr->get_centroid_interp_stencils(rte_phase);

  Vector<Real> cdr_densities(1 + m_plaskin->get_num_species());

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl  = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl       = m_amr->get_ebisl(rte_phase)[lvl];
    const Real dx                 = m_amr->get_dx()[lvl];
    const irreg_stencil& stencils = interp_stencils[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const EBCellFAB& E     = (*a_E[lvl])[dit()];


      // Do all cells
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();

	const RealVect e = RealVect(D_DECL(E(vof, 0), E(vof, 1), E(vof, 2)));
	for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
	  const int idx = solver_it.get_solver();
	  cdr_densities[idx] = (*(*a_cdr_states[idx])[lvl])[dit](vof, comp);
	}

	Vector<Real> sources = m_plaskin->compute_rte_source_terms(cdr_densities, e);
	for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
	  const int idx = solver_it.get_solver();
	  (*(*a_source[idx])[lvl])[dit()](vof, comp) = sources[idx];
	}
      }

      // Do irregular cells
      if(a_centering == centering::cell_center){
	ivs = ebisbox.getIrregIVS(box);
	for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();

	  // Reset these and apply stencil
	  RealVect e = RealVect::Zero;
	  cdr_densities.assign(0.0);
	
	  const VoFStencil& stencil = stencils[dit()](vof, comp);
	  for (int i = 0; i < stencil.size(); i++){
	    const VolIndex& ivof = stencil.vof(i);
	    const Real iweight   = stencil.weight(i);

	    for (int dir = 0; dir < SpaceDim; dir++){
	      e[dir] = E(ivof, dir)*iweight;
	    }
	  
	    for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
	      const int idx = solver_it.get_solver();
	      cdr_densities[idx] = +(*(*a_cdr_states[idx])[lvl])[dit](ivof, comp)*iweight;
	    }
	  }

	  Vector<Real> sources = m_plaskin->compute_rte_source_terms(cdr_densities, e);
	  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    (*(*a_source[idx])[lvl])[dit()](vof, comp) = sources[idx];
	  }
	}
      }
    }
  }

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_amr->average_down(*a_source[idx], rte_phase);
    m_amr->interp_ghost(*a_source[idx], rte_phase);
  }
}

void time_stepper::compute_rho(){
  CH_TIME("time_stepper::compute_rho()");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_rho()" << endl;
  }

  Vector<EBAMRCellData*> densities;
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver> solver = solver_it();
    densities.push_back(&(solver->get_state()));
  }

  this->compute_rho(m_poisson->get_source(),
		    densities,
		    centering::cell_center);
}

void time_stepper::compute_rho(MFAMRCellData&                a_rho,
			       const Vector<EBAMRCellData*>  a_densities,
			       const centering::which_center a_centering){
  CH_TIME("time_stepper::compute_rho(mfamrcell, vec(ebamrcell))");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_rho(mfamrcell, vec(ebamrcell))" << endl;
  }

  data_ops::set_value(a_rho, 0.0);

  EBAMRCellData rho_gas;
  m_amr->allocate_ptr(rho_gas);
  m_amr->alias(rho_gas, phase::gas, a_rho);

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(rho_gas, 0.0);

    // Add volumetric charge 
    for (cdr_iterator solver_it(*m_cdr, a_densities); solver_it.ok(); ++solver_it){
      const EBAMRCellData& density       = solver_it.get_data();
      const RefCountedPtr<species>& spec = solver_it.get_species();

      data_ops::incr(*rho_gas[lvl], *density[lvl], spec->get_charge());
    }

    // Scale by s_Qe/s_eps0
    data_ops::scale(*a_rho[lvl], units::s_Qe);
    data_ops::kappa_scale(*a_rho[lvl]);
  }

  // Transform to centroids
  if(a_centering == centering::cell_center){
    m_amr->interpolate_to_centroids(rho_gas, phase::gas);
  }


}

void time_stepper::instantiate_solvers(){
  CH_TIME("time_stepper::instantiate_solvers");
  if(m_verbosity > 5){
    pout() << "time_stepper::instantiate_solvers" << endl;
  }

  this->sanity_check();

  this->setup_cdr();
  this->setup_rte();
  this->setup_poisson();
  this->setup_sigma();
}

void time_stepper::initial_data(){
  CH_TIME("time_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "time_stepper::initial_data" << endl;
  }

  m_cdr->initial_data();
  if(!m_rte->is_stationary()){
    m_rte->initial_data();
  }
  m_sigma->initial_data();
}

void time_stepper::regrid(const int a_old_finest, const int a_new_finest){
  CH_TIME("time_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "time_stepper::regrid" << endl;
  }

  this->regrid_solvers(a_old_finest, a_new_finest);
  this->regrid_internals();
}

void time_stepper::regrid_solvers(const int a_old_finest, const int a_new_finest){
  CH_TIME("time_stepper::regrid_solvers");
  if(m_verbosity > 5){
    pout() << "time_stepper::regrid_solvers" << endl;
  }

  m_cdr->regrid(a_old_finest,     a_new_finest);
  m_poisson->regrid(a_old_finest, a_new_finest);
  m_rte->regrid(a_old_finest,     a_new_finest);
  m_sigma->regrid(a_old_finest,   a_new_finest);

  m_poisson->write_plot_file();
}

void time_stepper::sanity_check(){
  CH_TIME("time_stepper::sanity_check");
  if(m_verbosity > 5){
    pout() << "time_stepper::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_physdom.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_plaskin.isNull());
}

void time_stepper::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("time_stepper::set_amr");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_amr" << endl;
  }

  m_amr = a_amr;
}

void time_stepper::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("time_stepper::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_computational_geometry" << endl;
  }

  m_compgeom = a_compgeom;
}

void time_stepper::set_plasma_kinetics(const RefCountedPtr<plasma_kinetics>& a_plaskin){
  CH_TIME("time_stepper::set_plasma_kinetics");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_plasma_kinetics" << endl;
  }

  m_plaskin = a_plaskin;
}

void time_stepper::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("time_stepper::set_physical_domain");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_physical_domain" << endl;
  }

  m_physdom = a_physdom;
}

void time_stepper::set_potential(Real (*a_potential)(const Real a_time)){
  CH_TIME("time_stepper::set_potential");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_potential" << endl;
  }
  m_potential     = a_potential;
}

void time_stepper::set_verbosity(const int a_verbosity){
  CH_TIME("time_stepper::set_verbosity");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_verbosity" << endl;
  }
  
  m_verbosity = a_verbosity;


}

void time_stepper::set_solver_verbosity(const int a_verbosity){
  CH_TIME("time_stepper::set_solver_verbosity");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_solver_verbosity" << endl;
  }

  if(!m_cdr.isNull()){
    m_cdr->set_verbosity(a_verbosity);
  }
  if(!m_poisson.isNull()){
    m_poisson->set_verbosity(a_verbosity);
  }
  if(!m_rte.isNull()){
    m_rte->set_verbosity(a_verbosity);
  }
  if(!m_sigma.isNull()){
    m_sigma->set_verbosity(a_verbosity);
  }
}

void time_stepper::set_min_dt(const Real a_min_dt){
  m_min_dt = Max(a_min_dt, 0.0);
}

void time_stepper::set_max_dt(const Real a_max_dt){
  m_max_dt = Max(a_max_dt, 0.0);
}

void time_stepper::set_cfl(const Real a_cfl){
  m_cfl = a_cfl;
}

void time_stepper::set_relax_time(const Real a_relax_time){
  m_relax_time = a_relax_time;
}

void time_stepper::setup_cdr(){
  CH_TIME("time_stepper::setup_cdr");
  if(m_verbosity > 5){
    pout() << "time_stepper::setup_cdr" << endl;
  }

  m_cdr = RefCountedPtr<cdr_layout> (new cdr_layout(m_plaskin));
  m_cdr->set_verbosity(m_verbosity);
  m_cdr->set_amr(m_amr);
  m_cdr->set_computational_geometry(m_compgeom);
  m_cdr->set_physical_domain(m_physdom);
  m_cdr->sanity_check();
  m_cdr->allocate_internals();
}

void time_stepper::setup_poisson(){
  CH_TIME("time_stepper::setup_poisson");
  if(m_verbosity > 5){
    pout() << "time_stepper::setup_poisson" << endl;
  }

  m_poisson = RefCountedPtr<poisson_solver> (new poisson_multifluid_gmg());
  m_poisson->set_verbosity(m_verbosity);
  m_poisson->set_amr(m_amr);
  m_poisson->set_computational_geometry(m_compgeom);
  m_poisson->set_physical_domain(m_physdom);
  m_poisson->set_potential(m_potential);

  // Boundary conditions, but these should definitely come through an input function. 
  if(SpaceDim == 2){
    m_poisson->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
    m_poisson->set_neumann_wall_bc(0,   Side::Hi, 0.0);
    m_poisson->set_dirichlet_wall_bc(1, Side::Lo, potential::ground);
    m_poisson->set_dirichlet_wall_bc(1, Side::Hi, potential::live);
  }
  else if(SpaceDim == 3){
    m_poisson->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
    m_poisson->set_neumann_wall_bc(0,   Side::Hi, 0.0);
    m_poisson->set_neumann_wall_bc(1,   Side::Lo, 0.0);                  
    m_poisson->set_neumann_wall_bc(1,   Side::Hi, 0.0);
    m_poisson->set_dirichlet_wall_bc(2, Side::Lo, potential::ground);
    m_poisson->set_dirichlet_wall_bc(2, Side::Hi, potential::live);
  }

  m_poisson->sanity_check();
  m_poisson->allocate_internals();
}

void time_stepper::setup_rte(){
  CH_TIME("time_stepper::setup_rte");
  if(m_verbosity > 5){
    pout() << "time_stepper::setup_rte" << endl;
  }

  m_rte = RefCountedPtr<rte_layout> (new rte_layout(m_plaskin));
  m_rte->set_verbosity(m_verbosity);
  m_rte->set_amr(m_amr);
  m_rte->set_computational_geometry(m_compgeom);
  m_rte->set_physical_domain(m_physdom);
  m_rte->sanity_check();
  m_rte->allocate_internals();
}

void time_stepper::setup_sigma(){
  CH_TIME("time_stepper::setup_poisson");
  if(m_verbosity > 5){
    pout() << "time_stepper::setup_poisson" << endl;
  }

  m_sigma = RefCountedPtr<sigma_solver> (new sigma_solver());
  m_sigma->set_amr(m_amr);
  m_sigma->set_verbosity(m_verbosity);
  m_sigma->set_computational_geometry(m_compgeom);
  m_sigma->set_plasma_kinetics(m_plaskin);
  m_sigma->set_physical_domain(m_physdom);
  m_sigma->allocate_internals();
}

void time_stepper::solve_poisson(){
  CH_TIME("time_stepper::solve_poisson()");
  if(m_verbosity > 5){
    pout() << "time_stepper::solve_poisson()" << endl;
  }

  this->compute_rho();
  m_poisson->solve(m_poisson->get_state(),
		   m_poisson->get_source(),
		   m_sigma->get_state(),
		   false);

#if 1 // should be removed
  m_poisson->write_plot_file();
#endif
}

void time_stepper::solve_poisson(MFAMRCellData&                a_potential,
				 MFAMRCellData*                a_rhs,
				 const Vector<EBAMRCellData*>  a_densities,
				 const EBAMRIVData&            a_sigma,
				 const centering::which_center a_centering){
  CH_TIME("time_stepper::solve_poisson(full)");
  if(m_verbosity > 5){
    pout() << "time_stepper::solve_poisson(full)" << endl;
  }

  const int ncomp = 1;
  
  MFAMRCellData* p_rhs;
  MFAMRCellData rhs;
  if(a_rhs == NULL){
    m_amr->allocate(rhs, ncomp);
    p_rhs = &rhs;
  }
  else {
    p_rhs = a_rhs;
  }

  this->compute_rho(rhs, a_densities, a_centering);

  m_poisson->solve(a_potential, *p_rhs, a_sigma, false);
}

void time_stepper::solve_rte(const Real a_dt){
  CH_TIME("time_stepper::solve_rte()");
  if(m_verbosity > 5){
    pout() << "time_stepper::solve_rte()" << endl;
  }

  const phase::which_phase rte_phase = m_rte->get_phase();

  EBAMRCellData E;
  m_amr->allocate(E, rte_phase, SpaceDim);
  this->compute_E(E, rte_phase, m_poisson->get_state());

  Vector<EBAMRCellData*> states     = m_rte->get_states();
  Vector<EBAMRCellData*> rhs        = m_rte->get_sources();
  Vector<EBAMRCellData*> cdr_states = m_cdr->get_states();

  this->solve_rte(states, rhs, cdr_states, E, a_dt, centering::cell_center);
}

void time_stepper::solve_rte(Vector<EBAMRCellData*>&       a_states,
			     Vector<EBAMRCellData*>&       a_rhs,
			     const Vector<EBAMRCellData*>& a_cdr_states,
			     const EBAMRCellData&          a_E,
			     const Real                    a_dt,
			     const centering::which_center a_centering){
  CH_TIME("time_stepper::solve_rte(full)");
  if(m_verbosity > 5){
    pout() << "time_stepper::solve_rte(full)" << endl;
  }

  this->compute_photon_source_terms(a_rhs, a_cdr_states, a_E, a_centering);


  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    
    RefCountedPtr<rte_solver>& solver = solver_it();
    EBAMRCellData& state              = *a_states[idx];
    EBAMRCellData& rhs                = *a_rhs[idx];
    solver->advance(a_dt, state, rhs);
  }

}

void time_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("time_stepper::synchronize_solver_times");
  if(m_verbosity > 5){
    pout() << "time_stepper::synchronize_solver_times" << endl;
  }

  m_cdr->set_time(a_step,     a_time, a_dt);
  m_poisson->set_time(a_step, a_time, a_dt);
  m_rte->set_time(a_step,     a_time, a_dt);
  m_sigma->set_time(a_step,   a_time, a_dt);
}

Real time_stepper::compute_relaxation_time(){
  CH_TIME("time_stepper::compute_relaxation_time");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_relaxation_time" << endl;
  }

  const int finest_level = 0;
  const Real tolerance   = 1.E-2;

  EBAMRCellData E, J, dt;
  m_amr->allocate(E,  m_cdr->get_phase(), SpaceDim);
  m_amr->allocate(J,  m_cdr->get_phase(), SpaceDim);
  m_amr->allocate(dt, m_cdr->get_phase(), SpaceDim);

  data_ops::set_value(dt, 1.E99);

  this->compute_E(E, m_cdr->get_phase(), m_poisson->get_state());
  this->compute_J(J);

  // Find the largest electric field in each direction
  Vector<Real> max_E(SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++){
    Real max, min;
    data_ops::get_max_min(max, min, E, dir);
    max_E[dir] = Max(Abs(max), Abs(min));
  }
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_cdr->get_phase())[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box         = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);
      
      EBCellFAB& dt_fab  = (*dt[lvl])[dit()];
      const EBCellFAB& e = (*E[lvl])[dit()];
      const EBCellFAB& j = (*J[lvl])[dit()];

      //
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();

	for (int dir = 0; dir < SpaceDim; dir++){
	  if(Abs(e(vof, dir)) > tolerance*max_E[dir]){
	    dt_fab(vof, dir) = Abs(units::s_eps0*e(vof, dir)/j(vof,dir));
	  }
	}
      }
    }
  }
  
  // Find the smallest dt
  Real min_dt = 1.E99;
  for (int dir = 0; dir < SpaceDim; dir++){
    Real max, min;
    data_ops::get_max_min(max, min, dt, dir);
    min_dt = Min(min_dt, min);
  }

  // Communicate the result
#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("time_stepper::compute_relaxation_time() - communication error on norm");
  }
  min_dt = tmp;
#endif

  pout() << "returning " << min_dt << endl;


  return min_dt;
}
