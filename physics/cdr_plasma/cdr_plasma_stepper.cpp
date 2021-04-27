/*!
  @file cdr_plasma_stepper.cpp
  @brief Implementation of cdr_plasma_stepper.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "cdr_plasma_stepper.H"
#include "field_solver_multigrid.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "units.H"

#include <ParmParse.H>
#include <EBLevelDataOps.H>

// C-style VLA methods can have occasional memory issues. Use them at your own peril. 
#define USE_FAST_REACTIONS  1
#define USE_FAST_VELOCITIES 0
#define USE_FAST_DIFFUSION  1

using namespace physics::cdr_plasma;

Real cdr_plasma_stepper::s_constant_one(const RealVect a_pos){
  return 1.0;
}

cdr_plasma_stepper::cdr_plasma_stepper(){
  m_class_name = "cdr_plasma_stepper";
  m_verbosity  = -1;
  m_solver_verbosity = -1;
  m_phase = phase::gas;
  m_realm = realm::primal;
  m_subcycle = false;

  set_poisson_wall_func(s_constant_one);
}

cdr_plasma_stepper::cdr_plasma_stepper(RefCountedPtr<cdr_plasma_physics>& a_physics) : cdr_plasma_stepper() {
  m_physics = a_physics;
}

void cdr_plasma_stepper::post_initialize() {
}

void cdr_plasma_stepper::register_realms(){
  m_amr->register_realm(m_realm);
}

void cdr_plasma_stepper::post_regrid(){

}

void cdr_plasma_stepper::register_operators(){
  CH_TIME("cdr_plasma_stepper::register_operators");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::register_operators" << endl;
  }
  
  m_cdr->register_operators();
  m_poisson->register_operators();
  m_rte->register_operators();
  m_sigma->register_operators();
}

cdr_plasma_stepper::~cdr_plasma_stepper(){
}

int cdr_plasma_stepper::query_ghost(){
  return 3;
}

bool cdr_plasma_stepper::stationary_rte(){
  CH_TIME("cdr_plasma_stepper::stationary_rte");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::stationary_rte" << endl;
  }

  return m_rte->is_stationary();
}

bool cdr_plasma_stepper::solve_poisson(){
  CH_TIME("cdr_plasma_stepper::solve_poisson");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::solve_poisson" << endl;
  }

  this->compute_rho();
  const bool converged = m_poisson->solve(m_poisson->get_state(),
					  m_poisson->get_source(),
					  m_sigma->get_state(),
					  false);

  return converged;
}

bool cdr_plasma_stepper::solve_poisson(MFAMRCellData&                a_potential,
				       MFAMRCellData&                a_rhs,
				       const Vector<EBAMRCellData*>  a_densities,
				       const EBAMRIVData&            a_sigma,
				       const centering               a_centering){
  CH_TIME("cdr_plasma_stepper::solve_poisson(full)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::solve_poisson(full)" << endl;
  }

  this->compute_rho(a_rhs, a_densities, a_centering);

  const bool converged = m_poisson->solve(a_potential, a_rhs, a_sigma, false);

  return converged;
}

void cdr_plasma_stepper::allocate_internals(){
  CH_TIME("cdr_plasma_stepper::allocate_internals"); 
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::allocate_internals" << endl;
  }

  /// Do nothing
}

void cdr_plasma_stepper::deallocate_internals(){
  CH_TIME("cdr_plasma_stepper::deallocate_internals"); 
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::deallocate_internals" << endl;
  }

  /// Do nothing
}

void cdr_plasma_stepper::advance_reaction_network(const Real a_time, const Real a_dt){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network(solvers)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network(solvers)" << endl;
  }

  // Compute the electric field
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, m_cdr->get_phase(), SpaceDim);
  this->compute_E(E, m_cdr->get_phase(), m_poisson->get_state());


  Vector<EBAMRCellData*> particle_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> photon_sources   = m_rte->get_sources();
  Vector<EBAMRCellData*> particle_states  = m_cdr->get_states();
  Vector<EBAMRCellData*> photon_states    = m_rte->get_states();

  // Call the AMR version (without the gradient)
  advance_reaction_network(particle_sources,
			   photon_sources, 
			   particle_states,
			   photon_states,
			   E,
			   a_time,
			   a_dt);
}

void cdr_plasma_stepper::advance_reaction_network(Vector<EBAMRCellData*>&       a_particle_sources,
						  Vector<EBAMRCellData*>&       a_photon_sources,
						  const Vector<EBAMRCellData*>& a_particle_densities,
						  const Vector<EBAMRCellData*>& a_photon_densities,
						  const EBAMRCellData&          a_E,
						  const Real&                   a_time,
						  const Real&                   a_dt){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network(nograd)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network(nograd)" << endl;
  }

  const int num_species = m_physics->get_num_cdr_species();

  Vector<EBAMRCellData*> grad_cdr(num_species); // Holders for grad(cdr)
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    grad_cdr[idx] = new EBAMRCellData();                            // This storage must be deleted
    m_amr->allocate(*grad_cdr[idx], m_realm, m_cdr->get_phase(), SpaceDim);  // Allocate

    m_amr->compute_gradient(*grad_cdr[idx], *a_particle_densities[idx], m_realm, phase::gas); // Compute grad()
    m_amr->average_down(*grad_cdr[idx], m_realm, m_cdr->get_phase());        // Average down
    m_amr->interp_ghost(*grad_cdr[idx], m_realm, m_cdr->get_phase());        // Interpolate ghost cells
  }

  this->advance_reaction_network(a_particle_sources,
				 a_photon_sources,
				 a_particle_densities,
				 grad_cdr,
				 a_photon_densities,
				 a_E,
				 a_time,
				 a_dt);

  // Delete extra storage - didn't use smart pointers for this...
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->deallocate(*grad_cdr[idx]);
    delete grad_cdr[idx];
  }
}

void cdr_plasma_stepper::advance_reaction_network(Vector<EBAMRCellData*>&       a_particle_sources,
						  Vector<EBAMRCellData*>&       a_photon_sources,
						  const Vector<EBAMRCellData*>& a_particle_densities,
						  const Vector<EBAMRCellData*>& a_particle_gradients,
						  const Vector<EBAMRCellData*>& a_photon_densities,
						  const EBAMRCellData&          a_E,
						  const Real&                   a_time,
						  const Real&                   a_dt){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network(amr)" << endl;
  }


  const int num_species = m_physics->get_num_cdr_species();
  const int num_photons = m_physics->get_num_rte_species();

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    Vector<LevelData<EBCellFAB>* > particle_sources(num_species);
    Vector<LevelData<EBCellFAB>* > particle_densities(num_species);
    Vector<LevelData<EBCellFAB>* > particle_gradients(num_species);
    Vector<LevelData<EBCellFAB>* > photon_sources(num_photons);
    Vector<LevelData<EBCellFAB>* > photon_densities(num_photons);

    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      particle_sources[idx]   = (*a_particle_sources[idx])[lvl];
      particle_densities[idx] = (*a_particle_densities[idx])[lvl];
      particle_gradients[idx] = (*a_particle_gradients[idx])[lvl];
    }

    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      photon_sources[idx]   = (*a_photon_sources[idx])[lvl];
      photon_densities[idx] = (*a_photon_densities[idx])[lvl];
    }

    // Call the level versions
    advance_reaction_network(particle_sources,
			     photon_sources,
			     particle_densities,
			     particle_gradients,
			     photon_densities,
			     *a_E[lvl],
			     a_time,
			     a_dt,
			     lvl);
  }

#if 0 // R.M. May/2020: This is not a good place to do this kind of averaging and interpolating. If need it, do it elsewhere.
  // Average down species
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->average_down(*a_particle_sources[idx], m_realm, m_cdr->get_phase());
    m_amr->interp_ghost(*a_particle_sources[idx], m_realm, m_cdr->get_phase()); // This MAY be when if we extrapolate advection
  }

  // Average down photon solvers
  for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->average_down(*a_photon_sources[idx], m_realm, m_cdr->get_phase());
  }
#endif
}

void cdr_plasma_stepper::advance_reaction_network(Vector<LevelData<EBCellFAB>* >&       a_particle_sources,
						  Vector<LevelData<EBCellFAB>* >&       a_photon_sources,
						  const Vector<LevelData<EBCellFAB>* >& a_particle_densities,
						  const Vector<LevelData<EBCellFAB>* >& a_particle_gradients,
						  const Vector<LevelData<EBCellFAB> *>& a_photon_densities,
						  const LevelData<EBCellFAB>&           a_E,
						  const Real&                           a_time,
						  const Real&                           a_dt,
						  const int                             a_lvl){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network(level)" << endl;
  }

  const Real zero = 0.0;

  const int num_photons  = m_physics->get_num_rte_species();
  const int num_species  = m_physics->get_num_cdr_species();


  // Stencils for extrapolating things to cell centroids
  const irreg_amr_stencil<centroid_interp>& interp_stencils = m_amr->get_centroid_interp_stencils(m_realm, m_cdr->get_phase());

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[a_lvl];
  const Real dx                = m_amr->get_dx()[a_lvl];
  const RealVect origin        = m_amr->get_prob_lo();
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    Vector<EBCellFAB*> particle_sources(num_species);
    Vector<EBCellFAB*> particle_densities(num_species);
    Vector<EBCellFAB*> particle_gradients(num_species);
    Vector<EBCellFAB*> particle_velocities(num_species, NULL);
    Vector<EBCellFAB*> photon_sources(num_photons);
    Vector<EBCellFAB*> photon_densities(num_photons);

    
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<cdr_solver>& solver = solver_it();
      const int idx = solver_it.index();
      particle_sources[idx]   = &(*a_particle_sources[idx])[dit()];
      particle_densities[idx] = &(*a_particle_densities[idx])[dit()];
      particle_gradients[idx] = &(*a_particle_gradients[idx])[dit()];

      if(solver->is_mobile()){
	const EBAMRCellData& velo = solver->get_velo_cell();
	particle_velocities[idx] = &((*velo[a_lvl])[dit()]);
      }
    }

    
    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      photon_sources[idx]   = &(*a_photon_sources[idx])[dit()];
      photon_densities[idx] = &(*a_photon_densities[idx])[dit()];
    }
    
    // This does all cells
#if USE_FAST_REACTIONS
    ParmParse pp ("cache_blocking");

    int block_loops = 0;
    pp.query("tile_loops", block_loops);
    if(block_loops == 0){
      advance_reaction_network_reg_fast(particle_sources,
					photon_sources,
					particle_densities,
					particle_gradients,
					photon_densities,
					a_E[dit()],
					a_time,
					a_dt,
					dx,
					dbl.get(dit()));
    }
    else{
      Vector<int> tilesize(SpaceDim);
      pp.getarr("tile_size", tilesize, 0, SpaceDim);
      Vector<Box> boxes = m_amr->make_tiles(dbl.get(dit()), IntVect(D_DECL(tilesize[0], tilesize[1], tilesize[2])));
      for (int ibox = 0; ibox < boxes.size(); ibox++){
	advance_reaction_network_reg_fast(particle_sources,
					  photon_sources,
					  particle_densities,
					  particle_gradients,
					  photon_densities,
					  a_E[dit()],
					  a_time,
					  a_dt,
					  dx,
					  boxes[ibox]);
      }
    }
#else
    advance_reaction_network_reg(particle_sources,
				 photon_sources,
				 particle_densities,
				 particle_gradients,
				 photon_densities,
				 a_E[dit()],
				 a_time,
				 a_dt,
				 dx,
				 dbl.get(dit()));
#endif

    // This overwrites irregular celles
    advance_reaction_network_irreg(particle_sources,
				   photon_sources,
				   particle_densities,
				   particle_gradients,
				   particle_velocities,
				   photon_densities,
				   interp_stencils[a_lvl][dit()],
				   a_E[dit()],
				   a_time,
				   a_dt,
				   dx,
				   dbl.get(dit()),
				   a_lvl,
				   dit());
  }
}

void cdr_plasma_stepper::advance_reaction_network_reg(Vector<EBCellFAB*>&       a_particle_sources,
						      Vector<EBCellFAB*>&       a_photon_sources,
						      const Vector<EBCellFAB*>& a_particle_densities,
						      const Vector<EBCellFAB*>& a_particle_gradients,
						      const Vector<EBCellFAB*>& a_photon_densities,
						      const EBCellFAB&          a_E,
						      const Real&               a_time,
						      const Real&               a_dt,
						      const Real&               a_dx,
						      const Box&                a_box){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network_reg(patch)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network_reg(patch)" << endl;
  }

  const Real zero = 0.0;
    

  const int num_species  = m_physics->get_num_cdr_species();
  const int num_photons  = m_physics->get_num_rte_species();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     particle_densities(num_species);
  Vector<RealVect> particle_gradients(num_species);
  Vector<Real>     photon_sources(num_photons);
  Vector<Real>     photon_densities(num_photons);

  // Computed source terms onto here
  EBCellFAB part_src(ebisbox, a_E.getRegion(), num_species);
  EBCellFAB phot_src(ebisbox, a_E.getRegion(), num_photons);

  part_src.setVal(0.0);
  phot_src.setVal(0.0);

  const BaseFab<Real>& EFab     = a_E.getSingleValuedFAB();

  
  // Grid loop
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    if(ebisbox.isRegular(iv)){
    
      // Position and E
      pos    = origin + RealVect(iv)*a_dx;
      E      = RealVect(D_DECL(EFab(iv, 0),     EFab(iv, 1),     EFab(iv, 2)));

      // Fill vectors with densities
      for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_particle_densities[idx]).getSingleValuedFAB()(iv, 0);
	particle_densities[idx] = Max(zero, phi);
	particle_gradients[idx] = RealVect(D_DECL((*a_particle_gradients[idx]).getSingleValuedFAB()(iv, 0),
						  (*a_particle_gradients[idx]).getSingleValuedFAB()(iv, 1),
						  (*a_particle_gradients[idx]).getSingleValuedFAB()(iv, 2)));
      }
      for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_photon_densities[idx]).getSingleValuedFAB()(iv, 0);
	photon_densities[idx] = Max(zero, phi);
      }

      // Compute source terms
      const Real kappa = 1.0; // Kappa for regular cells
      m_physics->advance_reaction_network(particle_sources,
					  photon_sources,
					  particle_densities,
					  particle_gradients,
					  photon_densities,
					  E,
					  pos,
					  a_dx,
					  a_dt,
					  a_time,
					  kappa);

      // Put vector into temporary holders
      for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	part_src.getSingleValuedFAB()(iv,idx) = particle_sources[idx];
      }
    
      // Put vector into temporary holders
      for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	phot_src.getSingleValuedFAB()(iv,idx) = photon_sources[idx];
      }
    }
  }

  // Copy temporary storage back to solvers
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    (*a_particle_sources[idx]).setVal(0.0);
    (*a_particle_sources[idx]).plus(part_src, idx, 0, 1);

    // Covered cells are bogus. 
    (*a_particle_sources[idx]).setCoveredCellVal(0.0, 0);
  }

  // Copy temporary storage back to solvers
  for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    (*a_photon_sources[idx]).setVal(0.0);
    (*a_photon_sources[idx]).plus(phot_src, idx, 0, 1);

    // Covered cells are bogus. 
    (*a_photon_sources[idx]).setCoveredCellVal(0.0, 0);
  }
}

void cdr_plasma_stepper::advance_reaction_network_reg_fast(Vector<EBCellFAB*>&       a_particle_sources,
							   Vector<EBCellFAB*>&       a_photon_sources,
							   const Vector<EBCellFAB*>& a_particle_densities,
							   const Vector<EBCellFAB*>& a_particle_gradients,
							   const Vector<EBCellFAB*>& a_photon_densities,
							   const EBCellFAB&          a_E,
							   const Real&               a_time,
							   const Real&               a_dt,
							   const Real&               a_dx,
							   const Box&                a_box){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network_reg_fast(patch)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network_reg_fast(patch)" << endl;
  }

#if CH_SPACEDIM==2
  advance_reaction_network_reg_fast2D(a_particle_sources, a_photon_sources, a_particle_densities, a_particle_gradients,
				      a_photon_densities, a_E, a_time, a_dt, a_dx, a_box);
#elif CH_SPACEDIM==3
  advance_reaction_network_reg_fast3D(a_particle_sources, a_photon_sources, a_particle_densities, a_particle_gradients,
				      a_photon_densities, a_E, a_time, a_dt, a_dx, a_box);
#endif
  
}

void cdr_plasma_stepper::advance_reaction_network_reg_fast2D(Vector<EBCellFAB*>&       a_particle_sources,
							     Vector<EBCellFAB*>&       a_photon_sources,
							     const Vector<EBCellFAB*>& a_particle_densities,
							     const Vector<EBCellFAB*>& a_particle_gradients,
							     const Vector<EBCellFAB*>& a_photon_densities,
							     const EBCellFAB&          a_E,
							     const Real&               a_time,
							     const Real&               a_dt,
							     const Real&               a_dx,
							     const Box&                a_box){
#if CH_SPACEDIM==2
  CH_TIME("cdr_plasma_stepper::advance_reaction_network_reg_fast2D(patch)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network_reg_fast2D(patch)" << endl;
  }

  const Real zero = 0.0;
    

  const int num_species  = m_physics->get_num_cdr_species();
  const int num_photons  = m_physics->get_num_rte_species();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics. 
  RealVect         pos = RealVect::Zero;
  RealVect         E   = RealVect::Zero;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     particle_densities(num_species);
  Vector<RealVect> particle_gradients(num_species);
  Vector<Real>     photon_sources(num_photons);
  Vector<Real>     photon_densities(num_photons);

  // I need contiguous memory for what is about to happen, so begin by copying stuff onto smaller data holders
  
  // Temps for source terms and particle densities
  FArrayBox cdr_src(a_box, num_species);
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox cdr_gx(a_box, num_species);  // Gradient in x-direction
  FArrayBox cdr_gy(a_box, num_species);  // Gradient in y-direction
  cdr_phi.setVal(0.0);
  cdr_src.setVal(0.0);
  for (int i = 0; i < num_species; i++){
    cdr_phi.copy(a_particle_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
    cdr_gx.copy(a_particle_gradients[i]->getFArrayBox(),  a_box, 0, a_box, i, 1);
    cdr_gy.copy(a_particle_gradients[i]->getFArrayBox(),  a_box, 1, a_box, i, 1);
  }

  // Temps for photon source terms and densities
  FArrayBox rte_phi(a_box, num_photons);
  FArrayBox rte_src(a_box, num_photons);
  rte_phi.setVal(0.0);
  rte_src.setVal(0.0);
  for (int i = 0; i < num_photons; i++){
    rte_phi.copy(a_photon_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
  }

  // Temp for electric field
  FArrayBox Efab(a_box, SpaceDim);
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);

  // Pointer offsets
  const IntVect dims   = a_box.size();
  const IntVect lo     = a_box.smallEnd();
  const IntVect hi     = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int offset     = lo[0] + n0*lo[1];

  // C style variable-length array conversion magic
  auto vla_cdr_src = (Real (*__restrict__)[n1][n0]) (cdr_src.dataPtr() - offset);
  auto vla_rte_src = (Real (*__restrict__)[n1][n0]) (rte_src.dataPtr() - offset);
  auto vla_cdr_phi = (Real (*__restrict__)[n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_rte_phi = (Real (*__restrict__)[n1][n0]) (rte_phi.dataPtr() - offset);
  auto vla_E       = (Real (*__restrict__)[n1][n0]) (Efab.dataPtr()    - offset);
  auto vla_cdr_gx  = (Real (*__restrict__)[n1][n0]) (cdr_gx.dataPtr()  - offset);
  auto vla_cdr_gy  = (Real (*__restrict__)[n1][n0]) (cdr_gy.dataPtr()  - offset);

  for (int j = lo[1]; j <= hi[1]; ++j){
    for (int i = lo[0]; i <= hi[0]; ++i){
      
      // Particle densities
      for (int idx = 0; idx < num_species; ++idx){
	particle_densities[idx] = Max(0.0, vla_cdr_phi[idx][j][i]);
	particle_gradients[idx][0] = vla_cdr_gx[idx][j][i];
	particle_gradients[idx][1] = vla_cdr_gy[idx][j][i];
      }

      // Photon densities
      for (int idx = 0; idx < num_photons; ++idx){
	photon_densities[idx] = Max(0.0, vla_rte_phi[idx][j][i]);
      }

      E   = RealVect(vla_E[0][j][i], vla_E[1][j][i]);
      pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

      m_physics->advance_reaction_network(particle_sources,
					  photon_sources,
					  particle_densities,
					  particle_gradients,
					  photon_densities,
					  E,
					  pos,
					  a_dx,
					  a_dt,
					  a_time,
					  1.0);

      // Put result in correct palce
      for (int idx = 0; idx < num_species; ++idx){
	vla_cdr_src[idx][j][i] = particle_sources[idx];
      }

      // Put result in correct palce
      for (int idx = 0; idx < num_photons; ++idx){
	vla_rte_src[idx][j][i] = photon_sources[idx];
      }
    }
  }

  // Copy result back to solvers
  for (int i = 0; i < num_species; i++){
    FArrayBox& src = a_particle_sources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(cdr_src, a_box, i, a_box, 0, 1);
    a_particle_sources[i]->setCoveredCellVal(0.0, 0);
  }

  // Copy result back to solvers
  for (int i = 0; i < num_photons; i++){
    FArrayBox& src = a_photon_sources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(rte_src, a_box, i, a_box, 0, 1);
  }
#endif
}


void cdr_plasma_stepper::advance_reaction_network_reg_fast3D(Vector<EBCellFAB*>&       a_particle_sources,
							     Vector<EBCellFAB*>&       a_photon_sources,
							     const Vector<EBCellFAB*>& a_particle_densities,
							     const Vector<EBCellFAB*>& a_particle_gradients,
							     const Vector<EBCellFAB*>& a_photon_densities,
							     const EBCellFAB&          a_E,
							     const Real&               a_time,
							     const Real&               a_dt,
							     const Real&               a_dx,
							     const Box&                a_box){
#if CH_SPACEDIM==3
  CH_TIME("cdr_plasma_stepper::advance_reaction_network_reg_fast3D(patch)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network_reg_fast3D(patch)" << endl;
  }

  const Real zero = 0.0;

  const int num_species  = m_physics->get_num_cdr_species();
  const int num_photons  = m_physics->get_num_rte_species();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics. 
  RealVect         pos = RealVect::Zero;
  RealVect         E   = RealVect::Zero;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     particle_densities(num_species);
  Vector<RealVect> particle_gradients(num_species);
  Vector<Real>     photon_sources(num_photons);
  Vector<Real>     photon_densities(num_photons);

  // I need contiguous memory for what is about to happen, so begin by copying stuff onto smaller data holders
  
  // Temps for source terms and particle densities
  FArrayBox cdr_src(a_box, num_species);
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox cdr_gx(a_box, num_species);  // Gradient in x-direction
  FArrayBox cdr_gy(a_box, num_species);  // Gradient in y-direction
  FArrayBox cdr_gz(a_box, num_species);  // Gradient in z-direction
  cdr_phi.setVal(0.0);
  cdr_src.setVal(0.0);
  for (int i = 0; i < num_species; i++){
    cdr_phi.copy(a_particle_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
    cdr_gx.copy(a_particle_gradients[i]->getFArrayBox(),  a_box, 0, a_box, i, 1);
    cdr_gy.copy(a_particle_gradients[i]->getFArrayBox(),  a_box, 1, a_box, i, 1);
    cdr_gz.copy(a_particle_gradients[i]->getFArrayBox(),  a_box, 2, a_box, i, 1);
  }

  // Temps for photon source terms and densities
  FArrayBox rte_phi(a_box, num_photons);
  FArrayBox rte_src(a_box, num_photons);
  rte_phi.setVal(0.0);
  rte_src.setVal(0.0);
  for (int i = 0; i < num_photons; i++){
    rte_phi.copy(a_photon_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
  }

  // Temp for electric field
  FArrayBox Efab(a_box, SpaceDim);
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);

  // Pointer offsets
  const IntVect dims   = a_box.size();
  const IntVect lo     = a_box.smallEnd();
  const IntVect hi     = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int n2         = dims[2];
  const int offset     = lo[0] + n0*(lo[1] + n1*lo[2]);

  // C style variable-length array conversion magic
  auto vla_cdr_src = (Real (*__restrict__)[n2][n1][n0]) (cdr_src.dataPtr() - offset);
  auto vla_rte_src = (Real (*__restrict__)[n2][n1][n0]) (rte_src.dataPtr() - offset);
  auto vla_cdr_phi = (Real (*__restrict__)[n2][n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_rte_phi = (Real (*__restrict__)[n2][n1][n0]) (rte_phi.dataPtr() - offset);
  auto vla_E       = (Real (*__restrict__)[n2][n1][n0]) (Efab.dataPtr()    - offset);
  auto vla_cdr_gx  = (Real (*__restrict__)[n2][n1][n0]) (cdr_gx.dataPtr()  - offset);
  auto vla_cdr_gy  = (Real (*__restrict__)[n2][n1][n0]) (cdr_gy.dataPtr()  - offset);
  auto vla_cdr_gz  = (Real (*__restrict__)[n2][n1][n0]) (cdr_gz.dataPtr()  - offset);

  for (int k = lo[2]; k <= hi[2]; k++){
    for (int j = lo[1]; j <= hi[1]; ++j){
      for (int i = lo[0]; i <= hi[0]; ++i){
      
	// Particle densities
	for (int idx = 0; idx < num_species; ++idx){
	  particle_densities[idx]    = Max(0.0, vla_cdr_phi[idx][k][j][i]);
	  particle_gradients[idx][0] = vla_cdr_gx[idx][k][j][i];
	  particle_gradients[idx][1] = vla_cdr_gy[idx][k][j][i];
	  particle_gradients[idx][2] = vla_cdr_gz[idx][k][j][i];
	}

	// Photon densities
	for (int idx = 0; idx < num_photons; ++idx){
	  photon_densities[idx] = Max(0.0, vla_rte_phi[idx][k][j][i]);
	}

	E   = RealVect(vla_E[0][k][j][i], vla_E[1][k][j][i], vla_E[2][k][j][i]);
	pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

	m_physics->advance_reaction_network(particle_sources,
					    photon_sources,
					    particle_densities,
					    particle_gradients,
					    photon_densities,
					    E,
					    pos,
					    a_dx,
					    a_dt,
					    a_time,
					    1.0);

	// Put result in correct palce
	for (int idx = 0; idx < num_species; ++idx){
	  vla_cdr_src[idx][k][j][i] = particle_sources[idx];
	}

	// Put result in correct palce
	for (int idx = 0; idx < num_photons; ++idx){
	  vla_rte_src[idx][k][j][i] = photon_sources[idx];
	}
      }
    }
  }

  // Copy result back to solvers
  for (int i = 0; i < num_species; i++){
    FArrayBox& src = a_particle_sources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(cdr_src, a_box, i, a_box, 0, 1);
    a_particle_sources[i]->setCoveredCellVal(0.0, 0);
  }

  // Copy result back to solvers
  for (int i = 0; i < num_photons; i++){
    FArrayBox& src = a_photon_sources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(rte_src, a_box, i, a_box, 0, 1);
  }

#endif
}



void cdr_plasma_stepper::advance_reaction_network_irreg(Vector<EBCellFAB*>&          a_particle_sources,
							Vector<EBCellFAB*>&          a_photon_sources,
							const Vector<EBCellFAB*>&    a_particle_densities,
							const Vector<EBCellFAB*>&    a_particle_gradients,
							const Vector<EBCellFAB*>&    a_particle_velocities,
							const Vector<EBCellFAB*>&    a_photon_densities,
							const BaseIVFAB<VoFStencil>& a_interp_stencils,
							const EBCellFAB&             a_E,
							const Real&                  a_time,
							const Real&                  a_dt,
							const Real&                  a_dx,
							const Box&                   a_box,
							const int                    a_lvl,
							const DataIndex&             a_dit){
  if(m_interp_sources){
    advance_reaction_network_irreg_interp(a_particle_sources,
					  a_photon_sources,
					  a_particle_densities,
					  a_particle_gradients,
					  a_particle_velocities,
					  a_photon_densities,
					  a_interp_stencils,
					  a_E,
					  a_time,
					  a_dt,
					  a_dx,
					  a_box,
					  a_lvl,
					  a_dit);
  }
  else{
    advance_reaction_network_irreg_kappa(a_particle_sources,
					 a_photon_sources,
					 a_particle_densities,
					 a_particle_gradients,
					 a_photon_densities,
					 a_interp_stencils,
					 a_E,
					 a_time,
					 a_dt,
					 a_dx,
					 a_box,
					 a_lvl,
					 a_dit);
  }
}

void cdr_plasma_stepper::advance_reaction_network_irreg_interp(Vector<EBCellFAB*>&          a_particle_sources,
							       Vector<EBCellFAB*>&          a_photon_sources,
							       const Vector<EBCellFAB*>&    a_particle_densities,
							       const Vector<EBCellFAB*>&    a_particle_gradients,
							       const Vector<EBCellFAB*>&    a_particle_velocities,
							       const Vector<EBCellFAB*>&    a_photon_densities,
							       const BaseIVFAB<VoFStencil>& a_interp_stencils,
							       const EBCellFAB&             a_E,
							       const Real&                  a_time,
							       const Real&                  a_dt,
							       const Real&                  a_dx,
							       const Box&                   a_box,
							       const int                    a_lvl,
							       const DataIndex&             a_dit){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network_irreg_interp(patch)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network_irreg_interp(patch)" << endl;
  }

  const Real zero = 0.0;
  
  const int num_photons  = m_physics->get_num_rte_species();
  const int num_species  = m_physics->get_num_cdr_species();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     particle_sources(num_species, 0.);
  Vector<Real>     photon_sources(num_photons, 0.);
  Vector<Real>     particle_densities(num_species, 0.);
  Vector<RealVect> particle_gradients(num_species, RealVect::Zero);
  Vector<Real>     photon_densities(num_photons, 0.);

  VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof       = vofit();
    const Real kappa          = ebisbox.volFrac(vof);
    const VoFStencil& stencil = a_interp_stencils(vof, 0);

    pos = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);

    // Compute electric field on centroids
    E      = RealVect::Zero;
    for (int i = 0; i < stencil.size(); i++){
      const VolIndex& ivof = stencil.vof(i);
      const Real& iweight  = stencil.weight(i);
      for (int dir = 0; dir < SpaceDim; dir++){
	E[dir] += a_E(ivof, dir)*iweight;
      }
    }

    // Compute cdr_densities and their gradients on centroids
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<cdr_solver>& solver = solver_it();
      const int idx = solver_it.index();

      bool applyStencil = true;

#if 0
      if(solver->is_mobile()){
	const EBCellFAB& velo = *a_particle_velocities[idx];
	const RealVect v = RealVect(D_DECL(velo(vof,0), velo(vof,1), velo(vof,2)));
	if(PolyGeom::dot(v,ebisbox.normal(vof)) > 0.0){ // Flow away from the boundary.
	  applyStencil = false;
	}
	else{
	  applyStencil = true;
	}
      }
      else{
	applyStencil = true;
      }
#endif

      if(applyStencil){
	Real phi = 0.0;
	RealVect grad = RealVect::Zero;
	for (int i = 0; i < stencil.size(); i++){
	  const VolIndex& ivof = stencil.vof(i);
	  const Real& iweight  = stencil.weight(i);

	  phi += (*a_particle_densities[idx])(ivof,0)*iweight;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    grad[dir] += (*a_particle_gradients[idx])(ivof, dir)*iweight;
	  }
	}
      
	particle_densities[idx] = Max(zero, phi);
	particle_gradients[idx] = grad;
      }
      else{
	particle_densities[idx] = Max(zero, (*a_particle_densities[idx])(vof,0));
	particle_gradients[idx] = RealVect(D_DECL((*a_particle_gradients[idx])(vof,0),
						  (*a_particle_gradients[idx])(vof,1),
						  (*a_particle_gradients[idx])(vof,2)));
      }
    }

    // Compute RTE densities on the centroids
    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();

      Real phi = 0.0;
      for (int i = 0; i < stencil.size(); i++){
	const VolIndex& ivof = stencil.vof(i);
	const Real& iweight  = stencil.weight(i);
	phi += (*a_photon_densities[idx])(ivof, 0)*iweight;
      }
      photon_densities[idx] = Max(zero, phi);
    }

    // Compute source terms
    m_physics->advance_reaction_network(particle_sources,
					photon_sources,
					particle_densities,
					particle_gradients,
					photon_densities,
					E,
					pos,
					a_dx,
					a_dt,
					a_time,
					kappa);

    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_particle_sources[idx])(vof, 0) = particle_sources[idx];
    }
    
    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_photon_sources[idx])(vof, 0) = photon_sources[idx];
    }
  }
}

void cdr_plasma_stepper::advance_reaction_network_irreg_kappa(Vector<EBCellFAB*>&          a_particle_sources,
							      Vector<EBCellFAB*>&          a_photon_sources,
							      const Vector<EBCellFAB*>&    a_particle_densities,
							      const Vector<EBCellFAB*>&    a_particle_gradients,
							      const Vector<EBCellFAB*>&    a_photon_densities,
							      const BaseIVFAB<VoFStencil>& a_interp_stencils,
							      const EBCellFAB&             a_E,
							      const Real&                  a_time,
							      const Real&                  a_dt,
							      const Real&                  a_dx,
							      const Box&                   a_box,
							      const int                    a_lvl,
							      const DataIndex&             a_dit){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network_irreg_kappa(patch)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network_irreg_kappa(patch)" << endl;
  }

  const Real zero = 0.0;
  
  const int num_photons  = m_physics->get_num_rte_species();
  const int num_species  = m_physics->get_num_cdr_species();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     photon_sources(num_photons);
  Vector<Real>     particle_densities(num_species);
  Vector<RealVect> particle_gradients(num_species);
  Vector<Real>     photon_densities(num_photons);

  VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof       = vofit();
    const Real kappa          = ebisbox.volFrac(vof);
    const VoFStencil& stencil = a_interp_stencils(vof, 0);

    pos = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);

    // Compute electric field on centroids
    E = RealVect::Zero;
    for (int i = 0; i < stencil.size(); i++){
      const VolIndex& ivof = stencil.vof(i);
      const Real& iweight  = stencil.weight(i);
      for (int dir = 0; dir < SpaceDim; dir++){
	E[dir] += a_E(ivof, dir)*iweight;
      }
    }

    // Compute cdr_densities and their gradients on centroids
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();

      // Gradient on centroids for some reason
      RealVect grad = RealVect::Zero;
      for (int i = 0; i < stencil.size(); i++){
	const VolIndex& ivof = stencil.vof(i);
	const Real& iweight  = stencil.weight(i);
	for (int dir = 0; dir < SpaceDim; dir++){
	  grad[dir] += (*a_particle_gradients[idx])(ivof, dir)*iweight;
	}
      }
      
      particle_densities[idx] = Max(zero, (*a_particle_densities[idx])(vof,0));
      particle_gradients[idx] = grad;
    }

    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      photon_densities[idx] = Max(zero, (*a_photon_densities[idx])(vof, 0));
    }

    // Compute source terms
    m_physics->advance_reaction_network(particle_sources,
					photon_sources,
					particle_densities,
					particle_gradients,
					photon_densities,
					E,
					pos,
					a_dx,
					a_dt,
					a_time,
					kappa);

    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_particle_sources[idx])(vof, 0) = particle_sources[idx];
    }
    
    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_photon_sources[idx])(vof, 0) = photon_sources[idx];
    }
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_face(Vector<EBAMRFluxData*>&       a_diffco_face,
						 const Vector<EBAMRCellData*>& a_cdr_densities,
						 const EBAMRCellData&          a_E,
						 const Real&                   a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_face(full)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_face(full)" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->get_num_cdr_species();
  const int finest_level = m_amr->get_finest_level();

  // Allocate data for cell-centered diffusion coefficients
  Vector<EBAMRCellData> diffco(num_species);
  Vector<Real> cdr_densities(num_species);
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();

    m_amr->allocate(diffco[idx], m_realm, m_cdr->get_phase(), ncomp);
  }

  // Call the cell version
  compute_cdr_diffco_cell(diffco, a_cdr_densities, a_E, a_time);

  // Now compute face-centered things by taking the average of cell-centered things
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const int idx = solver_it.index();

    if(solver->is_diffusive()){ // Only need to do this for diffusive things
      data_ops::set_value(*a_diffco_face[idx], 1.0);
      m_amr->average_down(diffco[idx], m_realm, m_cdr->get_phase());
      m_amr->interp_ghost(diffco[idx], m_realm, m_cdr->get_phase()); 

      data_ops::average_cell_to_face_allcomps(*a_diffco_face[idx], diffco[idx], m_amr->get_domains());
    }
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_cell(Vector<EBAMRCellData>&       a_diffco_cell,
						 const Vector<EBAMRCellData*>& a_cdr_densities,
						 const EBAMRCellData&          a_E,
						 const Real&                   a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_cell(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_cell(amr)" << endl;
  }
  
  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->get_num_cdr_species();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<LevelData<EBCellFAB>* > diffco(num_species);
    Vector<LevelData<EBCellFAB>* > cdr_densities(num_species);
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx] = a_diffco_cell[idx][lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    // Cell the level version
    compute_cdr_diffco_cell(diffco, cdr_densities, *a_E[lvl], lvl, a_time);
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_cell(Vector<LevelData<EBCellFAB>* >&       a_diffco_cell,
						 const Vector<LevelData<EBCellFAB>* >& a_cdr_densities,
						 const LevelData<EBCellFAB>&           a_E,
						 const int                             a_lvl,
						 const Real&                           a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_cell(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_cell(level)" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->get_num_cdr_species();

  const irreg_amr_stencil<centroid_interp>& interp_stencils = m_amr->get_centroid_interp_stencils(m_realm, m_cdr->get_phase());

  // Call the level version
  const DisjointBoxLayout& dbl  = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[a_lvl];
  const Real dx                 = m_amr->get_dx()[a_lvl];
  const RealVect origin         = m_amr->get_prob_lo();
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = IntVectSet(box);
    const EBCellFAB& E     = a_E[dit()];

    Vector<EBCellFAB*> diffco(num_species);
    Vector<EBCellFAB*> cdr_densities(num_species);

    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx]        = &(*a_diffco_cell[idx])[dit()];
      cdr_densities[idx] = &(*a_cdr_densities[idx])[dit()];
    }

    // Regular cells
#if USE_FAST_DIFFUSION
    compute_cdr_diffco_cell_reg_fast(diffco, cdr_densities, E, dbl.get(dit()), m_amr->get_dx()[a_lvl], a_time);
#else
    compute_cdr_diffco_cell_reg(diffco, cdr_densities, E, dbl.get(dit()), m_amr->get_dx()[a_lvl], a_time);
#endif

    // Irregular cells
    compute_cdr_diffco_cell_irreg(diffco, cdr_densities, E, dbl.get(dit()), m_amr->get_dx()[a_lvl],
				  interp_stencils[a_lvl][dit()], a_time, a_lvl, dit());


  }
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_reg(Vector<EBCellFAB*>&       a_diffco_cell,
						     const Vector<EBCellFAB*>& a_cdr_densities,
						     const EBCellFAB&          a_E,
						     const Box                 a_box,
						     const Real                a_dx,
						     const Real                a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_cell_reg");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_cell_reg" << endl;
  }

  const int comp        = 0;
  const int num_species = m_physics->get_num_cdr_species();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     cdr_densities(num_species);

  // Computed coefficients terms onto here
  EBCellFAB tmp(ebisbox, a_E.getRegion(), num_species);

  const BaseFab<Real>& EFab  = a_E.getSingleValuedFAB();
  
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    
    pos = origin + iv*a_dx;
    E   = RealVect(D_DECL(EFab(iv, 0),     EFab(iv, 1),     EFab(iv, 2)));

    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      const Real phi = (*a_cdr_densities[idx]).getSingleValuedFAB()(iv, comp);
      cdr_densities[idx] = Max(0.0, phi);
    }

    const Vector<Real> coeffs = m_physics->compute_cdr_diffusion_coefficients(a_time,
    									      pos,
    									      E,
    									      cdr_densities);

    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->is_diffusive()){
    	tmp.getSingleValuedFAB()(iv,idx) = coeffs[idx];
      }
    }
  }

  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    if(solver_it()->is_diffusive()){
      (*a_diffco_cell[idx]).setVal(0.0);
      (*a_diffco_cell[idx]).plus(tmp, idx, 0, 1);

      // Covered cells are bogus. 
      (*a_diffco_cell[idx]).setCoveredCellVal(0.0, 0);
    }
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast(Vector<EBCellFAB*>&       a_diffco_cell,
							  const Vector<EBCellFAB*>& a_cdr_densities,
							  const EBCellFAB&          a_E,
							  const Box                 a_box,
							  const Real                a_dx,
							  const Real                a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast" << endl;
  }

#if CH_SPACEDIM==2
  compute_cdr_diffco_cell_reg_fast2D(a_diffco_cell, a_cdr_densities, a_E, a_box, a_dx, a_time);
#elif CH_SPACEDIM==3
  compute_cdr_diffco_cell_reg_fast3D(a_diffco_cell, a_cdr_densities, a_E, a_box, a_dx, a_time);
#endif
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast2D(Vector<EBCellFAB*>&       a_diffco_cell,
							    const Vector<EBCellFAB*>& a_cdr_densities,
							    const EBCellFAB&          a_E,
							    const Box                 a_box,
							    const Real                a_dx,
							    const Real                a_time){
#if CH_SPACEDIM==2
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast2D");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast2D" << endl;
  }

  const int comp        = 0;
  const int num_species = m_physics->get_num_cdr_species();
  const RealVect origin = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics
  RealVect     pos, E;
  Vector<Real> cdr_densities(num_species);

  // I need contiguous memory for the nasty that is about to happen. So begin by copying things onto smaller data holders
  FArrayBox cdr_dco(a_box, num_species);
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox Efab(a_box, SpaceDim);
  for (int i = 0; i < num_species; i++){
    cdr_phi.copy(a_cdr_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
  }
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);

  // Pointer offsets
  const IntVect dims = a_box.size();
  const IntVect lo   = a_box.smallEnd();
  const IntVect hi   = a_box.bigEnd();
  const int n0       = dims[0];
  const int n1       = dims[1];
  const int offset   = lo[0] + n0*lo[1];

  // C style variable-length array conversion magic
  auto vla_cdr_dco = (Real (*__restrict__)[n1][n0]) (cdr_dco.dataPtr() - offset);
  auto vla_cdr_phi = (Real (*__restrict__)[n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_E       = (Real (*__restrict__)[n1][n0]) (Efab.dataPtr()    - offset);

  for (int j = lo[1]; j <= hi[1]; ++j){
    for (int i = lo[0]; i <= hi[0]; ++i){
      
      // Particle densities
      for (int idx = 0; idx < num_species; ++idx){
	cdr_densities[idx] = Max(0.0, vla_cdr_phi[idx][j][i]);
      }

      E   = RealVect(vla_E[0][j][i], vla_E[1][j][i]);
      pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

      const Vector<Real> coeffs = m_physics->compute_cdr_diffusion_coefficients(a_time,
										pos,
										E,
										cdr_densities);

      // Put result in correct palce
      for (int idx = 0; idx < num_species; ++idx){
	vla_cdr_dco[idx][j][i] = coeffs[idx];
      }
    }
  }
  
  // Copy result back to solvers
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    if(solver->is_diffusive()){
      const int idx = solver_it.index();
      FArrayBox& dco = a_diffco_cell[idx]->getFArrayBox();
      dco.setVal(0.0);
      dco.copy(cdr_dco, a_box, idx, a_box, 0, 1);
      a_diffco_cell[idx]->setCoveredCellVal(0.0, 0);
    }
  }
#endif
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast3D(Vector<EBCellFAB*>&       a_diffco_cell,
							    const Vector<EBCellFAB*>& a_cdr_densities,
							    const EBCellFAB&          a_E,
							    const Box                 a_box,
							    const Real                a_dx,
							    const Real                a_time){
#if CH_SPACEDIM==3
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast3D");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast3D" << endl;
  }

  const int comp        = 0;
  const int num_species = m_physics->get_num_cdr_species();
  const RealVect origin = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics
  RealVect     pos, E;
  Vector<Real> cdr_densities(num_species);

  // I need contiguous memory for the nasty that is about to happen. So begin by copying things onto smaller data holders
  FArrayBox cdr_dco(a_box, num_species);
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox Efab(a_box, SpaceDim);
  for (int i = 0; i < num_species; i++){
    cdr_phi.copy(a_cdr_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
  }
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);

  // Pointer offsets
  const IntVect dims = a_box.size();
  const IntVect lo   = a_box.smallEnd();
  const IntVect hi   = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int n2         = dims[2];
  const int offset     = lo[0] + n0*(lo[1] + n1*lo[2]);


  // C style variable-length array conversion magic
  auto vla_cdr_dco = (Real (*__restrict__)[n2][n1][n0]) (cdr_dco.dataPtr() - offset);
  auto vla_cdr_phi = (Real (*__restrict__)[n2][n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_E       = (Real (*__restrict__)[n2][n1][n0]) (Efab.dataPtr()    - offset);

  for (int k = lo[2]; k <= hi[2]; k++){
    for (int j = lo[1]; j <= hi[1]; ++j){
      for (int i = lo[0]; i <= hi[0]; ++i){
      
	// Particle densities
	for (int idx = 0; idx < num_species; ++idx){
	  cdr_densities[idx] = Max(0.0, vla_cdr_phi[idx][k][j][i]);
	}
	
	E   = RealVect(vla_E[0][k][j][i], vla_E[1][k][j][i], vla_E[2][k][j][i]);
	pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

	const Vector<Real> coeffs = m_physics->compute_cdr_diffusion_coefficients(a_time,
										  pos,
										  E,
										  cdr_densities);
	
	// Put result in correct palce
	for (int idx = 0; idx < num_species; ++idx){
	  vla_cdr_dco[idx][k][j][i] = coeffs[idx];
	}
      }
    }
  }
    
  // Copy result back to solvers
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    if(solver->is_diffusive()){
      const int idx = solver_it.index();
      FArrayBox& dco = a_diffco_cell[idx]->getFArrayBox();
      dco.setVal(0.0);
      dco.copy(cdr_dco, a_box, idx, a_box, 0, 1);

      a_diffco_cell[idx]->setCoveredCellVal(0.0, 0);
    }
  }
#endif
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_irreg(Vector<EBCellFAB*>&          a_diffco_cell,
						       const Vector<EBCellFAB*>&    a_cdr_densities,
						       const EBCellFAB&             a_E,
						       const Box                    a_box,
						       const Real                   a_dx,
						       const BaseIVFAB<VoFStencil>& a_interp_stencils,
						       const Real&                  a_time,
						       const int                    a_lvl,
						       const DataIndex&             a_dit){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_cell_irreg");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_cell_irreg" << endl;
  }

  const int comp = 0;
  const int num_species  = m_physics->get_num_cdr_species();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->get_prob_lo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     cdr_densities(num_species);

  VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
      
    const RealVect pos = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);
    const RealVect E   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));

    Vector<Real> cdr_densities(num_species);
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      cdr_densities[idx] = (*a_cdr_densities[idx])(vof, comp);
    }
      
    const Vector<Real> coeffs = m_physics->compute_cdr_diffusion_coefficients(a_time,
									      pos,
									      E,
									      cdr_densities);
      
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->is_diffusive()){
	(*a_diffco_cell[idx])(vof, comp) = coeffs[idx];
      }
    }
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_eb(Vector<EBAMRIVData*>&       a_diffco_eb,
					       const Vector<EBAMRIVData*>& a_cdr_densities,
					       const EBAMRIVData&          a_E,
					       const Real&                 a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_eb(full)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_eb(full)" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->get_num_cdr_species();
  const int finest_level = m_amr->get_finest_level();
  const Real zero        = 0.0;

  Vector<LevelData<BaseIVFAB<Real> >* > diffco(num_species);
  Vector<LevelData<BaseIVFAB<Real> >* > cdr_densities(num_species);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<BaseIVFAB<Real> >* > diffco(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > cdr_densities(num_species);

    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx]        = (*a_diffco_eb[idx])[lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    // Call level versions
    compute_cdr_diffco_eb(diffco, cdr_densities, *a_E[lvl], a_time, lvl);
  }
}


void cdr_plasma_stepper::compute_cdr_diffco_eb(Vector<LevelData<BaseIVFAB<Real> >* >&       a_diffco_eb,
					       const Vector<LevelData<BaseIVFAB<Real> >* >& a_cdr_densities,
					       const LevelData<BaseIVFAB<Real> >&           a_E,
					       const Real&                                  a_time,
					       const int                                    a_lvl){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_eb(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_eb(level)" << endl;
  }

  const int comp = 0;
  const int num_species = m_physics->get_num_cdr_species();

  
  const DisjointBoxLayout& dbl  = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[a_lvl];
  const Real dx                 = m_amr->get_dx()[a_lvl];
  const RealVect origin         = m_amr->get_prob_lo();

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box            = dbl.get(dit());
    const EBISBox& ebisbox   = ebisl[dit()];
    const EBGraph& ebgraph   = ebisbox.getEBGraph();
    const IntVectSet ivs     = ebisbox.getIrregIVS(box);
    const BaseIVFAB<Real>& E = a_E[dit()];

    VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const RealVect cntr = ebisbox.bndryCentroid(vof);
      const RealVect e    = RealVect(D_DECL(E(vof,0), E(vof,1), E(vof,2)));
      const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin) + cntr*dx;
      
      Vector<Real> cdr_densities(num_species);
      for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_cdr_densities[idx])[dit()](vof, 0);
	cdr_densities[idx] = Max(0.0, phi);
      }


      const Vector<Real> diffco = m_physics->compute_cdr_diffusion_coefficients(a_time,
										pos,
										e,
										cdr_densities);
										  
      for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	if(solver_it()->is_diffusive()){
#if 0 // Original code
	  (*a_diffco_eb[idx])[dit()](vof, comp) = diffco[idx];
#else // Debug code
	  (*a_diffco_eb[idx])[dit()](vof, comp) = 0.0;//1.E-5;
#endif
	}
      }
    }
  }
}

void cdr_plasma_stepper::compute_cdr_fluxes(Vector<LevelData<BaseIVFAB<Real> >*>&       a_fluxes,
					    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_fluxes,
					    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_densities,
					    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_velocities,
					    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_gradients,
					    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_rte_fluxes,
					    const LevelData<BaseIVFAB<Real> >&          a_E,
					    const Real&                                 a_time,
					    const int                                   a_lvl){
  CH_TIME("cdr_plasma_stepper::compute_cdr_fluxes(full, level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_fluxes(full, level)" << endl;
  }

  const int num_species  = m_physics->get_num_cdr_species();
  const int num_photons  = m_physics->get_num_rte_species();
  const int comp         = 0;
  const int ncomp        = 1;

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_photons);

  // Grid stuff
  const DisjointBoxLayout& dbl  = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[a_lvl];
  const ProblemDomain& domain   = m_amr->get_domains()[a_lvl];
  const Real dx                 = m_amr->get_dx()[a_lvl];
  const MFLevelGrid& mflg       = *(m_amr->get_mflg(m_realm)[a_lvl]);
  const RealVect origin         = m_amr->get_prob_lo();

  // Patch loop
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box              = dbl.get(dit());
    const EBISBox& ebisbox      = ebisl[dit()];
    const EBGraph& ebgraph      = ebisbox.getEBGraph();
    const IntVectSet& diel_ivs  = mflg.interface_region(box, dit());
    const IntVectSet& elec_ivs  = ebisbox.getIrregIVS(box) - diel_ivs;

    // Loop over conductor cells
    for (VoFIterator vofit(elec_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();


      // Define the electric field
      const BaseIVFAB<Real>& E = a_E[dit()];
      const RealVect centroid  = ebisbox.bndryCentroid(vof);
      const RealVect normal    = ebisbox.normal(vof);
      const RealVect e         = RealVect(D_DECL(E(vof,0), E(vof,1), E(vof,2)));
      const RealVect pos       = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin) + centroid*dx;

      // Build ion densities and velocities
      for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](vof,comp);
	extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](vof,comp);
	extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](vof,comp);
	extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](vof,comp);
      }

      // Build photon intensities
      for (rte_iterator<rte_solver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[dit()](vof,comp);
      }

      const Vector<Real> fluxes = m_physics->compute_cdr_electrode_fluxes(a_time,
									  pos,
									  normal,
									  e,
									  extrap_cdr_densities,
									  extrap_cdr_velocities,
									  extrap_cdr_gradients,
									  extrap_rte_fluxes,
									  extrap_cdr_fluxes);
	
      // Put the fluxes in their respective place
      for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	(*a_fluxes[idx])[dit()](vof, comp) = fluxes[idx];
      }
    }

    // Loop over dielectric cells
    for (VoFIterator vofit(diel_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      // Define the electric field
      const BaseIVFAB<Real>& E = a_E[dit()];
      const RealVect centroid  = ebisbox.bndryCentroid(vof);
      const RealVect normal    = ebisbox.normal(vof);
      const RealVect e         = RealVect(D_DECL(E(vof,0), E(vof,1), E(vof,2)));
      const RealVect pos       = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin) + centroid*dx;
      const Real     time      = 0.0;

      // Build ion densities and velocities
      for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](vof,comp);
	extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](vof,comp);
	extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](vof,comp);
	extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](vof,comp);
      }

      // Build photon intensities
      for (rte_iterator<rte_solver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[dit()](vof,comp);
      }

      const Vector<Real> fluxes = m_physics->compute_cdr_dielectric_fluxes(a_time,
									   pos,
									   normal,
									   e,
									   extrap_cdr_densities,
									   extrap_cdr_velocities,
									   extrap_cdr_gradients,
									   extrap_rte_fluxes,
									   extrap_cdr_fluxes);
	
      // Put the fluxes in their respective place
      for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	(*a_fluxes[idx])[dit()](vof, comp) = fluxes[idx];
      }
    }
  }
}

void cdr_plasma_stepper::compute_cdr_fluxes(Vector<EBAMRIVData*>&       a_fluxes,
					    const Vector<EBAMRIVData*>& a_extrap_cdr_fluxes,
					    const Vector<EBAMRIVData*>& a_extrap_cdr_densities,
					    const Vector<EBAMRIVData*>& a_extrap_cdr_velocities,
					    const Vector<EBAMRIVData*>& a_extrap_cdr_gradients,
					    const Vector<EBAMRIVData*>& a_extrap_rte_fluxes,
					    const EBAMRIVData&          a_E,
					    const Real&                 a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_fluxes(full)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_fluxes(full)" << endl;
  }

  const int num_species  = m_physics->get_num_cdr_species();
  const int num_photons  = m_physics->get_num_rte_species();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){


    Vector<LevelData<BaseIVFAB<Real> >* > fluxes(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_fluxes(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_densities(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_velocities(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_gradients(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_rte_fluxes(num_photons);
    
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      fluxes[idx]                = (*a_fluxes[idx])[lvl];
      extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[lvl];
      extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[lvl];
      extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[lvl];
      extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[lvl];
    }


    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[lvl];
    }

    // Call the level version
    compute_cdr_fluxes(fluxes, extrap_cdr_fluxes, extrap_cdr_densities, extrap_cdr_velocities,
		       extrap_cdr_gradients, extrap_rte_fluxes, *a_E[lvl], a_time, lvl);
  }
}

void cdr_plasma_stepper::compute_cdr_domain_fluxes(Vector<EBAMRIFData*>&       a_fluxes,
						   const Vector<EBAMRIFData*>& a_extrap_cdr_fluxes,
						   const Vector<EBAMRIFData*>& a_extrap_cdr_densities,
						   const Vector<EBAMRIFData*>& a_extrap_cdr_velocities,
						   const Vector<EBAMRIFData*>& a_extrap_cdr_gradients,
						   const Vector<EBAMRIFData*>& a_extrap_rte_fluxes,
						   const EBAMRIFData&          a_E,
						   const Real&                 a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_domain_fluxes(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_domain_fluxes(level)" << endl;
  }

  const int num_species  = m_physics->get_num_cdr_species();
  const int num_photons  = m_physics->get_num_rte_species();
  const int finest_level = m_amr->get_finest_level();

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_photons);

  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<DomainFluxIFFAB>* > fluxes(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_fluxes(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_densities(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_velocities(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_gradients(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_rte_fluxes(num_photons);
    
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      fluxes[idx]                = (*a_fluxes[idx])[lvl];
      extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[lvl];
      extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[lvl];
      extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[lvl];
      extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[lvl];
    }

    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[lvl];
    }

    // Call the level version
    compute_cdr_domain_fluxes(fluxes, extrap_cdr_fluxes, extrap_cdr_densities, extrap_cdr_velocities,
			      extrap_cdr_gradients, extrap_rte_fluxes, *a_E[lvl], a_time, lvl);
    
  }
}

void cdr_plasma_stepper::compute_cdr_domain_fluxes(Vector<LevelData<DomainFluxIFFAB>*>        a_fluxes,
						   const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_fluxes,
						   const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_densities,
						   const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_velocities,
						   const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_gradients,
						   const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_rte_fluxes,
						   const LevelData<DomainFluxIFFAB>&          a_E,
						   const Real&                                a_time,
						   const int                                  a_lvl){
  CH_TIME("cdr_plasma_stepper::compute_cdr_domain_fluxes(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_domain_fluxes(level)" << endl;
  }

  const int num_species  = m_physics->get_num_cdr_species();
  const int num_photons  = m_physics->get_num_rte_species();
  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_photons);

  const DisjointBoxLayout& dbl  = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[a_lvl];
  const ProblemDomain& domain   = m_amr->get_domains()[a_lvl];
  const Real dx                 = m_amr->get_dx()[a_lvl];
  const RealVect origin         = m_amr->get_prob_lo();

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box         = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    
    const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
    
    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){
	const IntVectSet& ivs  = (*a_fluxes[0])[dit()](dir, sit()).getIVS();

	for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  const RealVect pos    = EBArith::getFaceLocation(face, dx*RealVect::Unit, origin);
	    
	  // Define the electric field
	  const RealVect E = RealVect(D_DECL((a_E)[dit()](dir, sit())(face,0),
					     (a_E)[dit()](dir, sit())(face,1),
					     (a_E)[dit()](dir, sit())(face,2)));

	  // Ion densities. 
	  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.index();
	    extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](dir, sit())(face,comp);
	  }

	  // Photon fluxes
	  for (rte_iterator<rte_solver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.index();
	    extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[dit()](dir, sit())(face,comp);
	  }

	  // Call cdr_plasma_physics
	  const Vector<Real> fluxes = m_physics->compute_cdr_domain_fluxes(a_time,
									   pos,
									   dir,
									   sit(),
									   E,
									   extrap_cdr_densities,
									   extrap_cdr_velocities,
									   extrap_cdr_gradients,
									   extrap_rte_fluxes,
									   extrap_cdr_fluxes);

	  // Put fluxes where they belong
	  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.index();
	    (*a_fluxes[idx])[dit()](dir, sit())(face, comp) = fluxes[idx];
	  }
	}
      }
    }
  }
}

void cdr_plasma_stepper::compute_gradients_at_eb(Vector<EBAMRIVData*>&         a_grad,
						 const phase::which_phase&     a_phase,
						 const Vector<EBAMRCellData*>& a_phi){
  CH_TIME("cdr_plasma_stepper::compute_gradients_at_eb");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_gradients_at_eb" << endl;
  }

  CH_assert(a_grad.size() == a_phi.size());

  EBAMRIVData   eb_gradient;
  EBAMRCellData gradient;
  m_amr->allocate(eb_gradient, m_realm, a_phase, SpaceDim);
  m_amr->allocate(gradient,    m_realm, a_phase, SpaceDim);

  for (int i = 0; i < a_phi.size(); i++){
    EBAMRIVData& grad_density    = *a_grad[i];
    const EBAMRCellData& density = *a_phi[i];

    CH_assert(grad_density[0]->nComp() == 1);
    CH_assert(density[0]->nComp()      == 1);
    
    m_amr->compute_gradient(gradient, density, m_realm, a_phase);       // Compute cell-centered gradient
    m_amr->average_down(gradient, m_realm, a_phase);                    // Average down - shouldn't be necesasry
    m_amr->interp_ghost(gradient, m_realm, a_phase);                    // Interpolate ghost cells (have to do this before interp)
    this->extrapolate_to_eb(eb_gradient, a_phase, gradient);            // Extrapolate to EB
    this->project_flux(grad_density, eb_gradient);                      // Project onto EB
  }
}

void cdr_plasma_stepper::compute_gradients_at_domain_faces(Vector<EBAMRIFData*>&         a_grad,
							   const phase::which_phase&     a_phase,
							   const Vector<EBAMRCellData*>& a_phi){
  CH_TIME("cdr_plasma_stepper::compute_gradients_at_eb");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_gradients_at_eb" << endl;
  }

  CH_assert(a_grad.size() == a_phi.size());

  EBAMRIFData   domain_gradient;
  EBAMRCellData gradient;
  m_amr->allocate(domain_gradient, m_realm, a_phase, SpaceDim);
  m_amr->allocate(gradient,        m_realm, a_phase, SpaceDim);

  for (int i = 0; i < a_phi.size(); i++){
    EBAMRIFData& grad_density    = *a_grad[i];
    const EBAMRCellData& density = *a_phi[i];

    CH_assert(grad_density[0]->nComp() == 1);
    CH_assert(density[0]->nComp()      == 1);
    
    m_amr->compute_gradient(gradient, density, m_realm, a_phase);                         
    m_amr->average_down(gradient, m_realm, a_phase);                             
    m_amr->interp_ghost(gradient, m_realm, a_phase);                             
    
    this->extrapolate_to_domain_faces(domain_gradient, a_phase, gradient);  // Extrapolate to EB
    this->project_domain(grad_density, domain_gradient);                    // Project normal compoent
  }
}

void cdr_plasma_stepper::extrapolate_vector_to_domain_faces(EBAMRIFData&             a_extrap,
							    const phase::which_phase a_phase,
							    const EBAMRCellData&     a_data){
  CH_TIME("cdr_plasma_stepper::extrapolate_to_domain_faces");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_to_domain_faces" << endl;
  }

  CH_assert(a_extrap[0]->nComp() == 1);
  CH_assert(a_data[0]->nComp()   == SpaceDim);

  EBAMRIFData domain_vector;
  m_amr->allocate(domain_vector, m_realm, a_phase, SpaceDim);
  this->extrapolate_to_domain_faces(domain_vector, a_phase, a_data);
  this->project_domain(a_extrap, domain_vector);
}

void cdr_plasma_stepper::extrapolate_vector_to_domain_faces(Vector<EBAMRIFData*>&         a_extrap,
							    const phase::which_phase      a_phase,
							    const Vector<EBAMRCellData*>& a_data){
  CH_TIME("cdr_plasma_stepper::extrapolate_vector_to_domain_faces");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_vector_to_domain_faces" << endl;
  }

  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const int idx = solver_it.index();
    extrapolate_vector_to_domain_faces(*a_extrap[idx], a_phase, *a_data[idx]);
  }
}

void cdr_plasma_stepper::extrapolate_velo_to_domain_faces(Vector<EBAMRIFData*>&         a_extrap,
							  const phase::which_phase      a_phase,
							  const Vector<EBAMRCellData*>& a_velocities){
  CH_TIME("cdr_plasma_stepper::extrapolate_velo_to_domain_faces");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_velo_to_domain_faces)" << endl;
  }

  //  for(int i = 0; i < a_extrap.size(); i++){
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->is_mobile()){
      extrapolate_vector_to_domain_faces(*a_extrap[idx], a_phase, *a_velocities[idx]);
    }
    else{
      data_ops::set_value(*a_extrap[idx], 0.0);
    }
  }
}

void cdr_plasma_stepper::compute_cdr_velocities(Vector<EBAMRCellData*>&       a_velocities,
						const Vector<EBAMRCellData*>& a_cdr_densities,
						const EBAMRCellData&          a_E,
						const Real&                   a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_velocities(full)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_velocities(full)" << endl;
  }
  
  const int num_species  = m_physics->get_num_cdr_species();
  const int finest_level = m_amr->get_finest_level();

  // Interpolate E to centroids
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, phase::gas, SpaceDim);
  data_ops::copy(E, a_E);
  m_amr->interpolate_to_centroids(E, m_realm, phase::gas);

  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<EBCellFAB>* > velocities(num_species);
    Vector<LevelData<EBCellFAB>* > cdr_densities(num_species);

    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      velocities[idx]    = (*a_velocities[idx])[lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    compute_cdr_velocities(velocities, cdr_densities, *E[lvl], lvl, a_time);
  }

  // Average down and interpolate ghost cells
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    if(solver_it()->is_mobile()){
      m_amr->average_down(*a_velocities[idx], m_realm, m_cdr->get_phase()); 
      m_amr->interp_ghost(*a_velocities[idx], m_realm, m_cdr->get_phase()); 
    }
  }
}

void cdr_plasma_stepper::compute_cdr_velocities(Vector<LevelData<EBCellFAB> *>&       a_velocities,
						const Vector<LevelData<EBCellFAB> *>& a_cdr_densities,
						const LevelData<EBCellFAB> &          a_E,
						const int                             a_lvl,
						const Real&                           a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_velocities(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_velocities(level)" << endl;
  }

  const phase::which_phase cdr_phase = m_cdr->get_phase();
  const int finest_level             = m_amr->get_finest_level();

  const DisjointBoxLayout& dbl  = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->get_ebisl(m_realm, cdr_phase)[a_lvl];
  const Real dx                 = m_amr->get_dx()[a_lvl];

  const int num_species = m_physics->get_num_cdr_species();
    
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    Vector<EBCellFAB*> vel(num_species);
    Vector<EBCellFAB*> phi(num_species);;
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->is_mobile()){
	vel[idx] = &(*a_velocities[idx])[dit()];
      }
      phi[idx] = &(*a_cdr_densities[idx])[dit()];
    }

    // Separate calls for regular and irregular
#if USE_FAST_VELOCITIES
    compute_cdr_velocities_reg_fast(vel,   phi, a_E[dit()], dbl.get(dit()), a_time, dx);
#else
    compute_cdr_velocities_reg(vel,   phi, a_E[dit()], dbl.get(dit()), a_time, dx);
#endif
    compute_cdr_velocities_irreg(vel, phi, a_E[dit()], dbl.get(dit()), a_time, dx, a_lvl, dit());
  }
}

void cdr_plasma_stepper::compute_cdr_velocities_reg(Vector<EBCellFAB*>&       a_velocities,
						    const Vector<EBCellFAB*>& a_cdr_densities,
						    const EBCellFAB&          a_E,
						    const Box&                a_box,
						    const Real&               a_time,
						    const Real&               a_dx){
  CH_TIME("cdr_plasma_stepper::compute_cdr_velocities_reg");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_velocities_reg" << endl;
  }

  const int comp         = 0;
  const RealVect origin  = m_amr->get_prob_lo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();


  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv    = bit();
    const RealVect pos  = origin + a_dx*iv;
    const RealVect e    = RealVect(D_DECL(E(iv, 0), E(iv, 1), E(iv, 2)));


    // Get densities
    Vector<Real> cdr_densities;
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      cdr_densities.push_back((*a_cdr_densities[idx]).getSingleValuedFAB()(iv, comp));
    }

    // Compute velocities
    Vector<RealVect> velocities = m_physics->compute_cdr_velocities(a_time, pos, e, cdr_densities);

    // Put velocities in the appropriate place. 
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver = solver_it();
      const int idx = solver_it.index();
      if(solver->is_mobile()){
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_velocities[idx]).getSingleValuedFAB()(iv, dir) = velocities[idx][dir];
	}
      }
    }
  }


  // Covered is bogus.
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->is_mobile()){
      const int idx = solver_it.index();
      for (int dir = 0; dir < SpaceDim; dir++){
	a_velocities[idx]->setCoveredCellVal(0.0, dir);
      }
    }
  }
}

void cdr_plasma_stepper::compute_cdr_velocities_reg_fast(Vector<EBCellFAB*>&       a_velocities,
							 const Vector<EBCellFAB*>& a_cdr_densities,
							 const EBCellFAB&          a_E,
							 const Box&                a_box,
							 const Real&               a_time,
							 const Real&               a_dx){
#if CH_SPACEDIM==2
  compute_cdr_velocities_reg_fast2D(a_velocities, a_cdr_densities, a_E, a_box, a_time, a_dx);
#elif CH_SPACEDIM==3
  compute_cdr_velocities_reg_fast3D(a_velocities, a_cdr_densities, a_E, a_box, a_time, a_dx);
#endif
}

void cdr_plasma_stepper::compute_cdr_velocities_reg_fast2D(Vector<EBCellFAB*>&       a_velocities,
							   const Vector<EBCellFAB*>& a_cdr_densities,
							   const EBCellFAB&          a_E,
							   const Box&                a_box,
							   const Real&               a_time,
							   const Real&               a_dx){
  CH_TIME("cdr_plasma_stepper::compute_cdr_velocities_reg_fast2D");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_velocities_reg_fast2D" << endl;
  }
  
#if CH_SPACEDIM==2
  const int num_species  = m_physics->get_num_cdr_species();
  const int comp         = 0;
  const RealVect origin  = m_amr->get_prob_lo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // I need contiguous memory for the nasty stuff that is about to happen, so begin by copying things onto smaller
  // data holders
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox cdr_vx(a_box,  num_species);
  FArrayBox cdr_vy(a_box,  num_species);
  for (int idx = 0; idx < num_species; idx++){
    cdr_phi.copy(a_cdr_densities[idx]->getFArrayBox(), a_box, 0, a_box, idx, 1);
  }

  FArrayBox Efab(a_box, SpaceDim);
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);


  // Pointer offsets
  const IntVect dims   = a_box.size();
  const IntVect lo     = a_box.smallEnd();
  const IntVect hi     = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int offset     = lo[0] + n0*lo[1];

  // C style variable-length array conversion magic
  auto vla_cdr_phi = (Real (*__restrict__)[n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_cdr_vx  = (Real (*__restrict__)[n1][n0]) (cdr_vx.dataPtr()  - offset);
  auto vla_cdr_vy  = (Real (*__restrict__)[n1][n0]) (cdr_vy.dataPtr()  - offset);
  auto vla_E       = (Real (*__restrict__)[n1][n0]) (Efab.dataPtr()    - offset);
  
  for (int j = lo[1]; j <= hi[1]; ++j){
    for (int i = lo[0]; i <= hi[0]; ++i){

      const RealVect pos = origin + RealVect(D_DECL(i,j,k))*a_dx;
      const RealVect E   = RealVect(vla_E[0][j][i], vla_E[1][j][i]);

      Vector<Real> cdr_densities(num_species);
      for (int idx = 0; idx < num_species; idx++){
	cdr_densities[idx] = Max(0.0, vla_cdr_phi[idx][j][i]);
      }

      const Vector<RealVect> velocities = m_physics->compute_cdr_velocities(a_time, pos, E, cdr_densities);

      // Put result in correct palce
      for (int idx = 0; idx < num_species; ++idx){
	vla_cdr_vx[idx][j][i] = velocities[idx][0];
	vla_cdr_vy[idx][j][i] = velocities[idx][1];
      }
    }
  }

  // Put the results back into solvers
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->is_mobile()){
      FArrayBox& v = a_velocities[idx]->getFArrayBox();

      v.copy(cdr_vx, a_box, idx, a_box, 0, 1);
      v.copy(cdr_vy, a_box, idx, a_box, 1, 1);

      for (int dir = 0; dir < SpaceDim; dir++){
	a_velocities[idx]->setCoveredCellVal(0.0, dir);
      }
    }
  }
#endif
}



void cdr_plasma_stepper::compute_cdr_velocities_reg_fast3D(Vector<EBCellFAB*>&       a_velocities,
							   const Vector<EBCellFAB*>& a_cdr_densities,
							   const EBCellFAB&          a_E,
							   const Box&                a_box,
							   const Real&               a_time,
							   const Real&               a_dx){
  CH_TIME("cdr_plasma_stepper::compute_cdr_velocities_reg_fast3D");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_velocities_reg_fast3D" << endl;
  }

#if CH_SPACEDIM==3
  const int num_species  = m_physics->get_num_cdr_species();
  const int comp         = 0;
  const RealVect origin  = m_amr->get_prob_lo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // I need contiguous memory for the nasty stuff that is about to happen, so begin by copying things onto smaller
  // data holders
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox cdr_vx(a_box,  num_species);
  FArrayBox cdr_vy(a_box,  num_species);
  FArrayBox cdr_vz(a_box,  num_species);
  for (int idx = 0; idx < num_species; idx++){
    cdr_phi.copy(a_cdr_densities[idx]->getFArrayBox(), a_box, 0, a_box, idx, 1);
  }

  FArrayBox Efab(a_box, SpaceDim);
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);


  // Pointer offsets
  const IntVect dims   = a_box.size();
  const IntVect lo     = a_box.smallEnd();
  const IntVect hi     = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int n2         = dims[2];
  const int offset     = lo[0] + n0*(lo[1] + n1*lo[2]);

  auto vla_cdr_phi = (Real (*__restrict__)[n2][n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_cdr_vx  = (Real (*__restrict__)[n2][n1][n0]) (cdr_vx.dataPtr()  - offset);
  auto vla_cdr_vy  = (Real (*__restrict__)[n2][n1][n0]) (cdr_vy.dataPtr()  - offset);
  auto vla_cdr_vz  = (Real (*__restrict__)[n2][n1][n0]) (cdr_vz.dataPtr()  - offset);
  auto vla_E       = (Real (*__restrict__)[n2][n1][n0]) (Efab.dataPtr()    - offset);
  
  for (int k = lo[2]; k <= hi[2]; k++){
    for (int j = lo[1]; j <= hi[1]; ++j){
      for (int i = lo[0]; i <= hi[0]; ++i){

	const RealVect pos = origin + RealVect(D_DECL(i,j,k))*a_dx;
	const RealVect E   = RealVect(vla_E[0][k][j][i], vla_E[1][k][j][i], vla_E[2][k][j][i]);

	Vector<Real> cdr_densities(num_species);
	for (int idx = 0; idx < num_species; idx++){
	  cdr_densities[idx] = Max(0.0, vla_cdr_phi[idx][k][j][i]);
	}

	const Vector<RealVect> velocities = m_physics->compute_cdr_velocities(a_time, pos, E, cdr_densities);

	// Put result in correct palce
	for (int idx = 0; idx < num_species; ++idx){
	  vla_cdr_vx[idx][k][j][i] = velocities[idx][0];
	  vla_cdr_vy[idx][k][j][i] = velocities[idx][1];
	  vla_cdr_vz[idx][k][j][i] = velocities[idx][2];
	}
      }
    }
  }
  
  // Put the results back into solvers
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->is_mobile()){
      FArrayBox& v = a_velocities[idx]->getFArrayBox();
      
      v.copy(cdr_vx, a_box, 0, a_box, 0, 1);
      v.copy(cdr_vy, a_box, 0, a_box, 1, 1);
      v.copy(cdr_vz, a_box, 0, a_box, 2, 1);

      for (int dir = 0; dir < SpaceDim; dir++){
	a_velocities[idx]->setCoveredCellVal(0.0, dir);
      }      
    }
  }
#endif
}


void cdr_plasma_stepper::compute_cdr_velocities_irreg(Vector<EBCellFAB*>&       a_velocities,
						      const Vector<EBCellFAB*>& a_cdr_densities,
						      const EBCellFAB&          a_E,
						      const Box&                a_box,
						      const Real&               a_time,
						      const Real&               a_dx,
						      const int                 a_lvl,
						      const DataIndex&          a_dit){
  CH_TIME("cdr_plasma_stepper::compute_cdr_velocities_irreg");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_velocities_irreg" << endl;
  }

  const int comp         = 0;
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->get_prob_lo();

  VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);

    // Get densities
    Vector<Real> cdr_densities;
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      cdr_densities.push_back((*a_cdr_densities[idx])(vof, comp));
    }
    
    // Compute velocities
    const Vector<RealVect> velocities = m_physics->compute_cdr_velocities(a_time, pos, e, cdr_densities);

    // Put velocities in the appropriate place. 
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      if(solver_it()->is_mobile()){
	const int idx = solver_it.index();
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_velocities[idx])(vof, dir) = velocities[idx][dir];
	}
      }
    }
  }
}

void cdr_plasma_stepper::pre_regrid(const int a_lmin, const int a_finest_level){
  CH_TIME("cdr_plasma_stepper::pre_regrid");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::pre_regrid" << endl;
  }

  // Solvers do pre-regridding shit. 
  m_cdr->pre_regrid(a_lmin, a_finest_level);
  m_poisson->pre_regrid(a_lmin, a_finest_level);
  m_rte->pre_regrid(a_lmin, a_finest_level);
  m_sigma->pre_regrid(a_lmin, a_finest_level);
}

void cdr_plasma_stepper::pre_regrid_internals(const int a_lbase, const int a_finest_level){
  CH_TIME("cdr_plasma_stepper::pre_regrid_internals");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::pre_regrid_internals" << endl;
  }
}

void cdr_plasma_stepper::compute_cdr_velocities(){
  CH_TIME("cdr_plasma_stepper::compute_cdr_velocities()");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_velocities()" << endl;
  }

  // Compute the electric field (again)
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, m_cdr->get_phase(), SpaceDim);
  this->compute_E(E, m_cdr->get_phase(), m_poisson->get_state());

  Vector<EBAMRCellData*> states     = m_cdr->get_states();
  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();

  this->compute_cdr_velocities(velocities, states, E, m_time);
}

void cdr_plasma_stepper::compute_cdr_diffusion(){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffusion");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffusion" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_physics->get_num_cdr_species();

  EBAMRCellData E_cell;
  EBAMRIVData   E_eb;
  m_amr->allocate(E_cell, m_realm, m_cdr->get_phase(), SpaceDim);
  m_amr->allocate(E_eb,   m_realm, m_cdr->get_phase(), SpaceDim);
  
  this->compute_E(E_cell, m_cdr->get_phase(), m_poisson->get_state());
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);

  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();

  cdr_plasma_stepper::compute_cdr_diffusion(E_cell, E_eb);
#if 0

  // Extrapolate states to the EB
  Vector<EBAMRIVData*> cdr_extrap(num_species);
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    cdr_extrap[idx] = new EBAMRIVData();  // This must be deleted
    m_amr->allocate(*cdr_extrap[idx], m_realm, m_cdr->get_phase(), ncomp);

    const irreg_amr_stencil<eb_centroid_interp>& stencil = m_amr->get_eb_centroid_interp_stencils(m_realm, m_cdr->get_phase());
    stencil.apply(*cdr_extrap[idx], *cdr_states[idx]);
  }
  
  Vector<EBAMRFluxData*> diffco_face = m_cdr->get_diffco_face();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->get_diffco_eb();

  this->compute_cdr_diffco_face(diffco_face, cdr_states, E_cell, m_time);
  this->compute_cdr_diffco_eb(diffco_eb,     cdr_extrap, E_eb,   m_time);

  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    delete cdr_extrap[idx];
  }
#endif
}

void cdr_plasma_stepper::compute_cdr_diffusion(const EBAMRCellData& a_E_cell, const EBAMRIVData& a_E_eb){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffusion");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffusion" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_physics->get_num_cdr_species();

  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();

  // Extrapolate states to the EB
  Vector<EBAMRIVData*> cdr_extrap(num_species);
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    cdr_extrap[idx] = new EBAMRIVData();  // This must be deleted
    m_amr->allocate(*cdr_extrap[idx], m_realm, m_cdr->get_phase(), ncomp);

    const irreg_amr_stencil<eb_centroid_interp>& stencil = m_amr->get_eb_centroid_interp_stencils(m_realm, m_cdr->get_phase());
    stencil.apply(*cdr_extrap[idx], *cdr_states[idx]);
  }
  
  Vector<EBAMRFluxData*> diffco_face = m_cdr->get_diffco_face();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->get_diffco_eb();

  this->compute_cdr_diffco_face(diffco_face, cdr_states, a_E_cell, m_time);
  this->compute_cdr_diffco_eb(diffco_eb,     cdr_extrap, a_E_eb,   m_time);

  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->deallocate(*cdr_extrap[idx]);
    delete cdr_extrap[idx];
  }
}

void cdr_plasma_stepper::compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential){
  CH_TIME("cdr_plasma_stepper::compute_E(mfamrcell, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_E(mfamrcell, mfamrcell)" << endl;
  }

  m_amr->compute_gradient(a_E, a_potential, m_realm);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E, m_realm);
  m_amr->interp_ghost(a_E, m_realm);
}

void cdr_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase){
  CH_TIME("cdr_plasma_stepper::compute_E(ebamrcell, phase)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_E(ebamrcell, phase)" << endl;
  }

  this->compute_E(a_E, a_phase, m_poisson->get_state());
}

void cdr_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential){
  CH_TIME("cdr_plasma_stepper::compute_E(ebamrcell, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_E(ebamrcell, phase, mfamrcell)" << endl;
  }

  EBAMRCellData pot_gas;
  m_amr->allocate_ptr(pot_gas);
  m_amr->alias(pot_gas, a_phase, a_potential);

  m_amr->compute_gradient(a_E, pot_gas, m_realm, a_phase);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E, m_realm, a_phase);
  m_amr->interp_ghost(a_E, m_realm, a_phase);
}

void cdr_plasma_stepper::compute_E(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("cdr_plasma_stepper::compute_E(ebamrflux, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_E(ebamrflux, phase, ebamrcell)" << endl;
  }

  CH_assert(a_E_face[0]->nComp() == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, a_phase)[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& E_cell = (*a_E_cell[lvl])[dit()];
      const EBISBox& ebisbox  = ebisl[dit()];
      const EBGraph& ebgraph  = ebisbox.getEBGraph();
      const Box& box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
      	EBFaceFAB& E_face = (*a_E_face[lvl])[dit()][dir];
	E_face.setVal(0.0);

      	EBLevelDataOps::averageCellToFace(E_face,
      					  E_cell,
      					  ebgraph,
      					  box,
      					  0,
      					  dir,
      					  domain,
      					  dir,
      					  dir);
      }

    }
    //    data_ops::average_cell_to_face_allcomps(*a_E_face[lvl], *a_E_cell[lvl], m_amr->get_domains()[lvl]);
    a_E_face[lvl]->exchange();
  }
}

void cdr_plasma_stepper::compute_E(EBAMRIVData& a_E_eb, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("cdr_plasma_stepper::compute_E(ebamriv, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_E(ebamriv, phase, ebamrcell)" << endl;
  }

  CH_assert(a_E_eb[0]->nComp()   == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const irreg_amr_stencil<eb_centroid_interp>& interp_stencil = m_amr->get_eb_centroid_interp_stencils(m_realm, a_phase);
  interp_stencil.apply(a_E_eb, a_E_cell);
}

void cdr_plasma_stepper::compute_Emax(Real& a_Emax, const phase::which_phase a_phase){
  CH_TIME("cdr_plasma_stepper::compute_Emax(Real, phase)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_Emax(Real, phase)" << endl;
  }

  EBAMRCellData E;
  m_amr->allocate(E, m_realm, a_phase, SpaceDim);

  this->compute_E(E, a_phase, m_poisson->get_state());
  m_amr->interpolate_to_centroids(E, m_realm, a_phase);

  Real max, min;
  data_ops::get_max_min_norm(max, min, E);

  a_Emax = max;
}

void cdr_plasma_stepper::compute_charge_flux(EBAMRIVData& a_flux, Vector<EBAMRIVData*>& a_cdr_fluxes){
  CH_TIME("cdr_plasma_stepper::compute_charge_flux");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_charge_flux" << endl;
  }

  MayDay::Abort("cdr_plasma_stepper::compute_charge_flux - I'm suspecting that this is deprecated code which is no longer used");

  data_ops::set_value(a_flux, 0.0);

  for (cdr_iterator<cdr_solver> solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<cdr_species>& spec      = solver_it.get_species();
    const EBAMRIVData& solver_flux          = *a_cdr_fluxes[solver_it.index()];

    data_ops::incr(a_flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  m_sigma->reset_cells(a_flux);
}

void cdr_plasma_stepper::compute_extrapolated_fluxes(Vector<EBAMRIVData*>&        a_fluxes,
						     const Vector<EBAMRCellData*> a_densities,
						     const Vector<EBAMRCellData*> a_velocities,
						     const phase::which_phase     a_phase){
  CH_TIME("cdr_plasma_stepper::compute_extrapolated_fluxes");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_extrapolated_fluxes" << endl;
  }

#if 0 // This code computes the cell-centered flux which is then extrapolated to the EB. Please don't remove this code. 
  EBAMRCellData cell_flux;
  EBAMRIVData   eb_flux;

  m_amr->allocate(cell_flux, m_realm, a_phase, SpaceDim);
  m_amr->allocate(eb_flux,   m_realm, a_phase, SpaceDim);

  for (int i = 0; i < a_fluxes.size(); i++){
    this->compute_flux(cell_flux, *a_densities[i], *a_velocities[i]);
    
    m_amr->average_down(cell_flux, m_realm, a_phase);
    m_amr->interp_ghost(cell_flux, m_realm, a_phase);

    this->extrapolate_to_eb(eb_flux, a_phase, cell_flux);
    
    this->project_flux(*a_fluxes[i], eb_flux);
  }
#else // New way of doing this. Extrapolate everything to the BC first. Then compute the flux and project it. We do this because
      // extrapolation stencils may have negative weights, and v*n may therefore nonphysically change sign. Better to compute
      // F = v_extrap*Max(0.0, phi_extrap) since we expect v to be "smooth" and phi_extrap to be a noisy bastard
  EBAMRIVData eb_flx; 
  EBAMRIVData eb_vel;
  EBAMRIVData eb_phi;

  m_amr->allocate(eb_flx, m_realm, a_phase, SpaceDim);
  m_amr->allocate(eb_vel, m_realm, a_phase, SpaceDim);
  m_amr->allocate(eb_phi, m_realm, a_phase, 1);

  const irreg_amr_stencil<eb_centroid_interp>& interp_stencils = m_amr->get_eb_centroid_interp_stencils(m_realm, a_phase);

  //  for (int i = 0; i < a_fluxes.size(); i++){
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->is_mobile()){
      interp_stencils.apply(eb_vel, *a_velocities[idx]);
      interp_stencils.apply(eb_phi, *a_densities[idx]);

      data_ops::floor(eb_phi, 0.0);

      data_ops::set_value(eb_flx, 0.0);
      data_ops::incr(eb_flx, eb_vel, 1.0);
      data_ops::multiply_scalar(eb_flx, eb_phi);

      this->project_flux(*a_fluxes[idx], eb_flx);

      m_amr->average_down(*a_fluxes[idx], m_realm, a_phase);
    }
    else{
      data_ops::set_value(*a_fluxes[idx], 0.0);
    }
  }
#endif
}

void cdr_plasma_stepper::compute_extrapolated_velocities(Vector<EBAMRIVData*>&        a_ebvelo,
							 const Vector<EBAMRCellData*> a_cell_vel,
							 const phase::which_phase     a_phase){
  CH_TIME("cdr_plasma_stepper::compute_extrapolated_velocities");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_extrapolated_velocities" << endl;
  }

  EBAMRIVData scratch;
  m_amr->allocate(scratch, m_realm, a_phase, SpaceDim);

  //  for (int i = 0; i < a_ebvelo.size(); i++){
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->is_mobile()){
      this->extrapolate_to_eb(scratch, a_phase, *a_cell_vel[idx]);
      this->project_flux(*a_ebvelo[idx], scratch);
    }
  }
}

void cdr_plasma_stepper::compute_extrapolated_domain_fluxes(Vector<EBAMRIFData*>&        a_fluxes,
							    const Vector<EBAMRCellData*> a_densities,
							    const Vector<EBAMRCellData*> a_velocities,
							    const phase::which_phase     a_phase){
  CH_TIME("cdr_plasma_stepper::compute_extrapolated_domain_fluxes");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_extrapolated_domain_fluxes" << endl;
  }

  EBAMRCellData cell_flux;
  EBAMRIFData domain_flux;

  m_amr->allocate(cell_flux,   m_realm, a_phase, SpaceDim);
  m_amr->allocate(domain_flux, m_realm, a_phase, SpaceDim);

  //  for (int i = 0; i < a_fluxes.size(); i++){
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->is_mobile()){
      data_ops::copy(cell_flux, *a_velocities[idx]);
      data_ops::multiply_scalar(cell_flux, *a_densities[idx]); // cell_flux = n*v

      // Extrapolate cell-centered to domain faces
      this->extrapolate_to_domain_faces(domain_flux, a_phase, cell_flux);

      // Project normal component onto domain face
      this->project_domain(*a_fluxes[idx], domain_flux);
    }
    else{
      data_ops::set_value(*a_fluxes[idx], 0.0);
    }
  }
}

void cdr_plasma_stepper::compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_density, const EBAMRCellData& a_velocity){
  CH_TIME("cdr_plasma_stepper::compute_flux(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_flux(amr)" << endl;
  }
  
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    CH_assert(a_flux[lvl]->nComp()     == SpaceDim);
    CH_assert(a_density[lvl]->nComp()  == 1);
    CH_assert(a_velocity[lvl]->nComp() == SpaceDim);

    // Call level versions
    compute_flux(*a_flux[lvl], *a_density[lvl], *a_velocity[lvl], lvl);
  }
}

void cdr_plasma_stepper::compute_flux(LevelData<EBCellFAB>&       a_flux,
				      const LevelData<EBCellFAB>& a_density,
				      const LevelData<EBCellFAB>& a_velocity,
				      const int                   a_lvl){
  CH_TIME("cdr_plasma_stepper::compute_flux(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_flux(amr)" << endl;
  }

  data_ops::set_value(a_flux, 0.0);
  data_ops::incr(a_flux, a_velocity, 1.0);
  data_ops::multiply_scalar(a_flux, a_density);
}

void cdr_plasma_stepper::compute_J(EBAMRCellData& a_J) const{
  CH_TIME("cdr_plasma_stepper::compute_J(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_J(amr)" << endl;
  }
  
  const int finest_level = m_amr->get_finest_level();

  // Call level versions
  for (int lvl = 0; lvl <= finest_level; lvl++){
    compute_J(*a_J[lvl], lvl);
  }

  m_amr->average_down(a_J, m_realm, m_cdr->get_phase());
  m_amr->interp_ghost(a_J, m_realm, m_cdr->get_phase());
}

void cdr_plasma_stepper::compute_J(LevelData<EBCellFAB>& a_J, const int a_lvl) const{
  CH_TIME("cdr_plasma_stepper::compute_J(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_J(level)" << endl;
  }

  data_ops::set_value(a_J, 0.0);
  
  const int density_comp = 0;

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[a_lvl];
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_species>& spec      = solver_it.get_species();

    if(solver->is_mobile()){
      const int q                       = spec->get_charge();
      const EBAMRCellData& density      = solver->get_state();
      const EBAMRCellData& velo         = solver->get_velo_cell();

      if(q != 0){
	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  const Box& box         = dbl.get(dit());
	  const EBISBox& ebisbox = ebisl[dit()];
	  const EBGraph& ebgraph = ebisbox.getEBGraph();
	  const IntVectSet ivs(box);
      
	  EBCellFAB& J       = a_J[dit()];
	  const EBCellFAB& n = (*density[a_lvl])[dit()];
	  const EBCellFAB& v = (*velo[a_lvl])[dit()];
      
	  EBCellFAB cdr_j(ebisbox, box, SpaceDim);
	  cdr_j.setVal(0.0);
	  for (int comp = 0; comp < SpaceDim; comp++){
	    cdr_j.plus(n, 0, comp, 1);
	  }
	  cdr_j *= v;
	  cdr_j *= q;
      
	  J += cdr_j;

	  // Should we monkey with irregular cells??? 
	}
      }
    }
  }
  
  data_ops::scale(a_J, units::s_Qe);
}

void cdr_plasma_stepper::compute_rho(){
  CH_TIME("cdr_plasma_stepper::compute_rho()");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_rho()" << endl;
  }

  Vector<EBAMRCellData*> densities;
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver> solver = solver_it();
    densities.push_back(&(solver->get_state()));
  }

  this->compute_rho(m_poisson->get_source(), densities, centering::cell_center);
}

void cdr_plasma_stepper::compute_rho(EBAMRCellData& a_rho, const phase::which_phase a_phase){
  CH_TIME("cdr_plasma_stepper::compute_rho(ebamrcelldata, phase)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_rho(ebamrcelldata, phase)" << endl;
  }

  CH_assert(a_phase == m_cdr->get_phase());

  data_ops::set_value(a_rho, 0.0);

  Vector<EBAMRCellData*> densities = m_cdr->get_states(); // Get densities from solver
  
  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    // Add volumetric charge 
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const EBAMRCellData& density       = *(densities[solver_it.index()]);
      const RefCountedPtr<cdr_species>& spec = solver_it.get_species();

      if(spec->get_charge() != 0){
	data_ops::incr(*a_rho[lvl], *density[lvl], spec->get_charge());
      }
    }

    // Scale by s_Qe/s_eps0
    data_ops::scale(*a_rho[lvl], units::s_Qe);
  }

  m_amr->average_down(a_rho, m_realm, a_phase);
  m_amr->interp_ghost(a_rho, m_realm, a_phase);
}

void cdr_plasma_stepper::compute_rho(MFAMRCellData&                 a_rho,
				     const Vector<EBAMRCellData*>&  a_densities,
				     const centering                a_centering){
  CH_TIME("cdr_plasma_stepper::compute_rho(mfamrcell, vec(ebamrcell))");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_rho(mfamrcell, vec(ebamrcell))" << endl;
  }

  data_ops::set_value(a_rho, 0.0);

  EBAMRCellData rho_gas;
  m_amr->allocate_ptr(rho_gas); 
  m_amr->alias(rho_gas, phase::gas, a_rho); 
  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    // Add volumetric charge 
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const EBAMRCellData& density       = *(a_densities[solver_it.index()]);
      const RefCountedPtr<cdr_species>& spec = solver_it.get_species();

      if(spec->get_charge() != 0){
	data_ops::incr(*rho_gas[lvl], *density[lvl], spec->get_charge());
      }
    }

    // Scale by s_Qe
    data_ops::scale(*a_rho[lvl], units::s_Qe);
  }

  m_amr->average_down(a_rho, m_realm);
  m_amr->interp_ghost(a_rho, m_realm);

  // Transform to centroids
  if(a_centering == centering::cell_center){
    m_amr->interpolate_to_centroids(rho_gas, m_realm, phase::gas);
  }
}

void cdr_plasma_stepper::deallocate_solver_internals(){
  CH_TIME("cdr_plasma_stepper::deallocate_solver_internals");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::deallocate_solver_internals" << endl;
  }

  m_cdr->deallocate_internals();
  m_rte->deallocate_internals();
  m_poisson->deallocate_internals();
  m_sigma->deallocate_internals();
}

void cdr_plasma_stepper::extrapolate_to_eb(Vector<EBAMRIVData*>&         a_extrap,
					   const phase::which_phase      a_phase,
					   const Vector<EBAMRCellData*>& a_data){
  CH_TIME("cdr_plasma_stepper::extrapolate_to_eb(vec)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_to_eb(vec)" << endl;
  }

  CH_assert(a_extrap.size() == a_data.size());

  for (int i = 0; i < a_extrap.size(); i++){
    this->extrapolate_to_eb(*a_extrap[i], a_phase, *a_data[i]);
  }
}

void cdr_plasma_stepper::extrapolate_to_eb(EBAMRIVData& a_extrap, const phase::which_phase a_phase, const EBAMRCellData& a_data){
  CH_TIME("cdr_plasma_stepper::extrapolate_to_eb");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_to_eb" << endl;
  }

  const irreg_amr_stencil<eb_centroid_interp>& stencils = m_amr->get_eb_centroid_interp_stencils(m_realm, a_phase);
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    extrapolate_to_eb(*a_extrap[lvl], a_phase, *a_data[lvl], lvl);
  }
}

void cdr_plasma_stepper::extrapolate_to_eb(LevelData<BaseIVFAB<Real> >& a_extrap,
					   const phase::which_phase     a_phase,
					   const LevelData<EBCellFAB>&  a_data,
					   const int                    a_lvl){
  CH_TIME("cdr_plasma_stepper::extrapolate_to_eb(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_to_eb(level)" << endl;
  }

  const irreg_amr_stencil<eb_centroid_interp>& stencils = m_amr->get_eb_centroid_interp_stencils(m_realm, a_phase);
  stencils.apply(a_extrap, a_data, a_lvl);
}

void cdr_plasma_stepper::extrapolate_to_domain_faces(EBAMRIFData&             a_extrap,
						     const phase::which_phase a_phase,
						     const EBAMRCellData&     a_data){
  CH_TIME("cdr_plasma_stepper::extrapolate_to_domain_faces(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_to_domain_faces(amr)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    extrapolate_to_domain_faces(*a_extrap[lvl], a_phase, *a_data[lvl], lvl);
  }
}

void cdr_plasma_stepper::extrapolate_to_domain_faces(LevelData<DomainFluxIFFAB>& a_extrap,
						     const phase::which_phase    a_phase,
						     const LevelData<EBCellFAB>& a_data,
						     const int                   a_lvl){
  CH_TIME("cdr_plasma_stepper::extrapolate_to_domain_faces(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_to_domain_faces(level)" << endl;
  }

  const int ncomp = a_data.nComp();
      
  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, a_phase)[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBCellFAB& data         = a_data[dit()];
    const EBISBox& ebisbox        = ebisl[dit()];
    const BaseFab<Real>& data_fab = data.getSingleValuedFAB();
      
    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){
	BaseIFFAB<Real>& extrap = a_extrap[dit()](dir, sit());

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
	}
      }
    }
  }
}

void cdr_plasma_stepper::extrapolate_to_domain_faces(Vector<EBAMRIFData*>&         a_extrap,
						     const phase::which_phase      a_phase,
						     const Vector<EBAMRCellData*>& a_data){
  CH_TIME("cdr_plasma_stepper::extrapolate_to_domain_faces(vec)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_to_domain_faces(vec)" << endl;
  }

  CH_assert(a_extrap.size() == a_data.size());

  for (int i = 0; i < a_extrap.size(); i++){
    this->extrapolate_to_domain_faces(*a_extrap[i], a_phase, *a_data[i]);
  }
}

void cdr_plasma_stepper::get_cdr_max(Real& a_cdr_max, std::string& a_solver_name){
  CH_TIME("cdr_plasma_stepper::get_cdr_max");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::get_cdr_max" << endl;
  }

  const int comp = 0;
  
  a_cdr_max = -1.E99;
  a_solver_name = "invalid solver";
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    Real max, min;
    data_ops::get_max_min(max, min, solver->get_state(), comp);

    if(max > a_cdr_max){
      a_cdr_max     = max;
      a_solver_name = solver->get_name();
    }
  }
}

void cdr_plasma_stepper::set_cdr(RefCountedPtr<cdr_layout<cdr_solver>>& a_cdr){
  CH_TIME("cdr_plasma_stepper::set_cdr");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_cdr" << endl;
  }
  m_cdr = a_cdr;
}

void cdr_plasma_stepper::set_poisson(RefCountedPtr<field_solver>& a_poisson){
  CH_TIME("cdr_plasma_stepper::set_poisson");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_poisson" << endl;
  }
  m_poisson = a_poisson;
}

void cdr_plasma_stepper::set_rte(RefCountedPtr<rte_layout<rte_solver>>& a_rte){
  CH_TIME("cdr_plasma_stepper::set_rte");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_rte" << endl;
  }
  m_rte = a_rte;
}

void cdr_plasma_stepper::setup_solvers(){
  CH_TIME("cdr_plasma_stepper::setup_solvers");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::setup_solvers" << endl;
  }
  parse_options();
  this->sanity_check();

  // Make solvers
  this->setup_cdr();
  this->setup_rte(); 
  this->setup_poisson();
  this->setup_sigma();

  this->set_solver_verbosity();

  // Allocate internal memory
  this->allocate_internals();
}

void cdr_plasma_stepper::initial_data(){
  CH_TIME("cdr_plasma_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::initial_data" << endl;
  }


  m_cdr->initial_data();        // Initial data comes in through cdr_species, in this case supplied by physics
  if(!m_rte->is_stationary()){
    m_rte->initial_data();
  }
  this->initial_sigma();

  // Solve Poisson equation
  this->solve_poisson();

  // Fill solvers with velocity and diffusion
  this->compute_cdr_velocities();
  this->compute_cdr_diffusion();

  // Do stationary RTE solve if we must
  if(this->stationary_rte()){                  // Solve RTE equations by using initial data and electric field
    const Real dummy_dt = 1.0;

    this->solve_rte(dummy_dt);                 // Argument does not matter, it's a stationary solver.
  }
}

void cdr_plasma_stepper::initial_sigma(){
  CH_TIME("cdr_plasma_stepper::initial_sigma");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::initial_sigma" << endl;
  }

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  EBAMRIVData& sigma = m_sigma->get_state();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, phase::gas)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& state = (*sigma[lvl])[dit()];

      const EBISBox& ebisbox = ebisl[dit()];
      const IntVectSet& ivs  = state.getIVS();
      const EBGraph& ebgraph = state.getEBGraph();
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = origin + vof.gridIndex()*dx + ebisbox.bndryCentroid(vof)*dx;
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_physics->initial_sigma(m_time, pos);
	}
      }
    }
  }

  m_amr->average_down(sigma, m_realm, phase::gas);
  m_sigma->reset_cells(sigma);
}

void cdr_plasma_stepper::project_flux(LevelData<BaseIVFAB<Real> >&       a_projected_flux,
				      const LevelData<BaseIVFAB<Real> >& a_flux,
				      const int                          a_lvl){
  CH_TIME("cdr_plasma_stepper::project_flux(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::project_flux(level)" << endl;
  }

  CH_assert(a_projected_flux.nComp() == 1);
  CH_assert(a_flux.nComp()           == SpaceDim);

  const DisjointBoxLayout& dbl  = m_amr->get_grids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box         = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    
    BaseIVFAB<Real>& proj_flux  = a_projected_flux[dit()];
    const BaseIVFAB<Real>& flux = a_flux[dit()];

    VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof     = vofit();
      const RealVect& normal  = ebisbox.normal(vof);
      const RealVect vec_flux = RealVect(D_DECL(flux(vof,0), flux(vof,1), flux(vof,2)));
      
      // For EB's, the geometrical normal vector is opposite of the finite volume method normal
      proj_flux(vof,0) = PolyGeom::dot(vec_flux, -normal);
    }
  }
}

void cdr_plasma_stepper::project_flux(EBAMRIVData& a_projected_flux, const EBAMRIVData& a_flux){
  CH_TIME("cdr_plasma_stepper::project_flux");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::project_flux" << endl;
  }

  const int comp         = 0;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    cdr_plasma_stepper::project_flux(*a_projected_flux[lvl], *a_flux[lvl], lvl);
  }
}

void cdr_plasma_stepper::project_domain(EBAMRIFData& a_projected_flux, const EBAMRIFData& a_flux){
  CH_TIME("cdr_plasma_stepper::project_domain");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::project_domain" << endl;
  }

  const int comp         = 0;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    CH_assert(a_projected_flux[lvl]->nComp() == 1);
    CH_assert(a_flux[lvl]->nComp()           == SpaceDim);

    const DisjointBoxLayout& dbl  = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl       = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

      const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
      for (int dir = 0; dir < SpaceDim; dir++){
	for (SideIterator sit; sit.ok(); ++sit){
	  BaseIFFAB<Real>& normal_comp    = (*a_projected_flux[lvl])[dit()](dir, sit());
	  const BaseIFFAB<Real>& gradient = (*a_flux[lvl])[dit()](dir, sit());

	  const int sgn          = sign(sit());
	  const EBGraph& ebgraph = normal_comp.getEBGraph();
	  const IntVectSet& ivs  = normal_comp.getIVS();

	  for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	    const FaceIndex& face = faceit();
	    normal_comp(face, comp) = gradient(face, dir);
	  }
	}
      }
    }
  }
}

void cdr_plasma_stepper::regrid(const int a_lmin, const int a_old_finest, const int a_new_finest){
  CH_TIME("cdr_plasma_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::regrid" << endl;
  }

  this->allocate_internals(); // Allocate memory for time stepper
  this->regrid_solvers(a_lmin, a_old_finest, a_new_finest);
  this->regrid_internals(a_lmin, a_old_finest, a_new_finest);

  // Solvers have been regridded. Now resolve the Poisson equation with the new data
  bool converged = this->solve_poisson();

  // If we don't converge, try new Poisson solver settings
  if(!converged){ 
    if(m_verbosity > 0){
      pout() << "driver::regrid - Poisson solver failed to converge. Trying to auto-tune new settings." << endl;
    }
	  
    m_poisson->auto_tune();
    converged = this->solve_poisson();

    if(!converged){
      if(m_verbosity > 0){
	pout() << "cdr_plasma_stepper::post_regrid - Poisson solver fails to converge" << endl;
      }
    }
  }

  // Compute stuff that is important for the CDR solvers. 
  this->compute_cdr_velocities();
  this->compute_cdr_diffusion();

  // If we're doing a stationary RTE solve, recompute source terms
  if(this->stationary_rte()){     // Solve RTE equations by using data that exists inside solvers
    const Real dummy_dt = 1.0;

    // Need new source terms for RTE equations
    this->advance_reaction_network(m_time, dummy_dt);
    this->solve_rte(dummy_dt);    // Argument does not matter, it's a stationary solver.
  }
}

void cdr_plasma_stepper::regrid_solvers(const int a_lmin, const int a_old_finest, const int a_new_finest){
  CH_TIME("cdr_plasma_stepper::regrid_solvers");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::regrid_solvers" << endl;
  }

  m_cdr->regrid(a_lmin,     a_old_finest, a_new_finest);
  m_poisson->regrid(a_lmin, a_old_finest, a_new_finest);
  m_rte->regrid(a_lmin,     a_old_finest, a_new_finest);
  m_sigma->regrid(a_lmin,   a_old_finest, a_new_finest);
}

void cdr_plasma_stepper::reset_dielectric_cells(EBAMRIVData& a_data){
  CH_TIME("cdr_plasma_stepper::reset_dielectric_cells");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::reset_dielectric_cells" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    const MFLevelGrid& mflg      = *m_amr->get_mflg(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      BaseIVFAB<Real>& data  = (*a_data[lvl])[dit()];
      const IntVectSet ivs   = data.getIVS() & mflg.interface_region(box, dit());
      const EBGraph& ebgraph = data.getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	for (int comp = 0; comp < data.nComp(); comp++){
	  data(vof, comp) = 0.0;
	}
      }
    }
  }
}

void cdr_plasma_stepper::sanity_check(){
  CH_TIME("cdr_plasma_stepper::sanity_check");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_physics.isNull());
}

void cdr_plasma_stepper::set_cdr_plasma_physics(const RefCountedPtr<cdr_plasma_physics>& a_physics){
  CH_TIME("cdr_plasma_stepper::set_cdr_plasma_physics");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_cdr_plasma_physics" << endl;
  }

  m_physics = a_physics;
}

void cdr_plasma_stepper::set_potential(Real (*a_potential)(const Real a_time)){
  CH_TIME("cdr_plasma_stepper::set_potential");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_potential" << endl;
  }

  m_potential = a_potential;
}

void cdr_plasma_stepper::parse_verbosity(){
  CH_TIME("cdr_plasma_stepper::parse_verbosity");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::parse_verbosity" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());
  pp.get("verbosity", m_verbosity);
}

void cdr_plasma_stepper::parse_solver_verbosity(){
  CH_TIME("cdr_plasma_stepper::parse_solver_verbosity");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::parse_solver_verbosity" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());
  pp.get("solver_verbosity", m_solver_verbosity);
}

void cdr_plasma_stepper::parse_cfl(){

  ParmParse pp(m_class_name.c_str());
  pp.get("cfl", m_cfl);
  if(m_cfl < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_cfl - CFL cannot be negative!");
  }
}

void cdr_plasma_stepper::parse_relax_time(){
  ParmParse pp(m_class_name.c_str());
  pp.get("relax_time", m_relax_time);
  if(m_relax_time < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_relax_time - relaxation time cannot be negative");
  }
}

void cdr_plasma_stepper::parse_source_growth(){
  ParmParse pp(m_class_name.c_str());
  pp.query("source_growth", m_src_growth);
  if(m_src_growth < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_source_growth - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_source_tolerance(){
  ParmParse pp(m_class_name.c_str());
  pp.query("source_tolerance", m_src_tolerance);
  if(m_src_tolerance < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_source_tolerance - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_min_dt(){
  ParmParse pp(m_class_name.c_str());
  pp.get("min_dt", m_min_dt);
  if(m_min_dt < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_min_dt - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_max_dt(){
  ParmParse pp(m_class_name.c_str());
  pp.get("max_dt", m_max_dt);
  if(m_max_dt < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_max_dt - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_fast_rte(){
  ParmParse pp(m_class_name.c_str());
  pp.get("fast_rte", m_fast_rte);
  if(m_fast_rte <= 0){
    MayDay::Abort("cdr_plasma_stepper::parse_fast_rte - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_fast_poisson(){
  ParmParse pp(m_class_name.c_str());
  pp.get("fast_poisson", m_fast_poisson);
  if(m_fast_poisson <= 0){
    MayDay::Abort("cdr_plasma_stepper::parse_fast_poisson - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_source_comp(){

  std::string str;
  ParmParse pp(m_class_name.c_str());
  pp.get("source_comp", str);
  if(str == "interp"){
    m_interp_sources = true;
  }
  else if(str == "cell_ave"){
    m_interp_sources = false;
  }
  else{
    MayDay::Abort("cdr_plasma_stepper::parse_source_comp - unknown type requested");
  }
}

void cdr_plasma_stepper::set_solver_verbosity(){
  CH_TIME("cdr_plasma_stepper::set_solver_verbosity");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_solver_verbosity" << endl;
  }

  if(!m_cdr.isNull()){
    m_cdr->set_verbosity(m_solver_verbosity);
  }
  if(!m_poisson.isNull()){
    m_poisson->set_verbosity(m_solver_verbosity);
  }
  if(!m_rte.isNull()){
    m_rte->set_verbosity(m_solver_verbosity);
  }
  if(!m_sigma.isNull()){
    m_sigma->set_verbosity(m_solver_verbosity);
  }
}

void cdr_plasma_stepper::set_fast_rte(const int a_fast_rte){
  m_fast_rte = a_fast_rte;

  ParmParse pp("cdr_plasma_stepper");
  pp.query("fast_rte", m_fast_rte);
  if(m_fast_rte <= 0){
    m_fast_rte = a_fast_rte;
  }
}

void cdr_plasma_stepper::set_fast_poisson(const int a_fast_poisson){
  m_fast_poisson = a_fast_poisson;

  ParmParse pp("cdr_plasma_stepper");
  pp.query("fast_poisson", m_fast_poisson);
  if(m_fast_poisson <= 0){
    m_fast_poisson = a_fast_poisson;
  }
}

void cdr_plasma_stepper::set_source_computation(){
  m_interp_sources = true;

  std::string str;
  ParmParse pp("cdr_plasma_stepper");
  if(pp.contains("source_comp")){
    pp.get("source_comp", str);
    if(str == "interp"){
      m_interp_sources = true;
    }
    else if(str == "cell_ave"){
      m_interp_sources = false;
    }
    else{
      MayDay::Abort("cdr_plasma_stepper::set_source_computation - unknown type requested");
    }
  }
}

void cdr_plasma_stepper::set_min_dt(const Real a_min_dt){
  const Real zero = 0.0;
  m_min_dt = Max(a_min_dt, zero);

  ParmParse pp("cdr_plasma_stepper");
  pp.query("min_dt", m_min_dt);
  if(m_min_dt < 0.0){
    m_min_dt = a_min_dt;
  }
}

void cdr_plasma_stepper::set_max_dt(const Real a_max_dt){
  const Real zero = 0.0;
  m_max_dt = Max(a_max_dt, zero);

  ParmParse pp("cdr_plasma_stepper");
  pp.query("max_dt", m_max_dt);
  if(m_max_dt < 0.0){
    m_max_dt = a_max_dt;
  }
}

void cdr_plasma_stepper::set_cfl(const Real a_cfl){
  m_cfl = a_cfl;

  ParmParse pp("cdr_plasma_stepper");
  pp.query("cfl", m_cfl);
  if(m_cfl < 0.0){
    m_cfl = a_cfl;
  }
}

void cdr_plasma_stepper::set_relax_time(const Real a_relax_time){
  m_relax_time = a_relax_time;

  ParmParse pp("cdr_plasma_stepper");
  pp.query("relax_time", m_relax_time);
  if(m_relax_time < 0.0){
    m_relax_time = a_relax_time;
  }
}

void cdr_plasma_stepper::set_source_growth(const Real a_src_growth){
  m_src_growth = a_src_growth;

  ParmParse pp("cdr_plasma_stepper");
  pp.query("source_growth", m_src_growth);
  if(m_src_growth < 0.0){
    m_src_growth = a_src_growth;
  }
}

void cdr_plasma_stepper::set_source_growth_tolerance(const Real a_src_tolerance){
  m_src_tolerance = a_src_tolerance;

  ParmParse pp("cdr_plasma_stepper");
  pp.query("source_tolerance", m_src_tolerance);
  if(m_src_tolerance < 0.0){
    m_src_tolerance = a_src_tolerance;
  }
}

void cdr_plasma_stepper::setup_cdr(){
  CH_TIME("cdr_plasma_stepper::setup_cdr");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::setup_cdr" << endl;
  }

  m_cdr->set_verbosity(m_solver_verbosity);
  m_cdr->parse_options();
  m_cdr->set_amr(m_amr);
  m_cdr->set_computational_geometry(m_compgeom);
  m_cdr->set_phase(phase::gas);
  m_cdr->sanity_check();
  m_cdr->set_realm(m_realm);
}

void cdr_plasma_stepper::allocate() {
  m_cdr->allocate_internals();
  m_poisson->allocate_internals();
  m_rte->allocate_internals();
  m_sigma->allocate_internals();
}

void cdr_plasma_stepper::setup_poisson(){
  CH_TIME("cdr_plasma_stepper::setup_poisson");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::setup_poisson" << endl;
  }

  m_poisson->set_verbosity(m_solver_verbosity);
  m_poisson->parse_options();
  m_poisson->set_amr(m_amr);
  m_poisson->set_computational_geometry(m_compgeom);
  m_poisson->set_realm(m_realm);

  m_poisson->set_poisson_wall_func(0, Side::Lo, m_wall_func_x_lo); // Set function-based Poisson on xlo
  m_poisson->set_poisson_wall_func(0, Side::Hi, m_wall_func_x_hi); // Set function-based Poisson on xhi
  m_poisson->set_poisson_wall_func(1, Side::Lo, m_wall_func_y_lo); // Set function-based Poisson on ylo
  m_poisson->set_poisson_wall_func(1, Side::Hi, m_wall_func_y_hi); // Set function-based Poisson on yhi
#if CH_SPACEDIM==3
  m_poisson->set_poisson_wall_func(2, Side::Lo, m_wall_func_z_lo); // Set function-based Poisson on zlo
  m_poisson->set_poisson_wall_func(2, Side::Hi, m_wall_func_z_hi); // Set function-based Poisson on zhi
#endif
  m_poisson->set_potential(m_potential); // Needs to happen AFTER set_poisson_wall_func

  m_poisson->sanity_check();
  //  m_poisson->allocate_internals();
}

void cdr_plasma_stepper::setup_rte(){
  CH_TIME("cdr_plasma_stepper::setup_rte");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::setup_rte" << endl;
  }

  m_rte->set_verbosity(m_solver_verbosity);
  m_rte->parse_options();
  m_rte->set_phase(phase::gas);
  m_rte->set_amr(m_amr);
  m_rte->set_computational_geometry(m_compgeom);
  m_rte->sanity_check();
  m_rte->set_realm(m_realm);
  //  m_rte->allocate_internals();
}

void cdr_plasma_stepper::setup_sigma(){
  CH_TIME("cdr_plasma_stepper::setup_sigma");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::setup_sigma" << endl;
  }

  m_sigma = RefCountedPtr<sigma_solver> (new sigma_solver());
  m_sigma->set_amr(m_amr);
  m_sigma->set_verbosity(m_solver_verbosity);
  m_sigma->set_computational_geometry(m_compgeom);
  m_sigma->set_realm(m_realm);
  //  m_sigma->allocate_internals();
}

void cdr_plasma_stepper::solver_dump(){
  CH_TIME("cdr_plasma_stepper::solver_dump");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::solver_dump" << endl;
  }

  m_cdr->write_plot_file();
  m_poisson->write_plot_file();
  m_rte->write_plot_file();
}

void cdr_plasma_stepper::solve_rte(const Real a_dt){
  CH_TIME("cdr_plasma_stepper::solve_rte()");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::solve_rte()" << endl;
  }

  const phase::which_phase rte_phase = m_rte->get_phase();

  EBAMRCellData E;
  m_amr->allocate(E, m_realm, rte_phase, SpaceDim);
  this->compute_E(E, rte_phase, m_poisson->get_state());

  Vector<EBAMRCellData*> states     = m_rte->get_states();
  Vector<EBAMRCellData*> rhs        = m_rte->get_sources();
  Vector<EBAMRCellData*> cdr_states = m_cdr->get_states();

  this->solve_rte(states, rhs, cdr_states, E, m_time, a_dt, centering::cell_center);
}

void cdr_plasma_stepper::solve_rte(Vector<EBAMRCellData*>&       a_rte_states,
				   Vector<EBAMRCellData*>&       a_rte_sources,
				   const Vector<EBAMRCellData*>& a_cdr_states,
				   const EBAMRCellData&          a_E,
				   const Real                    a_time,
				   const Real                    a_dt,
				   const centering               a_centering){
  CH_TIME("cdr_plasma_stepper::solve_rte(full)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::solve_rte(full)" << endl;
  }

  //  this->compute_rte_sources(a_rte_sources, a_cdr_states, a_E, a_time, a_centering);
  //  advance_reaction_network(a_time, a_dt);

  for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    
    RefCountedPtr<rte_solver>& solver = solver_it();
    EBAMRCellData& state              = *a_rte_states[idx];
    EBAMRCellData& rhs                = *a_rte_sources[idx];
    solver->advance(a_dt, state, rhs);
  }
}

void cdr_plasma_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("cdr_plasma_stepper::synchronize_solver_times");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::synchronize_solver_times" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;

  m_cdr->set_time(a_step,     a_time, a_dt);
  m_poisson->set_time(a_step, a_time, a_dt);
  m_rte->set_time(a_step,     a_time, a_dt);
  m_sigma->set_time(a_step,   a_time, a_dt);
}

Real cdr_plasma_stepper::compute_electrode_current(){
  CH_TIME("cdr_plasma_stepper::compute_electrode_current");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_electrode_current" << endl;
  }

  // Need to copy onto temporary storage because 
  EBAMRIVData charge_flux;
  m_amr->allocate(charge_flux, m_realm, m_cdr->get_phase(), 1);
  data_ops::set_value(charge_flux, 0.0);
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<cdr_species>& spec      = solver_it.get_species();
    const EBAMRIVData& solver_flux          = solver->get_ebflux();

    data_ops::incr(charge_flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  this->reset_dielectric_cells(charge_flux);
  m_amr->conservative_average(charge_flux, m_realm, m_cdr->get_phase());

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->get_dx()[compute_lvl];
  for (DataIterator dit = m_amr->get_grids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
    const BaseIVFAB<Real>& flx = (*charge_flux[compute_lvl])[dit()];

    const IntVectSet ivs = flx.getIVS() & m_amr->get_grids(m_realm)[compute_lvl].get(dit());
    for (VoFIterator vofit(ivs, flx.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& vof   = vofit();
      const Real& bndryFrac = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[compute_lvl][dit()].bndryArea(vof);
      const Real& flux = flx(vof, 0);
      sum += flux*bndryFrac;
    }
  }

  sum *= pow(dx, SpaceDim-1);


#ifdef CH_MPI
  const Real sum1 = sum;
  sum = EBLevelDataOps::parallelSum(sum1);  
#endif

  return sum;
}

Real cdr_plasma_stepper::compute_dielectric_current(){
  CH_TIME("cdr_plasma_stepper::compute_dielectric_current");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_dielectric_current" << endl;
  }

  // Need to copy onto temporary storage because 
  EBAMRIVData charge_flux;
  m_amr->allocate(charge_flux, m_realm, m_cdr->get_phase(), 1);
  data_ops::set_value(charge_flux, 0.0);
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<cdr_species>& spec      = solver_it.get_species();
    const EBAMRIVData& solver_flux          = solver->get_ebflux();

    data_ops::incr(charge_flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  m_sigma->reset_cells(charge_flux);
  m_amr->conservative_average(charge_flux, m_realm, m_cdr->get_phase());

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->get_dx()[compute_lvl];
  for (DataIterator dit = m_amr->get_grids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
    const BaseIVFAB<Real>& flx = (*charge_flux[compute_lvl])[dit()];

    const IntVectSet ivs = flx.getIVS() & m_amr->get_grids(m_realm)[compute_lvl].get(dit());
    for (VoFIterator vofit(ivs, flx.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& vof   = vofit();
      const Real& bndryFrac = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[compute_lvl][dit()].bndryArea(vof);
      const Real& flux = flx(vof, 0);
      sum += flux*bndryFrac;
    }
  }

  sum *= pow(dx, SpaceDim-1);


#ifdef CH_MPI
  const Real sum1 = sum;
  sum = EBLevelDataOps::parallelSum(sum1);  
#endif

  return sum;
}

Real cdr_plasma_stepper::compute_domain_current(){
  CH_TIME("cdr_plasma_stepper::compute_domain_current");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_domain_current" << endl;
  }

  const int comp = 0;

  // Need to copy onto temporary storage because 
  EBAMRIFData charge_flux;
  m_amr->allocate(charge_flux, m_realm, m_cdr->get_phase(), 1);
  data_ops::set_value(charge_flux, 0.0);
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<cdr_species>& spec      = solver_it.get_species();
    const EBAMRIFData& solver_flux          = solver->get_domainflux();

    data_ops::incr(charge_flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->get_dx()[compute_lvl];
  for (DataIterator dit = m_amr->get_grids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
    const DomainFluxIFFAB& flux = (*charge_flux[compute_lvl])[dit()];

    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){
	const BaseIFFAB<Real>& fluxdir = flux(dir, sit());

	FaceStop::WhichFaces stopcrit = FaceStop::AllBoundaryOnly;
	const IntVectSet& ivs  = fluxdir.getIVS();
	const EBGraph& ebgraph = fluxdir.getEBGraph();

	for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  sum += sign(sit())*fluxdir(faceit(), comp);
	}
      }
    }
  }

  sum *= pow(dx, SpaceDim-1);


#ifdef CH_MPI
  const Real sum1 = sum;
  sum = EBLevelDataOps::parallelSum(sum1);  
#endif

  return sum;
}

Real cdr_plasma_stepper::compute_ohmic_induction_current(){
  CH_TIME("cdr_plasma_stepper::compute_relaxation_time");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_relaxation_time" << endl;
  }

  Real current = 0.0;

  EBAMRCellData J, E, JdotE;
  m_amr->allocate(J,     m_realm, m_cdr->get_phase(), SpaceDim);
  m_amr->allocate(E,     m_realm, m_cdr->get_phase(), SpaceDim);
  m_amr->allocate(JdotE, m_realm, m_cdr->get_phase(), SpaceDim);

  this->compute_E(E, m_cdr->get_phase(), m_poisson->get_state());
  this->compute_J(J);

  // Compute J.dot.E 
  data_ops::dot_prod(JdotE, J,E);
  m_amr->average_down(JdotE, m_realm, m_cdr->get_phase());

  // Only compue on coarsest level
  const int coar = 0;
  const Real dx  = m_amr->get_dx()[coar];
  data_ops::kappa_sum(current, *JdotE[coar]);
  current *= pow(dx, SpaceDim);

  return current;
}

Real cdr_plasma_stepper::compute_relaxation_time(){
  CH_TIME("cdr_plasma_stepper::compute_relaxation_time");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_relaxation_time" << endl;
  }

  const int comp         = 0;
  const int finest_level = 0;
  const Real SAFETY      = 1.E-20;

  Real t1 = MPI_Wtime();
  EBAMRCellData E, J, dt;
  m_amr->allocate(E,  m_realm, m_cdr->get_phase(), SpaceDim);
  m_amr->allocate(J,  m_realm, m_cdr->get_phase(), SpaceDim);
  m_amr->allocate(dt, m_realm, m_cdr->get_phase(), 1);

  data_ops::set_value(dt, 1.234567E89);

  this->compute_E(E, m_cdr->get_phase(), m_poisson->get_state());
  this->compute_J(J);

  // Find the largest electric field in each direction
  Vector<Real> max_E(SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++){

    Real max, min;
    data_ops::get_max_min(max, min, E, dir);
    max_E[dir] = Max(Abs(max), Abs(min));
  }

  const int finest_relax_level = finest_level;
  
  for (int lvl = 0; lvl <= finest_relax_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_cdr->get_phase())[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box         = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);
      
      EBCellFAB& dt_fab  = (*dt[lvl])[dit()];
      const EBCellFAB& e = (*E[lvl])[dit()];
      const EBCellFAB& j = (*J[lvl])[dit()];

      EBCellFAB e_magnitude(ebisbox, box, 1);
      EBCellFAB j_magnitude(ebisbox, box, 1);

      // Compute magnitudes, increment with safety factor to avoid division by zero
      e_magnitude.setVal(0.0);
      j_magnitude.setVal(0.0);
      
      data_ops::vector_length(e_magnitude, e, box);
      data_ops::vector_length(j_magnitude, j, box);
      j_magnitude += SAFETY;      

      dt_fab.setVal(units::s_eps0);
      dt_fab *= e_magnitude;
      dt_fab /= j_magnitude;

      // Now do the irregular cells
      VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect ee = RealVect(D_DECL(e(vof, 0), e(vof, 1), e(vof, 2)));
	const RealVect jj = RealVect(D_DECL(j(vof, 0), j(vof, 1), j(vof, 2)));

	dt_fab(vof, comp) = Abs(units::s_eps0*ee.vectorLength()/(1.E-20 + jj.vectorLength()));
      }
    }
  }

  // Find the smallest dt
  Real min_dt = 1.E99;
  Real max, min;
  data_ops::get_max_min(max, min, dt, comp);
  min_dt = Min(min_dt, min);

  return min_dt;
}

Real cdr_plasma_stepper::get_time(){
  return m_time;
}

Real cdr_plasma_stepper::get_dt(){
  return m_dt;
}

Real cdr_plasma_stepper::get_cfl_dt(){
  return m_dt_cfl;
}

RefCountedPtr<cdr_layout<cdr_solver>>& cdr_plasma_stepper::get_cdr(){
  return m_cdr;
}

RefCountedPtr<field_solver>& cdr_plasma_stepper::get_poisson(){
  return m_poisson;
}

RefCountedPtr<rte_layout<rte_solver>>& cdr_plasma_stepper::get_rte(){
  return m_rte;
}

RefCountedPtr<sigma_solver>& cdr_plasma_stepper::get_sigma(){
  return m_sigma;
}

// New functions for driver
void cdr_plasma_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const{
  CH_TIME("driver::write_checkpoint_data");
  if(m_verbosity > 3){
    pout() << "driver::write_checkpoint_data" << endl;
  }

  // CDR solvers checkpoint their data
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->write_checkpoint_level(a_handle, a_lvl);
  }

  // RTE solvers checkpoint their data
  for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<rte_solver>& solver = solver_it();
    solver->write_checkpoint_level(a_handle, a_lvl);
  }

  m_poisson->write_checkpoint_level(a_handle, a_lvl);
  m_sigma->write_checkpoint_level(a_handle, a_lvl);
}

void cdr_plasma_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("driver::read_checkpoint_data");
  if(m_verbosity > 3){
    pout() << "driver::read_checkpoint_data" << endl;
  }

  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->read_checkpoint_level(a_handle, a_lvl);
  }

  for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->read_checkpoint_level(a_handle, a_lvl);
  }

  m_poisson->read_checkpoint_level(a_handle, a_lvl);
  m_sigma->read_checkpoint_level(a_handle, a_lvl);
}

int cdr_plasma_stepper::get_num_plot_vars() const{
  CH_TIME("cdr_plasma_stepper::get_num_plot_vars");
  if(m_verbosity > 3){
    pout() << "cdr_plasma_stepper::get_num_plot_vars" << endl;
  }
  int ncomp = 0;
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    ncomp += solver->get_num_plotvars();
  }
  
  for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    ncomp += solver->get_num_plotvars();
  }

  ncomp += m_poisson->get_num_plotvars();
  ncomp += m_sigma->get_num_plotvars();
  ncomp += SpaceDim; // For plotting the current density

  return ncomp;
}

void cdr_plasma_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("cdr_plasma_stepper::write_plot_data");
  if(m_verbosity > 3){
    pout() << "cdr_plasma_stepper::write_plot_data" << endl;
  }

  // Poisson solver copies over its output data
  a_plotvar_names.append(m_poisson->get_plotvar_names());
  m_poisson->write_plot_data(a_output, a_icomp);

  // Surface charge solver writes
  a_plotvar_names.append(m_sigma->get_plotvar_names());
  m_sigma->write_plot_data(a_output, a_icomp);

  // CDR solvers copy their output data
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    a_plotvar_names.append(solver->get_plotvar_names());
    solver->write_plot_data(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    a_plotvar_names.append(solver->get_plotvar_names());
    solver->write_plot_data(a_output, a_icomp);
  }

  // Write the current to the output
  this->write_J(a_output, a_icomp);
  a_plotvar_names.push_back("x-J");
  a_plotvar_names.push_back("y-J");
  if(SpaceDim == 3){
    a_plotvar_names.push_back("z-J");
  }
  
}

void cdr_plasma_stepper::write_J(EBAMRCellData& a_output, int& a_icomp) const{
  CH_TIME("cdr_plasma_stepper::write_J");
  if(m_verbosity > 3){
    pout() << "cdr_plasma_stepper::write_J" << endl;
  }

  // Allocates storage and computes J
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, phase::gas, SpaceDim);
  this->compute_J(scratch);

  const Interval src_interv(0, SpaceDim-1);
  const Interval dst_interv(a_icomp, a_icomp + SpaceDim -1);
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }
  a_icomp += SpaceDim;
}

void cdr_plasma_stepper::post_checkpoint_setup(){
  CH_TIME("cdr_plasma_stepper::post_checkpoint_setup");
  if(m_verbosity > 3){
    pout() << "cdr_plasma_stepper::post_checkpoint_setup" << endl;
  }

  this->solve_poisson();       // Solve Poisson equation by 
  if(this->stationary_rte()){  // Solve RTE equations if stationary solvers
    const Real dummy_dt = 0.0;
    this->solve_rte(dummy_dt); // Argument does not matter, it's a stationary solver.
  }
  this->allocate_internals();  // Prepare internal storage for time stepper

  // Fill solvers with important stuff
  this->compute_cdr_velocities();
  this->compute_cdr_diffusion();
}

void cdr_plasma_stepper::set_poisson_wall_func(const int a_dir, const Side::LoHiSide a_side, Real (*a_func)(const RealVect a_pos)){
  CH_TIME("cdr_plasma_stepper::set_poisson_wall_func(dir, side, func)");
  if(m_verbosity > 4){
    pout() << "cdr_plasma_stepper::set_poisson_wall_func(dir, side, func)" << endl;
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

void cdr_plasma_stepper::set_poisson_wall_func(Real (*a_func)(const RealVect a_pos)){
  CH_TIME("cdr_plasma_stepper::set_poisson_wall_func(func)");
  if(m_verbosity > 4){
    pout() << "cdr_plasma_stepper::set_poisson_wall_func(func)" << endl;
  }

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      this->set_poisson_wall_func(dir, sit(), a_func);
    }
  }
}

void cdr_plasma_stepper::print_step_report(){
  CH_TIME("cdr_plasma_stepper::print_step_report");
  if(m_verbosity > 4){
    pout() << "cdr_plasma_stepper::print_step_report" << endl;
  }

  // Compute the maximum electric field
  Real Emax;
  this->compute_Emax(Emax, phase::gas);

  //
  Real nmax;
  std::string solver_max;
  this->get_cdr_max(nmax, solver_max);

  const Real cfl_dt = this->get_cfl_dt();

  std::string str;
  if(m_timecode == time_code::advection){
    str = " (Restricted by advection)";
  }
  else if(m_timecode == time_code::advection_diffusion){
    str = " (Restricted by advection-diffusion)";
  }
  else if(m_timecode == time_code::error){
    str = " (Restricted by error)";
  }
  else if(m_timecode == time_code::diffusion){
    str = " (Restricted by diffusion)";
  }
  else if(m_timecode == time_code::source){
    MayDay::Abort("driver::step_report - shouldn't happen, source term has been taken out of the design");
    str = " (Restricted by source term)";
  }
  else if(m_timecode == time_code::relaxation_time){
    str = " (Restricted by relaxation time)";
  }
  else if(m_timecode == time_code::restricted){
    str = " (Restricted by time stepper)";
  }
  else if(m_timecode == time_code::hardcap){
    str = " (Restricted by a hardcap)";
  }
  pout() << "                                   mode  = " << str << endl
         << "                                   cfl   = " << m_dt/cfl_dt << endl
	 << "                                   Emax  = " << Emax << endl
	 << "                                   n_max = " << nmax << "(" + solver_max + ")" << endl;
}
