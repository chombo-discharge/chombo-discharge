/*!
  @file cdr_plasma_stepper.cpp
  @brief Implementation of cdr_plasma_stepper.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "cdr_plasma_stepper.H"
#include "CD_FieldSolverMultigrid.H"
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>
#include "units.H"

#include <ParmParse.H>
#include <EBLevelDataOps.H>

// C-style VLA methods can have occasional memory issues. Use them at your own peril. 
#define USE_FAST_REACTIONS  1
#define USE_FAST_VELOCITIES 0
#define USE_FAST_DIFFUSION  1

#include "CD_NamespaceHeader.H"
using namespace physics::cdr_plasma;

Real cdr_plasma_stepper::s_constant_one(const RealVect a_pos){
  return 1.0;
}

cdr_plasma_stepper::cdr_plasma_stepper(){
  m_className = "cdr_plasma_stepper";
  m_verbosity  = -1;
  m_solver_verbosity = -1;
  m_phase = phase::gas;
  m_realm = Realm::Primal;
  m_subcycle = false;
}

cdr_plasma_stepper::cdr_plasma_stepper(RefCountedPtr<cdr_plasma_physics>& a_physics) : cdr_plasma_stepper() {
  m_physics = a_physics;
}

void cdr_plasma_stepper::postInitialize() {
}

void cdr_plasma_stepper::registerRealms(){
  m_amr->registerRealm(m_realm);
}

void cdr_plasma_stepper::postRegrid(){

}

void cdr_plasma_stepper::registerOperators(){
  CH_TIME("cdr_plasma_stepper::registerOperators");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::registerOperators" << endl;
  }
  
  m_cdr->registerOperators();
  m_fieldSolver->registerOperators();
  m_rte->registerOperators();
  m_sigma->registerOperators();
}

cdr_plasma_stepper::~cdr_plasma_stepper(){
}

int cdr_plasma_stepper::queryGhost(){
  return 3;
}

bool cdr_plasma_stepper::stationary_rte(){
  CH_TIME("cdr_plasma_stepper::stationary_rte");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::stationary_rte" << endl;
  }

  return m_rte->isStationary();
}

bool cdr_plasma_stepper::solve_poisson(){
  CH_TIME("cdr_plasma_stepper::solve_poisson");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::solve_poisson" << endl;
  }

  this->compute_rho();
  const bool converged = m_fieldSolver->solve(m_fieldSolver->getPotential(),
					  m_fieldSolver->getRho(),
					  m_sigma->getPhi(),
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

  const bool converged = m_fieldSolver->solve(a_potential, a_rhs, a_sigma, false);

  return converged;
}

void cdr_plasma_stepper::allocateInternals(){
  CH_TIME("cdr_plasma_stepper::allocateInternals"); 
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::allocateInternals" << endl;
  }

  /// Do nothing
}

void cdr_plasma_stepper::deallocateInternals(){
  CH_TIME("cdr_plasma_stepper::deallocateInternals"); 
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::deallocateInternals" << endl;
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
  m_amr->allocate(E, m_realm, m_cdr->getPhase(), SpaceDim);
  this->compute_E(E, m_cdr->getPhase(), m_fieldSolver->getPotential());


  Vector<EBAMRCellData*> particle_sources = m_cdr->getSources();
  Vector<EBAMRCellData*> Photon_sources   = m_rte->getSources();
  Vector<EBAMRCellData*> particle_states  = m_cdr->getPhis();
  Vector<EBAMRCellData*> Photon_states    = m_rte->getPhis();

  // Call the AMR version (without the gradient)
  advance_reaction_network(particle_sources,
			   Photon_sources, 
			   particle_states,
			   Photon_states,
			   E,
			   a_time,
			   a_dt);
}

void cdr_plasma_stepper::advance_reaction_network(Vector<EBAMRCellData*>&       a_particle_sources,
						  Vector<EBAMRCellData*>&       a_Photon_sources,
						  const Vector<EBAMRCellData*>& a_particle_densities,
						  const Vector<EBAMRCellData*>& a_Photon_densities,
						  const EBAMRCellData&          a_E,
						  const Real&                   a_time,
						  const Real&                   a_dt){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network(nograd)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network(nograd)" << endl;
  }

  const int num_species = m_physics->get_num_CdrSpecies();

  Vector<EBAMRCellData*> grad_cdr(num_species); // Holders for grad(cdr)
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    grad_cdr[idx] = new EBAMRCellData();                            // This storage must be deleted
    m_amr->allocate(*grad_cdr[idx], m_realm, m_cdr->getPhase(), SpaceDim);  // Allocate

    m_amr->computeGradient(*grad_cdr[idx], *a_particle_densities[idx], m_realm, phase::gas); // Compute grad()
    m_amr->averageDown(*grad_cdr[idx], m_realm, m_cdr->getPhase());        // Average down
    m_amr->interpGhost(*grad_cdr[idx], m_realm, m_cdr->getPhase());        // Interpolate ghost cells
  }

  this->advance_reaction_network(a_particle_sources,
				 a_Photon_sources,
				 a_particle_densities,
				 grad_cdr,
				 a_Photon_densities,
				 a_E,
				 a_time,
				 a_dt);

  // Delete extra storage - didn't use smart pointers for this...
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->deallocate(*grad_cdr[idx]);
    delete grad_cdr[idx];
  }
}

void cdr_plasma_stepper::advance_reaction_network(Vector<EBAMRCellData*>&       a_particle_sources,
						  Vector<EBAMRCellData*>&       a_Photon_sources,
						  const Vector<EBAMRCellData*>& a_particle_densities,
						  const Vector<EBAMRCellData*>& a_particle_gradients,
						  const Vector<EBAMRCellData*>& a_Photon_densities,
						  const EBAMRCellData&          a_E,
						  const Real&                   a_time,
						  const Real&                   a_dt){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network(amr)" << endl;
  }


  const int num_species = m_physics->get_num_CdrSpecies();
  const int num_Photons = m_physics->get_num_RtSpecies();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    Vector<LevelData<EBCellFAB>* > particle_sources(num_species);
    Vector<LevelData<EBCellFAB>* > particle_densities(num_species);
    Vector<LevelData<EBCellFAB>* > particle_gradients(num_species);
    Vector<LevelData<EBCellFAB>* > Photon_sources(num_Photons);
    Vector<LevelData<EBCellFAB>* > Photon_densities(num_Photons);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      particle_sources[idx]   = (*a_particle_sources[idx])[lvl];
      particle_densities[idx] = (*a_particle_densities[idx])[lvl];
      particle_gradients[idx] = (*a_particle_gradients[idx])[lvl];
    }

    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      Photon_sources[idx]   = (*a_Photon_sources[idx])[lvl];
      Photon_densities[idx] = (*a_Photon_densities[idx])[lvl];
    }

    // Call the level versions
    advance_reaction_network(particle_sources,
			     Photon_sources,
			     particle_densities,
			     particle_gradients,
			     Photon_densities,
			     *a_E[lvl],
			     a_time,
			     a_dt,
			     lvl);
  }

#if 0 // R.M. May/2020: This is not a good place to do this kind of averaging and interpolating. If need it, do it elsewhere.
  // Average down species
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->averageDown(*a_particle_sources[idx], m_realm, m_cdr->getPhase());
    m_amr->interpGhost(*a_particle_sources[idx], m_realm, m_cdr->getPhase()); // This MAY be when if we extrapolate advection
  }

  // Average down Photon solvers
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->averageDown(*a_Photon_sources[idx], m_realm, m_cdr->getPhase());
  }
#endif
}

void cdr_plasma_stepper::advance_reaction_network(Vector<LevelData<EBCellFAB>* >&       a_particle_sources,
						  Vector<LevelData<EBCellFAB>* >&       a_Photon_sources,
						  const Vector<LevelData<EBCellFAB>* >& a_particle_densities,
						  const Vector<LevelData<EBCellFAB>* >& a_particle_gradients,
						  const Vector<LevelData<EBCellFAB> *>& a_Photon_densities,
						  const LevelData<EBCellFAB>&           a_E,
						  const Real&                           a_time,
						  const Real&                           a_dt,
						  const int                             a_lvl){
  CH_TIME("cdr_plasma_stepper::advance_reaction_network(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::advance_reaction_network(level)" << endl;
  }

  const Real zero = 0.0;

  const int num_Photons  = m_physics->get_num_RtSpecies();
  const int num_species  = m_physics->get_num_CdrSpecies();


  // Stencils for extrapolating things to cell centroids
  const IrregAmrStencil<CentroidInterpolationStencil>& interp_stencils = m_amr->getCentroidInterpolationStencils(m_realm, m_cdr->getPhase());

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real dx                = m_amr->getDx()[a_lvl];
  const RealVect origin        = m_amr->getProbLo();
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    Vector<EBCellFAB*> particle_sources(num_species);
    Vector<EBCellFAB*> particle_densities(num_species);
    Vector<EBCellFAB*> particle_gradients(num_species);
    Vector<EBCellFAB*> particle_velocities(num_species, NULL);
    Vector<EBCellFAB*> Photon_sources(num_Photons);
    Vector<EBCellFAB*> Photon_densities(num_Photons);

    
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<CdrSolver>& solver = solver_it();
      const int idx = solver_it.index();
      particle_sources[idx]   = &(*a_particle_sources[idx])[dit()];
      particle_densities[idx] = &(*a_particle_densities[idx])[dit()];
      particle_gradients[idx] = &(*a_particle_gradients[idx])[dit()];

      if(solver->isMobile()){
	const EBAMRCellData& velo = solver->getCellCenteredVelocity();
	particle_velocities[idx] = &((*velo[a_lvl])[dit()]);
      }
    }

    
    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      Photon_sources[idx]   = &(*a_Photon_sources[idx])[dit()];
      Photon_densities[idx] = &(*a_Photon_densities[idx])[dit()];
    }
    
    // This does all cells
#if USE_FAST_REACTIONS
    advance_reaction_network_reg_fast(particle_sources,
				      Photon_sources,
				      particle_densities,
				      particle_gradients,
				      Photon_densities,
				      a_E[dit()],
				      a_time,
				      a_dt,
				      dx,
				      dbl.get(dit()));

#else
    advance_reaction_network_reg(particle_sources,
				 Photon_sources,
				 particle_densities,
				 particle_gradients,
				 Photon_densities,
				 a_E[dit()],
				 a_time,
				 a_dt,
				 dx,
				 dbl.get(dit()));
#endif

    // This overwrites irregular celles
    advance_reaction_network_irreg(particle_sources,
				   Photon_sources,
				   particle_densities,
				   particle_gradients,
				   particle_velocities,
				   Photon_densities,
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
						      Vector<EBCellFAB*>&       a_Photon_sources,
						      const Vector<EBCellFAB*>& a_particle_densities,
						      const Vector<EBCellFAB*>& a_particle_gradients,
						      const Vector<EBCellFAB*>& a_Photon_densities,
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
    

  const int num_species  = m_physics->get_num_CdrSpecies();
  const int num_Photons  = m_physics->get_num_RtSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     particle_densities(num_species);
  Vector<RealVect> particle_gradients(num_species);
  Vector<Real>     Photon_sources(num_Photons);
  Vector<Real>     Photon_densities(num_Photons);

  // Computed source terms onto here
  EBCellFAB part_src(ebisbox, a_E.getRegion(), num_species);
  EBCellFAB phot_src(ebisbox, a_E.getRegion(), num_Photons);

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
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_particle_densities[idx]).getSingleValuedFAB()(iv, 0);
	particle_densities[idx] = Max(zero, phi);
	particle_gradients[idx] = RealVect(D_DECL((*a_particle_gradients[idx]).getSingleValuedFAB()(iv, 0),
						  (*a_particle_gradients[idx]).getSingleValuedFAB()(iv, 1),
						  (*a_particle_gradients[idx]).getSingleValuedFAB()(iv, 2)));
      }
      for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_Photon_densities[idx]).getSingleValuedFAB()(iv, 0);
	Photon_densities[idx] = Max(zero, phi);
      }

      // Compute source terms
      const Real kappa = 1.0; // Kappa for regular cells
      m_physics->advance_reaction_network(particle_sources,
					  Photon_sources,
					  particle_densities,
					  particle_gradients,
					  Photon_densities,
					  E,
					  pos,
					  a_dx,
					  a_dt,
					  a_time,
					  kappa);

      // Put vector into temporary holders
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	part_src.getSingleValuedFAB()(iv,idx) = particle_sources[idx];
      }
    
      // Put vector into temporary holders
      for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	phot_src.getSingleValuedFAB()(iv,idx) = Photon_sources[idx];
      }
    }
  }

  // Copy temporary storage back to solvers
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    (*a_particle_sources[idx]).setVal(0.0);
    (*a_particle_sources[idx]).plus(part_src, idx, 0, 1);

    // Covered cells are bogus. 
    (*a_particle_sources[idx]).setCoveredCellVal(0.0, 0);
  }

  // Copy temporary storage back to solvers
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    (*a_Photon_sources[idx]).setVal(0.0);
    (*a_Photon_sources[idx]).plus(phot_src, idx, 0, 1);

    // Covered cells are bogus. 
    (*a_Photon_sources[idx]).setCoveredCellVal(0.0, 0);
  }
}

void cdr_plasma_stepper::advance_reaction_network_reg_fast(Vector<EBCellFAB*>&       a_particle_sources,
							   Vector<EBCellFAB*>&       a_Photon_sources,
							   const Vector<EBCellFAB*>& a_particle_densities,
							   const Vector<EBCellFAB*>& a_particle_gradients,
							   const Vector<EBCellFAB*>& a_Photon_densities,
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
  advance_reaction_network_reg_fast2D(a_particle_sources, a_Photon_sources, a_particle_densities, a_particle_gradients,
				      a_Photon_densities, a_E, a_time, a_dt, a_dx, a_box);
#elif CH_SPACEDIM==3
  advance_reaction_network_reg_fast3D(a_particle_sources, a_Photon_sources, a_particle_densities, a_particle_gradients,
				      a_Photon_densities, a_E, a_time, a_dt, a_dx, a_box);
#endif
  
}

void cdr_plasma_stepper::advance_reaction_network_reg_fast2D(Vector<EBCellFAB*>&       a_particle_sources,
							     Vector<EBCellFAB*>&       a_Photon_sources,
							     const Vector<EBCellFAB*>& a_particle_densities,
							     const Vector<EBCellFAB*>& a_particle_gradients,
							     const Vector<EBCellFAB*>& a_Photon_densities,
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
    

  const int num_species  = m_physics->get_num_CdrSpecies();
  const int num_Photons  = m_physics->get_num_RtSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into cdr_plasma_physics. 
  RealVect         pos = RealVect::Zero;
  RealVect         E   = RealVect::Zero;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     particle_densities(num_species);
  Vector<RealVect> particle_gradients(num_species);
  Vector<Real>     Photon_sources(num_Photons);
  Vector<Real>     Photon_densities(num_Photons);

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

  // Temps for Photon source terms and densities
  FArrayBox rte_phi(a_box, num_Photons);
  FArrayBox rte_src(a_box, num_Photons);
  rte_phi.setVal(0.0);
  rte_src.setVal(0.0);
  for (int i = 0; i < num_Photons; i++){
    rte_phi.copy(a_Photon_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
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
      for (int idx = 0; idx < num_Photons; ++idx){
	Photon_densities[idx] = Max(0.0, vla_rte_phi[idx][j][i]);
      }

      E   = RealVect(vla_E[0][j][i], vla_E[1][j][i]);
      pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

      m_physics->advance_reaction_network(particle_sources,
					  Photon_sources,
					  particle_densities,
					  particle_gradients,
					  Photon_densities,
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
      for (int idx = 0; idx < num_Photons; ++idx){
	vla_rte_src[idx][j][i] = Photon_sources[idx];
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
  for (int i = 0; i < num_Photons; i++){
    FArrayBox& src = a_Photon_sources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(rte_src, a_box, i, a_box, 0, 1);
  }
#endif
}


void cdr_plasma_stepper::advance_reaction_network_reg_fast3D(Vector<EBCellFAB*>&       a_particle_sources,
							     Vector<EBCellFAB*>&       a_Photon_sources,
							     const Vector<EBCellFAB*>& a_particle_densities,
							     const Vector<EBCellFAB*>& a_particle_gradients,
							     const Vector<EBCellFAB*>& a_Photon_densities,
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

  const int num_species  = m_physics->get_num_CdrSpecies();
  const int num_Photons  = m_physics->get_num_RtSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into cdr_plasma_physics. 
  RealVect         pos = RealVect::Zero;
  RealVect         E   = RealVect::Zero;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     particle_densities(num_species);
  Vector<RealVect> particle_gradients(num_species);
  Vector<Real>     Photon_sources(num_Photons);
  Vector<Real>     Photon_densities(num_Photons);

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

  // Temps for Photon source terms and densities
  FArrayBox rte_phi(a_box, num_Photons);
  FArrayBox rte_src(a_box, num_Photons);
  rte_phi.setVal(0.0);
  rte_src.setVal(0.0);
  for (int i = 0; i < num_Photons; i++){
    rte_phi.copy(a_Photon_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
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
	for (int idx = 0; idx < num_Photons; ++idx){
	  Photon_densities[idx] = Max(0.0, vla_rte_phi[idx][k][j][i]);
	}

	E   = RealVect(vla_E[0][k][j][i], vla_E[1][k][j][i], vla_E[2][k][j][i]);
	pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

	m_physics->advance_reaction_network(particle_sources,
					    Photon_sources,
					    particle_densities,
					    particle_gradients,
					    Photon_densities,
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
	for (int idx = 0; idx < num_Photons; ++idx){
	  vla_rte_src[idx][k][j][i] = Photon_sources[idx];
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
  for (int i = 0; i < num_Photons; i++){
    FArrayBox& src = a_Photon_sources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(rte_src, a_box, i, a_box, 0, 1);
  }

#endif
}



void cdr_plasma_stepper::advance_reaction_network_irreg(Vector<EBCellFAB*>&          a_particle_sources,
							Vector<EBCellFAB*>&          a_Photon_sources,
							const Vector<EBCellFAB*>&    a_particle_densities,
							const Vector<EBCellFAB*>&    a_particle_gradients,
							const Vector<EBCellFAB*>&    a_particle_velocities,
							const Vector<EBCellFAB*>&    a_Photon_densities,
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
					  a_Photon_sources,
					  a_particle_densities,
					  a_particle_gradients,
					  a_particle_velocities,
					  a_Photon_densities,
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
					 a_Photon_sources,
					 a_particle_densities,
					 a_particle_gradients,
					 a_Photon_densities,
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
							       Vector<EBCellFAB*>&          a_Photon_sources,
							       const Vector<EBCellFAB*>&    a_particle_densities,
							       const Vector<EBCellFAB*>&    a_particle_gradients,
							       const Vector<EBCellFAB*>&    a_particle_velocities,
							       const Vector<EBCellFAB*>&    a_Photon_densities,
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
  
  const int num_Photons  = m_physics->get_num_RtSpecies();
  const int num_species  = m_physics->get_num_CdrSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     particle_sources(num_species, 0.);
  Vector<Real>     Photon_sources(num_Photons, 0.);
  Vector<Real>     particle_densities(num_species, 0.);
  Vector<RealVect> particle_gradients(num_species, RealVect::Zero);
  Vector<Real>     Photon_densities(num_Photons, 0.);

  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
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
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<CdrSolver>& solver = solver_it();
      const int idx = solver_it.index();

      bool applyStencil = true;

#if 0
      if(solver->isMobile()){
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
    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();

      Real phi = 0.0;
      for (int i = 0; i < stencil.size(); i++){
	const VolIndex& ivof = stencil.vof(i);
	const Real& iweight  = stencil.weight(i);
	phi += (*a_Photon_densities[idx])(ivof, 0)*iweight;
      }
      Photon_densities[idx] = Max(zero, phi);
    }

    // Compute source terms
    m_physics->advance_reaction_network(particle_sources,
					Photon_sources,
					particle_densities,
					particle_gradients,
					Photon_densities,
					E,
					pos,
					a_dx,
					a_dt,
					a_time,
					kappa);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_particle_sources[idx])(vof, 0) = particle_sources[idx];
    }
    
    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_Photon_sources[idx])(vof, 0) = Photon_sources[idx];
    }
  }
}

void cdr_plasma_stepper::advance_reaction_network_irreg_kappa(Vector<EBCellFAB*>&          a_particle_sources,
							      Vector<EBCellFAB*>&          a_Photon_sources,
							      const Vector<EBCellFAB*>&    a_particle_densities,
							      const Vector<EBCellFAB*>&    a_particle_gradients,
							      const Vector<EBCellFAB*>&    a_Photon_densities,
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
  
  const int num_Photons  = m_physics->get_num_RtSpecies();
  const int num_species  = m_physics->get_num_CdrSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     Photon_sources(num_Photons);
  Vector<Real>     particle_densities(num_species);
  Vector<RealVect> particle_gradients(num_species);
  Vector<Real>     Photon_densities(num_Photons);

  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
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
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
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

    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      Photon_densities[idx] = Max(zero, (*a_Photon_densities[idx])(vof, 0));
    }

    // Compute source terms
    m_physics->advance_reaction_network(particle_sources,
					Photon_sources,
					particle_densities,
					particle_gradients,
					Photon_densities,
					E,
					pos,
					a_dx,
					a_dt,
					a_time,
					kappa);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_particle_sources[idx])(vof, 0) = particle_sources[idx];
    }
    
    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_Photon_sources[idx])(vof, 0) = Photon_sources[idx];
    }
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_face(Vector<EBAMRFluxData*>&       a_diffusionCoefficient_face,
						 const Vector<EBAMRCellData*>& a_cdr_densities,
						 const EBAMRCellData&          a_E,
						 const Real&                   a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_face(full)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_face(full)" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->get_num_CdrSpecies();
  const int finest_level = m_amr->getFinestLevel();

  // Allocate data for cell-centered diffusion coefficients
  Vector<EBAMRCellData> diffco(num_species);
  Vector<Real> cdr_densities(num_species);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();

    m_amr->allocate(diffco[idx], m_realm, m_cdr->getPhase(), ncomp);
  }

  // Call the cell version
  compute_cdr_diffco_cell(diffco, a_cdr_densities, a_E, a_time);

  // Now compute face-centered things by taking the average of cell-centered things
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();

    if(solver->isDiffusive()){ // Only need to do this for diffusive things
      data_ops::set_value(*a_diffusionCoefficient_face[idx], 1.0);
      m_amr->averageDown(diffco[idx], m_realm, m_cdr->getPhase());
      m_amr->interpGhost(diffco[idx], m_realm, m_cdr->getPhase()); 

      data_ops::average_cell_to_face_allcomps(*a_diffusionCoefficient_face[idx], diffco[idx], m_amr->getDomains());
    }
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_cell(Vector<EBAMRCellData>&       a_diffusionCoefficient_cell,
						 const Vector<EBAMRCellData*>& a_cdr_densities,
						 const EBAMRCellData&          a_E,
						 const Real&                   a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_cell(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_cell(amr)" << endl;
  }
  
  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->get_num_CdrSpecies();
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<LevelData<EBCellFAB>* > diffco(num_species);
    Vector<LevelData<EBCellFAB>* > cdr_densities(num_species);
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx] = a_diffusionCoefficient_cell[idx][lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    // Cell the level version
    compute_cdr_diffco_cell(diffco, cdr_densities, *a_E[lvl], lvl, a_time);
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_cell(Vector<LevelData<EBCellFAB>* >&       a_diffusionCoefficient_cell,
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
  const int num_species  = m_physics->get_num_CdrSpecies();

  const IrregAmrStencil<CentroidInterpolationStencil>& interp_stencils = m_amr->getCentroidInterpolationStencils(m_realm, m_cdr->getPhase());

  // Call the level version
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const RealVect origin         = m_amr->getProbLo();
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = IntVectSet(box);
    const EBCellFAB& E     = a_E[dit()];

    Vector<EBCellFAB*> diffco(num_species);
    Vector<EBCellFAB*> cdr_densities(num_species);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx]        = &(*a_diffusionCoefficient_cell[idx])[dit()];
      cdr_densities[idx] = &(*a_cdr_densities[idx])[dit()];
    }

    // Regular cells
#if USE_FAST_DIFFUSION
    compute_cdr_diffco_cell_reg_fast(diffco, cdr_densities, E, dbl.get(dit()), m_amr->getDx()[a_lvl], a_time);
#else
    compute_cdr_diffco_cell_reg(diffco, cdr_densities, E, dbl.get(dit()), m_amr->getDx()[a_lvl], a_time);
#endif

    // Irregular cells
    compute_cdr_diffco_cell_irreg(diffco, cdr_densities, E, dbl.get(dit()), m_amr->getDx()[a_lvl],
				  interp_stencils[a_lvl][dit()], a_time, a_lvl, dit());


  }
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_reg(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
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
  const int num_species = m_physics->get_num_CdrSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

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

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      const Real phi = (*a_cdr_densities[idx]).getSingleValuedFAB()(iv, comp);
      cdr_densities[idx] = Max(0.0, phi);
    }

    const Vector<Real> coeffs = m_physics->compute_cdr_diffusion_coefficients(a_time,
									      pos,
									      E,
									      cdr_densities);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->isDiffusive()){
	tmp.getSingleValuedFAB()(iv,idx) = coeffs[idx];
      }
    }
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    if(solver_it()->isDiffusive()){
      (*a_diffusionCoefficient_cell[idx]).setVal(0.0);
      (*a_diffusionCoefficient_cell[idx]).plus(tmp, idx, 0, 1);

      // Covered cells are bogus. 
      (*a_diffusionCoefficient_cell[idx]).setCoveredCellVal(0.0, 0);
    }
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
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
  compute_cdr_diffco_cell_reg_fast2D(a_diffusionCoefficient_cell, a_cdr_densities, a_E, a_box, a_dx, a_time);
#elif CH_SPACEDIM==3
  compute_cdr_diffco_cell_reg_fast3D(a_diffusionCoefficient_cell, a_cdr_densities, a_E, a_box, a_dx, a_time);
#endif
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast2D(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
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
  const int num_species = m_physics->get_num_CdrSpecies();
  const RealVect origin = m_amr->getProbLo();

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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    if(solver->isDiffusive()){
      const int idx = solver_it.index();
      FArrayBox& dco = a_diffusionCoefficient_cell[idx]->getFArrayBox();
      dco.setVal(0.0);
      dco.copy(cdr_dco, a_box, idx, a_box, 0, 1);
      a_diffusionCoefficient_cell[idx]->setCoveredCellVal(0.0, 0);
    }
  }
#endif
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_reg_fast3D(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
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
  const int num_species = m_physics->get_num_CdrSpecies();
  const RealVect origin = m_amr->getProbLo();

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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    if(solver->isDiffusive()){
      const int idx = solver_it.index();
      FArrayBox& dco = a_diffusionCoefficient_cell[idx]->getFArrayBox();
      dco.setVal(0.0);
      dco.copy(cdr_dco, a_box, idx, a_box, 0, 1);

      a_diffusionCoefficient_cell[idx]->setCoveredCellVal(0.0, 0);
    }
  }
#endif
}

void cdr_plasma_stepper::compute_cdr_diffco_cell_irreg(Vector<EBCellFAB*>&          a_diffusionCoefficient_cell,
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
  const int num_species  = m_physics->get_num_CdrSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into cdr_plasma_physics
  RealVect         pos, E;
  Vector<Real>     cdr_densities(num_species);

  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
      
    const RealVect pos = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);
    const RealVect E   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));

    Vector<Real> cdr_densities(num_species);
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      cdr_densities[idx] = (*a_cdr_densities[idx])(vof, comp);
    }
      
    const Vector<Real> coeffs = m_physics->compute_cdr_diffusion_coefficients(a_time,
									      pos,
									      E,
									      cdr_densities);
      
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->isDiffusive()){
	(*a_diffusionCoefficient_cell[idx])(vof, comp) = coeffs[idx];
      }
    }
  }
}

void cdr_plasma_stepper::compute_cdr_diffco_eb(Vector<EBAMRIVData*>&       a_ebDiffusionCoefficient,
					       const Vector<EBAMRIVData*>& a_cdr_densities,
					       const EBAMRIVData&          a_E,
					       const Real&                 a_time){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_eb(full)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_eb(full)" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->get_num_CdrSpecies();
  const int finest_level = m_amr->getFinestLevel();
  const Real zero        = 0.0;

  Vector<LevelData<BaseIVFAB<Real> >* > diffco(num_species);
  Vector<LevelData<BaseIVFAB<Real> >* > cdr_densities(num_species);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<BaseIVFAB<Real> >* > diffco(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > cdr_densities(num_species);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx]        = (*a_ebDiffusionCoefficient[idx])[lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    // Call level versions
    compute_cdr_diffco_eb(diffco, cdr_densities, *a_E[lvl], a_time, lvl);
  }
}


void cdr_plasma_stepper::compute_cdr_diffco_eb(Vector<LevelData<BaseIVFAB<Real> >* >&       a_ebDiffusionCoefficient,
					       const Vector<LevelData<BaseIVFAB<Real> >* >& a_cdr_densities,
					       const LevelData<BaseIVFAB<Real> >&           a_E,
					       const Real&                                  a_time,
					       const int                                    a_lvl){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffco_eb(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffco_eb(level)" << endl;
  }

  const int comp = 0;
  const int num_species = m_physics->get_num_CdrSpecies();

  
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const RealVect origin         = m_amr->getProbLo();

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box            = dbl.get(dit());
    const EBISBox& ebisbox   = ebisl[dit()];
    const EBGraph& ebgraph   = ebisbox.getEBGraph();
    const IntVectSet ivs     = ebisbox.getIrregIVS(box);
    const BaseIVFAB<Real>& E = a_E[dit()];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const RealVect cntr = ebisbox.bndryCentroid(vof);
      const RealVect e    = RealVect(D_DECL(E(vof,0), E(vof,1), E(vof,2)));
      const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin) + cntr*dx;
      
      Vector<Real> cdr_densities(num_species);
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_cdr_densities[idx])[dit()](vof, 0);
	cdr_densities[idx] = Max(0.0, phi);
      }


      const Vector<Real> diffco = m_physics->compute_cdr_diffusion_coefficients(a_time,
										pos,
										e,
										cdr_densities);
										  
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	if(solver_it()->isDiffusive()){
#if 0 // Original code
	  (*a_ebDiffusionCoefficient[idx])[dit()](vof, comp) = diffco[idx];
#else // Debug code
	  (*a_ebDiffusionCoefficient[idx])[dit()](vof, comp) = 0.0;//1.E-5;
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

  const int num_species  = m_physics->get_num_CdrSpecies();
  const int num_Photons  = m_physics->get_num_RtSpecies();
  const int comp         = 0;
  const int ncomp        = 1;

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_Photons);

  // Grid stuff
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const ProblemDomain& domain   = m_amr->getDomains()[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const MFLevelGrid& mflg       = *(m_amr->getMFLevelGrid(m_realm)[a_lvl]);
  const RealVect origin         = m_amr->getProbLo();

  // Patch loop
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box              = dbl.get(dit());
    const EBISBox& ebisbox      = ebisl[dit()];
    const EBGraph& ebgraph      = ebisbox.getEBGraph();
    const IntVectSet& diel_ivs  = mflg.interfaceRegion(box, dit());
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
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](vof,comp);
	extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](vof,comp);
	extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](vof,comp);
	extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](vof,comp);
      }

      // Build Photon intensities
      for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
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
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
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
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](vof,comp);
	extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](vof,comp);
	extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](vof,comp);
	extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](vof,comp);
      }

      // Build Photon intensities
      for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
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
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
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

  const int num_species  = m_physics->get_num_CdrSpecies();
  const int num_Photons  = m_physics->get_num_RtSpecies();
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){


    Vector<LevelData<BaseIVFAB<Real> >* > fluxes(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_fluxes(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_densities(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_velocities(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_gradients(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_rte_fluxes(num_Photons);
    
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      fluxes[idx]                = (*a_fluxes[idx])[lvl];
      extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[lvl];
      extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[lvl];
      extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[lvl];
      extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[lvl];
    }


    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
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

  const int num_species  = m_physics->get_num_CdrSpecies();
  const int num_Photons  = m_physics->get_num_RtSpecies();
  const int finest_level = m_amr->getFinestLevel();

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_Photons);

  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<DomainFluxIFFAB>* > fluxes(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_fluxes(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_densities(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_velocities(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_gradients(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_rte_fluxes(num_Photons);
    
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      fluxes[idx]                = (*a_fluxes[idx])[lvl];
      extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[lvl];
      extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[lvl];
      extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[lvl];
      extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[lvl];
    }

    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
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

  const int num_species  = m_physics->get_num_CdrSpecies();
  const int num_Photons  = m_physics->get_num_RtSpecies();
  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_Photons);

  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const ProblemDomain& domain   = m_amr->getDomains()[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const RealVect origin         = m_amr->getProbLo();

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
	  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.index();
	    extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](dir, sit())(face,comp);
	  }

	  // Photon fluxes
	  for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
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
	  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.index();
	    (*a_fluxes[idx])[dit()](dir, sit())(face, comp) = fluxes[idx];
	  }
	}
      }
    }
  }
}

void cdr_plasma_stepper::computeGradients_at_eb(Vector<EBAMRIVData*>&         a_grad,
						 const phase::which_phase&     a_phase,
						 const Vector<EBAMRCellData*>& a_phi){
  CH_TIME("cdr_plasma_stepper::computeGradients_at_eb");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::computeGradients_at_eb" << endl;
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
    
    m_amr->computeGradient(gradient, density, m_realm, a_phase);       // Compute cell-centered gradient
    m_amr->averageDown(gradient, m_realm, a_phase);                    // Average down - shouldn't be necesasry
    m_amr->interpGhost(gradient, m_realm, a_phase);                    // Interpolate ghost cells (have to do this before interp)
    this->extrapolate_to_eb(eb_gradient, a_phase, gradient);            // Extrapolate to EB
    this->project_flux(grad_density, eb_gradient);                      // Project onto EB
  }
}

void cdr_plasma_stepper::computeGradients_at_domain_faces(Vector<EBAMRIFData*>&         a_grad,
							   const phase::which_phase&     a_phase,
							   const Vector<EBAMRCellData*>& a_phi){
  CH_TIME("cdr_plasma_stepper::computeGradients_at_eb");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::computeGradients_at_eb" << endl;
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
    
    m_amr->computeGradient(gradient, density, m_realm, a_phase);                         
    m_amr->averageDown(gradient, m_realm, a_phase);                             
    m_amr->interpGhost(gradient, m_realm, a_phase);                             
    
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

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
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
  
  const int num_species  = m_physics->get_num_CdrSpecies();
  const int finest_level = m_amr->getFinestLevel();

  // Interpolate E to centroids
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, phase::gas, SpaceDim);
  data_ops::copy(E, a_E);
  m_amr->interpToCentroids(E, m_realm, phase::gas);

  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<EBCellFAB>* > velocities(num_species);
    Vector<LevelData<EBCellFAB>* > cdr_densities(num_species);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      velocities[idx]    = (*a_velocities[idx])[lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    compute_cdr_velocities(velocities, cdr_densities, *E[lvl], lvl, a_time);
  }

  // Average down and interpolate ghost cells
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    if(solver_it()->isMobile()){
      m_amr->averageDown(*a_velocities[idx], m_realm, m_cdr->getPhase()); 
      m_amr->interpGhost(*a_velocities[idx], m_realm, m_cdr->getPhase()); 
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

  const phase::which_phase cdr_phase = m_cdr->getPhase();
  const int finest_level             = m_amr->getFinestLevel();

  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, cdr_phase)[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];

  const int num_species = m_physics->get_num_CdrSpecies();
    
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    Vector<EBCellFAB*> vel(num_species);
    Vector<EBCellFAB*> phi(num_species);;
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->isMobile()){
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
  const RealVect origin  = m_amr->getProbLo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();


  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv    = bit();
    const RealVect pos  = origin + a_dx*iv;
    const RealVect e    = RealVect(D_DECL(E(iv, 0), E(iv, 1), E(iv, 2)));


    // Get densities
    Vector<Real> cdr_densities;
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      cdr_densities.push_back((*a_cdr_densities[idx]).getSingleValuedFAB()(iv, comp));
    }

    // Compute velocities
    Vector<RealVect> velocities = m_physics->compute_cdr_velocities(a_time, pos, e, cdr_densities);

    // Put velocities in the appropriate place. 
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver = solver_it();
      const int idx = solver_it.index();
      if(solver->isMobile()){
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_velocities[idx]).getSingleValuedFAB()(iv, dir) = velocities[idx][dir];
	}
      }
    }
  }


  // Covered is bogus.
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->isMobile()){
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
  const int num_species  = m_physics->get_num_CdrSpecies();
  const int comp         = 0;
  const RealVect origin  = m_amr->getProbLo();
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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
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
  const int num_species  = m_physics->get_num_CdrSpecies();
  const int comp         = 0;
  const RealVect origin  = m_amr->getProbLo();
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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
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
  const RealVect origin  = m_amr->getProbLo();

  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);

    // Get densities
    Vector<Real> cdr_densities;
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      cdr_densities.push_back((*a_cdr_densities[idx])(vof, comp));
    }
    
    // Compute velocities
    const Vector<RealVect> velocities = m_physics->compute_cdr_velocities(a_time, pos, e, cdr_densities);

    // Put velocities in the appropriate place. 
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      if(solver_it()->isMobile()){
	const int idx = solver_it.index();
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_velocities[idx])(vof, dir) = velocities[idx][dir];
	}
      }
    }
  }
}

void cdr_plasma_stepper::preRegrid(const int a_lmin, const int a_finestLevel){
  CH_TIME("cdr_plasma_stepper::preRegrid");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::preRegrid" << endl;
  }

  // Solvers do pre-regridding shit. 
  m_cdr->preRegrid(a_lmin, a_finestLevel);
  m_fieldSolver->preRegrid(a_lmin, a_finestLevel);
  m_rte->preRegrid(a_lmin, a_finestLevel);
  m_sigma->preRegrid(a_lmin, a_finestLevel);
}

void cdr_plasma_stepper::preRegridInternals(const int a_lbase, const int a_finestLevel){
  CH_TIME("cdr_plasma_stepper::preRegridInternals");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::preRegridInternals" << endl;
  }
}

void cdr_plasma_stepper::compute_cdr_velocities(){
  CH_TIME("cdr_plasma_stepper::compute_cdr_velocities()");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_velocities()" << endl;
  }

  // Compute the electric field (again)
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, m_cdr->getPhase(), SpaceDim);
  this->compute_E(E, m_cdr->getPhase(), m_fieldSolver->getPotential());

  Vector<EBAMRCellData*> states     = m_cdr->getPhis();
  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();

  this->compute_cdr_velocities(velocities, states, E, m_time);
}

void cdr_plasma_stepper::compute_cdr_diffusion(){
  CH_TIME("cdr_plasma_stepper::compute_cdr_diffusion");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_cdr_diffusion" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_physics->get_num_CdrSpecies();

  EBAMRCellData E_cell;
  EBAMRIVData   E_eb;
  m_amr->allocate(E_cell, m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(E_eb,   m_realm, m_cdr->getPhase(), SpaceDim);
  
  this->compute_E(E_cell, m_cdr->getPhase(), m_fieldSolver->getPotential());
  this->compute_E(E_eb,   m_cdr->getPhase(), E_cell);

  Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();

  cdr_plasma_stepper::compute_cdr_diffusion(E_cell, E_eb);
#if 0

  // Extrapolate states to the EB
  Vector<EBAMRIVData*> cdr_extrap(num_species);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    cdr_extrap[idx] = new EBAMRIVData();  // This must be deleted
    m_amr->allocate(*cdr_extrap[idx], m_realm, m_cdr->getPhase(), ncomp);

    const IrregAmrStencil<EbCentroidInterpolationStencil>& stencil = m_amr->getEbCentroidInterpolationStencilStencils(m_realm, m_cdr->getPhase());
    stencil.apply(*cdr_extrap[idx], *cdr_states[idx]);
  }
  
  Vector<EBAMRFluxData*> diffco_face = m_cdr->getFaceCenteredDiffusionCoefficient();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->getEbCenteredDiffusionCoefficient();

  this->compute_cdr_diffco_face(diffco_face, cdr_states, E_cell, m_time);
  this->compute_cdr_diffco_eb(diffco_eb,     cdr_extrap, E_eb,   m_time);

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
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
  const int num_species = m_physics->get_num_CdrSpecies();

  Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();

  // Extrapolate states to the EB
  Vector<EBAMRIVData*> cdr_extrap(num_species);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    cdr_extrap[idx] = new EBAMRIVData();  // This must be deleted
    m_amr->allocate(*cdr_extrap[idx], m_realm, m_cdr->getPhase(), ncomp);

    const IrregAmrStencil<EbCentroidInterpolationStencil>& stencil = m_amr->getEbCentroidInterpolationStencilStencils(m_realm, m_cdr->getPhase());
    stencil.apply(*cdr_extrap[idx], *cdr_states[idx]);
  }
  
  Vector<EBAMRFluxData*> diffco_face = m_cdr->getFaceCenteredDiffusionCoefficient();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->getEbCenteredDiffusionCoefficient();

  this->compute_cdr_diffco_face(diffco_face, cdr_states, a_E_cell, m_time);
  this->compute_cdr_diffco_eb(diffco_eb,     cdr_extrap, a_E_eb,   m_time);

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
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

  m_amr->computeGradient(a_E, a_potential, m_realm);
  data_ops::scale(a_E, -1.0);

  m_amr->averageDown(a_E, m_realm);
  m_amr->interpGhost(a_E, m_realm);
}

void cdr_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase){
  CH_TIME("cdr_plasma_stepper::compute_E(ebamrcell, phase)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_E(ebamrcell, phase)" << endl;
  }

  this->compute_E(a_E, a_phase, m_fieldSolver->getPotential());
}

void cdr_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential){
  CH_TIME("cdr_plasma_stepper::compute_E(ebamrcell, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_E(ebamrcell, phase, mfamrcell)" << endl;
  }

  EBAMRCellData pot_gas;
  m_amr->allocatePointer(pot_gas);
  m_amr->alias(pot_gas, a_phase, a_potential);

  m_amr->computeGradient(a_E, pot_gas, m_realm, a_phase);
  data_ops::scale(a_E, -1.0);

  m_amr->averageDown(a_E, m_realm, a_phase);
  m_amr->interpGhost(a_E, m_realm, a_phase);
}

void cdr_plasma_stepper::compute_E(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("cdr_plasma_stepper::compute_E(ebamrflux, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_E(ebamrflux, phase, ebamrcell)" << endl;
  }

  CH_assert(a_E_face[0]->nComp() == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, a_phase)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];

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
    //    data_ops::average_cell_to_face_allcomps(*a_E_face[lvl], *a_E_cell[lvl], m_amr->getDomains()[lvl]);
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

  const IrregAmrStencil<EbCentroidInterpolationStencil>& interp_stencil = m_amr->getEbCentroidInterpolationStencilStencils(m_realm, a_phase);
  interp_stencil.apply(a_E_eb, a_E_cell);
}

void cdr_plasma_stepper::compute_Emax(Real& a_Emax, const phase::which_phase a_phase){
  CH_TIME("cdr_plasma_stepper::compute_Emax(Real, phase)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_Emax(Real, phase)" << endl;
  }

  EBAMRCellData E;
  m_amr->allocate(E, m_realm, a_phase, SpaceDim);

  this->compute_E(E, a_phase, m_fieldSolver->getPotential());
  m_amr->interpToCentroids(E, m_realm, a_phase);

  Real max, min;
  data_ops::get_max_min_norm(max, min, E);

  a_Emax = max;
}

void cdr_plasma_stepper::computeCharge_flux(EBAMRIVData& a_flux, Vector<EBAMRIVData*>& a_cdr_fluxes){
  CH_TIME("cdr_plasma_stepper::computeCharge_flux");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::computeCharge_flux" << endl;
  }

  MayDay::Abort("cdr_plasma_stepper::computeCharge_flux - I'm suspecting that this is deprecated code which is no longer used");

  data_ops::set_value(a_flux, 0.0);

  for (CdrIterator<CdrSolver> solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = *a_cdr_fluxes[solver_it.index()];

    data_ops::incr(a_flux, solver_flux, spec->getChargeNumber()*units::s_Qe);
  }

  m_sigma->resetCells(a_flux);
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
    this->computeFlux(cell_flux, *a_densities[i], *a_velocities[i]);
    
    m_amr->averageDown(cell_flux, m_realm, a_phase);
    m_amr->interpGhost(cell_flux, m_realm, a_phase);

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

  const IrregAmrStencil<EbCentroidInterpolationStencil>& interp_stencils = m_amr->getEbCentroidInterpolationStencilStencils(m_realm, a_phase);

  //  for (int i = 0; i < a_fluxes.size(); i++){
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      interp_stencils.apply(eb_vel, *a_velocities[idx]);
      interp_stencils.apply(eb_phi, *a_densities[idx]);

      data_ops::floor(eb_phi, 0.0);

      data_ops::set_value(eb_flx, 0.0);
      data_ops::incr(eb_flx, eb_vel, 1.0);
      data_ops::multiply_scalar(eb_flx, eb_phi);

      this->project_flux(*a_fluxes[idx], eb_flx);

      m_amr->averageDown(*a_fluxes[idx], m_realm, a_phase);
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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
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

void cdr_plasma_stepper::computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_density, const EBAMRCellData& a_velocity){
  CH_TIME("cdr_plasma_stepper::computeFlux(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::computeFlux(amr)" << endl;
  }
  
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    CH_assert(a_flux[lvl]->nComp()     == SpaceDim);
    CH_assert(a_density[lvl]->nComp()  == 1);
    CH_assert(a_velocity[lvl]->nComp() == SpaceDim);

    // Call level versions
    computeFlux(*a_flux[lvl], *a_density[lvl], *a_velocity[lvl], lvl);
  }
}

void cdr_plasma_stepper::computeFlux(LevelData<EBCellFAB>&       a_flux,
				      const LevelData<EBCellFAB>& a_density,
				      const LevelData<EBCellFAB>& a_velocity,
				      const int                   a_lvl){
  CH_TIME("cdr_plasma_stepper::computeFlux(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::computeFlux(amr)" << endl;
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
  
  const int finest_level = m_amr->getFinestLevel();

  // Call level versions
  for (int lvl = 0; lvl <= finest_level; lvl++){
    compute_J(*a_J[lvl], lvl);
  }

  m_amr->averageDown(a_J, m_realm, m_cdr->getPhase());
  m_amr->interpGhost(a_J, m_realm, m_cdr->getPhase());
}

void cdr_plasma_stepper::compute_J(LevelData<EBCellFAB>& a_J, const int a_lvl) const{
  CH_TIME("cdr_plasma_stepper::compute_J(level)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_J(level)" << endl;
  }

  data_ops::set_value(a_J, 0.0);
  
  const int density_comp = 0;

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();

    if(solver->isMobile()){
      const int q                       = spec->getChargeNumber();
      const EBAMRCellData& density      = solver->getPhi();
      const EBAMRCellData& velo         = solver->getCellCenteredVelocity();

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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver> solver = solver_it();
    densities.push_back(&(solver->getPhi()));
  }

  this->compute_rho(m_fieldSolver->getRho(), densities, centering::cell_center);
}

void cdr_plasma_stepper::compute_rho(EBAMRCellData& a_rho, const phase::which_phase a_phase){
  CH_TIME("cdr_plasma_stepper::compute_rho(ebamrcelldata, phase)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_rho(ebamrcelldata, phase)" << endl;
  }

  CH_assert(a_phase == m_cdr->getPhase());

  data_ops::set_value(a_rho, 0.0);

  Vector<EBAMRCellData*> densities = m_cdr->getPhis(); // Get densities from solver
  
  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    // Add volumetric charge 
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const EBAMRCellData& density       = *(densities[solver_it.index()]);
      const RefCountedPtr<CdrSpecies>& spec = solver_it.getSpecies();

      if(spec->getChargeNumber() != 0){
	data_ops::incr(*a_rho[lvl], *density[lvl], spec->getChargeNumber());
      }
    }

    // Scale by s_Qe/s_eps0
    data_ops::scale(*a_rho[lvl], units::s_Qe);
  }

  m_amr->averageDown(a_rho, m_realm, a_phase);
  m_amr->interpGhost(a_rho, m_realm, a_phase);
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
  m_amr->allocatePointer(rho_gas); 
  m_amr->alias(rho_gas, phase::gas, a_rho); 
  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    // Add volumetric charge 
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const EBAMRCellData& density       = *(a_densities[solver_it.index()]);
      const RefCountedPtr<CdrSpecies>& spec = solver_it.getSpecies();

      if(spec->getChargeNumber() != 0){
	data_ops::incr(*rho_gas[lvl], *density[lvl], spec->getChargeNumber());
      }
    }

    // Scale by s_Qe
    data_ops::scale(*a_rho[lvl], units::s_Qe);
  }

  m_amr->averageDown(a_rho, m_realm);
  m_amr->interpGhost(a_rho, m_realm);

  // Transform to centroids
  if(a_centering == centering::cell_center){
    m_amr->interpToCentroids(rho_gas, m_realm, phase::gas);
  }
}

void cdr_plasma_stepper::deallocate_solver_internals(){
  CH_TIME("cdr_plasma_stepper::deallocate_solver_internals");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::deallocate_solver_internals" << endl;
  }

  m_cdr->deallocateInternals();
  m_rte->deallocateInternals();
  m_fieldSolver->deallocateInternals();
  m_sigma->deallocateInternals();
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

  const IrregAmrStencil<EbCentroidInterpolationStencil>& stencils = m_amr->getEbCentroidInterpolationStencilStencils(m_realm, a_phase);
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
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

  const IrregAmrStencil<EbCentroidInterpolationStencil>& stencils = m_amr->getEbCentroidInterpolationStencilStencils(m_realm, a_phase);
  stencils.apply(a_extrap, a_data, a_lvl);
}

void cdr_plasma_stepper::extrapolate_to_domain_faces(EBAMRIFData&             a_extrap,
						     const phase::which_phase a_phase,
						     const EBAMRCellData&     a_data){
  CH_TIME("cdr_plasma_stepper::extrapolate_to_domain_faces(amr)");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::extrapolate_to_domain_faces(amr)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
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
      
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, a_phase)[a_lvl];
  
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
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    Real max, min;
    data_ops::get_max_min(max, min, solver->getPhi(), comp);

    if(max > a_cdr_max){
      a_cdr_max     = max;
      a_solver_name = solver->getName();
    }
  }
}

void cdr_plasma_stepper::set_cdr(RefCountedPtr<cdr_layout<CdrSolver>>& a_cdr){
  CH_TIME("cdr_plasma_stepper::set_cdr");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_cdr" << endl;
  }
  m_cdr = a_cdr;
}

void cdr_plasma_stepper::set_poisson(RefCountedPtr<FieldSolver>& a_poisson){
  CH_TIME("cdr_plasma_stepper::set_poisson");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_poisson" << endl;
  }
  m_fieldSolver = a_poisson;
}

void cdr_plasma_stepper::set_rte(RefCountedPtr<RtLayout<RtSolver>>& a_rte){
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
  parseOptions();
  this->sanityCheck();

  // Make solvers
  this->setup_cdr();
  this->setup_rte(); 
  this->setup_poisson();
  this->setup_sigma();

  this->set_solver_verbosity();

  // Allocate internal memory
  this->allocateInternals();
}

void cdr_plasma_stepper::initialData(){
  CH_TIME("cdr_plasma_stepper::initialData");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::initialData" << endl;
  }


  m_cdr->initialData();        // Initial data comes in through CdrSpecies, in this case supplied by physics
  if(!m_rte->isStationary()){
    m_rte->initialData();
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

  const RealVect origin  = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  EBAMRIVData& sigma = m_sigma->getPhi();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, phase::gas)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    
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

  m_amr->averageDown(sigma, m_realm, phase::gas);
  m_sigma->resetCells(sigma);
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

  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box         = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    
    BaseIVFAB<Real>& proj_flux  = a_projected_flux[dit()];
    const BaseIVFAB<Real>& flux = a_flux[dit()];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
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
  const int finest_level = m_amr->getFinestLevel();

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
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    CH_assert(a_projected_flux[lvl]->nComp() == 1);
    CH_assert(a_flux[lvl]->nComp()           == SpaceDim);

    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[lvl];

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

  this->allocateInternals(); // Allocate memory for time stepper
  this->regrid_solvers(a_lmin, a_old_finest, a_new_finest);
  this->regridInternals(a_lmin, a_old_finest, a_new_finest);

  // Solvers have been regridded. Now resolve the Poisson equation with the new data
  bool converged = this->solve_poisson();

  // If we don't converge, try new Poisson solver settings
  if(!converged){ 
    if(m_verbosity > 0){
      pout() << "Driver::regrid - Poisson solver failed to converge." << endl;
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
  m_fieldSolver->regrid(a_lmin, a_old_finest, a_new_finest);
  m_rte->regrid(a_lmin,     a_old_finest, a_new_finest);
  m_sigma->regrid(a_lmin,   a_old_finest, a_new_finest);
}

void cdr_plasma_stepper::reset_dielectric_cells(EBAMRIVData& a_data){
  CH_TIME("cdr_plasma_stepper::reset_dielectric_cells");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::reset_dielectric_cells" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const MFLevelGrid& mflg      = *m_amr->getMFLevelGrid(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      BaseIVFAB<Real>& data  = (*a_data[lvl])[dit()];
      const IntVectSet ivs   = data.getIVS() & mflg.interfaceRegion(box, dit());
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

void cdr_plasma_stepper::sanityCheck(){
  CH_TIME("cdr_plasma_stepper::sanityCheck");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::sanityCheck" << endl;
  }

  CH_assert(!m_computationalGeometry.isNull());
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

void cdr_plasma_stepper::set_potential(std::function<Real(const Real a_time)> a_potential){
  CH_TIME("cdr_plasma_stepper::set_potential");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::set_potential" << endl;
  }

  m_potential = a_potential;
}

void cdr_plasma_stepper::parseVerbosity(){
  CH_TIME("cdr_plasma_stepper::parseVerbosity");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::parseVerbosity" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  pp.get("verbosity", m_verbosity);
}

void cdr_plasma_stepper::parse_solver_verbosity(){
  CH_TIME("cdr_plasma_stepper::parse_solver_verbosity");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::parse_solver_verbosity" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  pp.get("solver_verbosity", m_solver_verbosity);
}

void cdr_plasma_stepper::parse_cfl(){

  ParmParse pp(m_className.c_str());
  pp.get("cfl", m_cfl);
  if(m_cfl < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_cfl - CFL cannot be negative!");
  }
}

void cdr_plasma_stepper::parse_relax_time(){
  ParmParse pp(m_className.c_str());
  pp.get("relax_time", m_relax_time);
  if(m_relax_time < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_relax_time - relaxation time cannot be negative");
  }
}

void cdr_plasma_stepper::parse_source_growth(){
  ParmParse pp(m_className.c_str());
  pp.query("source_growth", m_src_growth);
  if(m_src_growth < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_source_growth - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_source_tolerance(){
  ParmParse pp(m_className.c_str());
  pp.query("source_tolerance", m_src_tolerance);
  if(m_src_tolerance < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_source_tolerance - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_min_dt(){
  ParmParse pp(m_className.c_str());
  pp.get("min_dt", m_min_dt);
  if(m_min_dt < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_min_dt - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_max_dt(){
  ParmParse pp(m_className.c_str());
  pp.get("max_dt", m_max_dt);
  if(m_max_dt < 0.0){
    MayDay::Abort("cdr_plasma_stepper::parse_max_dt - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_fast_rte(){
  ParmParse pp(m_className.c_str());
  pp.get("fast_rte", m_fast_rte);
  if(m_fast_rte <= 0){
    MayDay::Abort("cdr_plasma_stepper::parse_fast_rte - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_fast_poisson(){
  ParmParse pp(m_className.c_str());
  pp.get("fast_poisson", m_fast_poisson);
  if(m_fast_poisson <= 0){
    MayDay::Abort("cdr_plasma_stepper::parse_fast_poisson - value cannot be negative");
  }
}

void cdr_plasma_stepper::parse_source_comp(){

  std::string str;
  ParmParse pp(m_className.c_str());
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
    m_cdr->setVerbosity(m_solver_verbosity);
  }
  if(!m_fieldSolver.isNull()){
    m_fieldSolver->setVerbosity(m_solver_verbosity);
  }
  if(!m_rte.isNull()){
    m_rte->setVerbosity(m_solver_verbosity);
  }
  if(!m_sigma.isNull()){
    m_sigma->setVerbosity(m_solver_verbosity);
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

void cdr_plasma_stepper::setSource_computation(){
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
      MayDay::Abort("cdr_plasma_stepper::setSource_computation - unknown type requested");
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

void cdr_plasma_stepper::setSource_growth(const Real a_src_growth){
  m_src_growth = a_src_growth;

  ParmParse pp("cdr_plasma_stepper");
  pp.query("source_growth", m_src_growth);
  if(m_src_growth < 0.0){
    m_src_growth = a_src_growth;
  }
}

void cdr_plasma_stepper::setSource_growth_tolerance(const Real a_src_tolerance){
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

  m_cdr->setVerbosity(m_solver_verbosity);
  m_cdr->parseOptions();
  m_cdr->setAmr(m_amr);
  m_cdr->setComputationalGeometry(m_computationalGeometry);
  m_cdr->setPhase(phase::gas);
  m_cdr->sanityCheck();
  m_cdr->setRealm(m_realm);
}

void cdr_plasma_stepper::allocate() {
  m_cdr->allocateInternals();
  m_fieldSolver->allocateInternals();
  m_rte->allocateInternals();
  m_sigma->allocateInternals();
}

void cdr_plasma_stepper::setup_poisson(){
  CH_TIME("cdr_plasma_stepper::setup_poisson");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::setup_poisson" << endl;
  }

  m_fieldSolver->setVerbosity(m_solver_verbosity);
  m_fieldSolver->parseOptions();
  m_fieldSolver->setAmr(m_amr);
  m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
  m_fieldSolver->setRealm(m_realm);
  m_fieldSolver->setVoltage(m_potential); // Needs to happen AFTER set_poisson_wall_func
}

void cdr_plasma_stepper::setup_rte(){
  CH_TIME("cdr_plasma_stepper::setup_rte");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::setup_rte" << endl;
  }

  m_rte->setVerbosity(m_solver_verbosity);
  m_rte->parseOptions();
  m_rte->setPhase(phase::gas);
  m_rte->setAmr(m_amr);
  m_rte->setComputationalGeometry(m_computationalGeometry);
  m_rte->sanityCheck();
  m_rte->setRealm(m_realm);
  //  m_rte->allocateInternals();
}

void cdr_plasma_stepper::setup_sigma(){
  CH_TIME("cdr_plasma_stepper::setup_sigma");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::setup_sigma" << endl;
  }

  m_sigma = RefCountedPtr<SigmaSolver> (new SigmaSolver());
  m_sigma->setAmr(m_amr);
  m_sigma->setVerbosity(m_solver_verbosity);
  m_sigma->setComputationalGeometry(m_computationalGeometry);
  m_sigma->setRealm(m_realm);
  //  m_sigma->allocateInternals();
}

void cdr_plasma_stepper::solver_dump(){
  CH_TIME("cdr_plasma_stepper::solver_dump");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::solver_dump" << endl;
  }

  m_cdr->writePlotFile();
  m_fieldSolver->writePlotFile();
  m_rte->writePlotFile();
}

void cdr_plasma_stepper::solve_rte(const Real a_dt){
  CH_TIME("cdr_plasma_stepper::solve_rte()");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::solve_rte()" << endl;
  }

  const phase::which_phase rte_phase = m_rte->getPhase();

  EBAMRCellData E;
  m_amr->allocate(E, m_realm, rte_phase, SpaceDim);
  this->compute_E(E, rte_phase, m_fieldSolver->getPotential());

  Vector<EBAMRCellData*> states     = m_rte->getPhis();
  Vector<EBAMRCellData*> rhs        = m_rte->getSources();
  Vector<EBAMRCellData*> cdr_states = m_cdr->getPhis();

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

  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    
    RefCountedPtr<RtSolver>& solver = solver_it();
    EBAMRCellData& state              = *a_rte_states[idx];
    EBAMRCellData& rhs                = *a_rte_sources[idx];
    solver->advance(a_dt, state, rhs);
  }
}

void cdr_plasma_stepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("cdr_plasma_stepper::synchronizeSolverTimes");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;

  m_cdr->setTime(a_step,     a_time, a_dt);
  m_fieldSolver->setTime(a_step, a_time, a_dt);
  m_rte->setTime(a_step,     a_time, a_dt);
  m_sigma->setTime(a_step,   a_time, a_dt);
}

Real cdr_plasma_stepper::compute_electrode_current(){
  CH_TIME("cdr_plasma_stepper::compute_electrode_current");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_stepper::compute_electrode_current" << endl;
  }

  // Need to copy onto temporary storage because 
  EBAMRIVData charge_flux;
  m_amr->allocate(charge_flux, m_realm, m_cdr->getPhase(), 1);
  data_ops::set_value(charge_flux, 0.0);
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    data_ops::incr(charge_flux, solver_flux, spec->getChargeNumber()*units::s_Qe);
  }

  this->reset_dielectric_cells(charge_flux);
  m_amr->conservativeAverage(charge_flux, m_realm, m_cdr->getPhase());

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->getDx()[compute_lvl];
  for (DataIterator dit = m_amr->getGrids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
    const BaseIVFAB<Real>& flx = (*charge_flux[compute_lvl])[dit()];

    const IntVectSet ivs = flx.getIVS() & m_amr->getGrids(m_realm)[compute_lvl].get(dit());
    for (VoFIterator vofit(ivs, flx.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& vof   = vofit();
      const Real& bndryFrac = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[compute_lvl][dit()].bndryArea(vof);
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
  m_amr->allocate(charge_flux, m_realm, m_cdr->getPhase(), 1);
  data_ops::set_value(charge_flux, 0.0);
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    data_ops::incr(charge_flux, solver_flux, spec->getChargeNumber()*units::s_Qe);
  }

  m_sigma->resetCells(charge_flux);
  m_amr->conservativeAverage(charge_flux, m_realm, m_cdr->getPhase());

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->getDx()[compute_lvl];
  for (DataIterator dit = m_amr->getGrids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
    const BaseIVFAB<Real>& flx = (*charge_flux[compute_lvl])[dit()];

    const IntVectSet ivs = flx.getIVS() & m_amr->getGrids(m_realm)[compute_lvl].get(dit());
    for (VoFIterator vofit(ivs, flx.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& vof   = vofit();
      const Real& bndryFrac = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[compute_lvl][dit()].bndryArea(vof);
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
  m_amr->allocate(charge_flux, m_realm, m_cdr->getPhase(), 1);
  data_ops::set_value(charge_flux, 0.0);
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();
    const EBAMRIFData& solver_flux          = solver->getDomainFlux();

    data_ops::incr(charge_flux, solver_flux, spec->getChargeNumber()*units::s_Qe);
  }

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->getDx()[compute_lvl];
  for (DataIterator dit = m_amr->getGrids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
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
  m_amr->allocate(J,     m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(E,     m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(JdotE, m_realm, m_cdr->getPhase(), SpaceDim);

  this->compute_E(E, m_cdr->getPhase(), m_fieldSolver->getPotential());
  this->compute_J(J);

  // Compute J.dot.E 
  data_ops::dot_prod(JdotE, J,E);
  m_amr->averageDown(JdotE, m_realm, m_cdr->getPhase());

  // Only compue on coarsest level
  const int coar = 0;
  const Real dx  = m_amr->getDx()[coar];
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
  m_amr->allocate(E,  m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(J,  m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(dt, m_realm, m_cdr->getPhase(), 1);

  data_ops::set_value(dt, 1.234567E89);

  this->compute_E(E, m_cdr->getPhase(), m_fieldSolver->getPotential());
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
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[lvl];

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
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
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

Real cdr_plasma_stepper::getTime(){
  return m_time;
}

Real cdr_plasma_stepper::get_dt(){
  return m_dt;
}

Real cdr_plasma_stepper::get_cfl_dt(){
  return m_dt_cfl;
}

RefCountedPtr<cdr_layout<CdrSolver>>& cdr_plasma_stepper::get_cdr(){
  return m_cdr;
}

RefCountedPtr<FieldSolver>& cdr_plasma_stepper::get_poisson(){
  return m_fieldSolver;
}

RefCountedPtr<RtLayout<RtSolver>>& cdr_plasma_stepper::get_rte(){
  return m_rte;
}

RefCountedPtr<SigmaSolver>& cdr_plasma_stepper::get_sigma(){
  return m_sigma;
}

// New functions for Driver
void cdr_plasma_stepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const{
  CH_TIME("Driver::writeCheckpointData");
  if(m_verbosity > 3){
    pout() << "Driver::writeCheckpointData" << endl;
  }

  // CDR solvers checkpoint their data
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  // RTE solvers checkpoint their data
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<RtSolver>& solver = solver_it();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->writeCheckpointLevel(a_handle, a_lvl);
  m_sigma->writeCheckpointLevel(a_handle, a_lvl);
}

void cdr_plasma_stepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("Driver::readCheckpointData");
  if(m_verbosity > 3){
    pout() << "Driver::readCheckpointData" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver = solver_it();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->readCheckpointLevel(a_handle, a_lvl);
  m_sigma->readCheckpointLevel(a_handle, a_lvl);
}

int cdr_plasma_stepper::getNumberOfPlotVariables() const{
  CH_TIME("cdr_plasma_stepper::getNumberOfPlotVariables");
  if(m_verbosity > 3){
    pout() << "cdr_plasma_stepper::getNumberOfPlotVariables" << endl;
  }
  int ncomp = 0;
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    ncomp += solver->getNumberOfPlotVariables();
  }
  
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver = solver_it();
    ncomp += solver->getNumberOfPlotVariables();
  }

  ncomp += m_fieldSolver->getNumberOfPlotVariables();
  ncomp += m_sigma->getNumberOfPlotVariables();
  ncomp += SpaceDim; // For plotting the current density

  return ncomp;
}

void cdr_plasma_stepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  CH_TIME("cdr_plasma_stepper::writePlotData");
  if(m_verbosity > 3){
    pout() << "cdr_plasma_stepper::writePlotData" << endl;
  }

  // Poisson solver copies over its output data
  a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  m_fieldSolver->writePlotData(a_output, a_icomp);

  // Surface charge solver writes
  a_plotVariableNames.append(m_sigma->getPlotVariableNames());
  m_sigma->writePlotData(a_output, a_icomp);

  // CDR solvers copy their output data
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver = solver_it();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // Write the current to the output
  this->write_J(a_output, a_icomp);
  a_plotVariableNames.push_back("x-J");
  a_plotVariableNames.push_back("y-J");
  if(SpaceDim == 3){
    a_plotVariableNames.push_back("z-J");
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
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }
  a_icomp += SpaceDim;
}

void cdr_plasma_stepper::postCheckpointSetup(){
  CH_TIME("cdr_plasma_stepper::postCheckpointSetup");
  if(m_verbosity > 3){
    pout() << "cdr_plasma_stepper::postCheckpointSetup" << endl;
  }

  this->solve_poisson();       // Solve Poisson equation by 
  if(this->stationary_rte()){  // Solve RTE equations if stationary solvers
    const Real dummy_dt = 0.0;
    this->solve_rte(dummy_dt); // Argument does not matter, it's a stationary solver.
  }
  this->allocateInternals();  // Prepare internal storage for time stepper

  // Fill solvers with important stuff
  this->compute_cdr_velocities();
  this->compute_cdr_diffusion();
}

void cdr_plasma_stepper::printStepReport(){
  CH_TIME("cdr_plasma_stepper::printStepReport");
  if(m_verbosity > 4){
    pout() << "cdr_plasma_stepper::printStepReport" << endl;
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
  if(m_timeCode == TimeCode::Advection){
    str = " (Restricted by advection)";
  }
  else if(m_timeCode == TimeCode::AdvectionDiffusion){
    str = " (Restricted by advection-diffusion)";
  }
  else if(m_timeCode == TimeCode::Error){
    str = " (Restricted by error)";
  }
  else if(m_timeCode == TimeCode::Diffusion){
    str = " (Restricted by diffusion)";
  }
  else if(m_timeCode == TimeCode::Source){
    MayDay::Abort("Driver::stepReport - shouldn't happen, source term has been taken out of the design");
    str = " (Restricted by source term)";
  }
  else if(m_timeCode == TimeCode::RelaxationTime){
    str = " (Restricted by relaxation time)";
  }
  else if(m_timeCode == TimeCode::Restricted){
    str = " (Restricted by time stepper)";
  }
  else if(m_timeCode == TimeCode::Hardcap){
    str = " (Restricted by a hardcap)";
  }
  pout() << "                                   mode  = " << str << endl
	 << "                                   cfl   = " << m_dt/cfl_dt << endl
	 << "                                   Emax  = " << Emax << endl
	 << "                                   n_max = " << nmax << "(" + solver_max + ")" << endl;
}
#include "CD_NamespaceFooter.H"
