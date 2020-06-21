/*!
  @file   ito_plasma_godunov.cpp
  @author Robert Marskar
  @date   June 2020
  @brief  Implementation of ito_plasma_godunov
*/

#include "ito_plasma_godunov.H"

#include <ParmParse.H>

#define DEBUG 1

using namespace physics::ito_plasma;

ito_plasma_godunov::ito_plasma_godunov(){
  m_name = "ito_plasma_godunov";
}

ito_plasma_godunov::ito_plasma_godunov(RefCountedPtr<ito_plasma_physics>& a_physics){
  m_name    = "ito_plasma_godunov";
  m_physics = a_physics;
}

ito_plasma_godunov::~ito_plasma_godunov(){

}

void ito_plasma_godunov::parse_options() {
  CH_TIME("ito_plasma_godunov::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("merge_interval", m_merge_interval);
}

void ito_plasma_godunov::allocate_internals(){
  CH_TIME("ito_plasma_godunov::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }
}

Real ito_plasma_godunov::advance(const Real a_dt) {
  CH_TIME("ito_plasma_godunov::advance");
  if(m_verbosity > 5){
    pout() << m_name + "::advance" << endl;
  }

  Real t_total = -MPI_Wtime();
  
  Real t_advect = 0.0;
  Real t_diffuse = 0.0;
  Real t_isect = 0.0;
  Real t_remap = 0.0;
  Real t_deposit = 0.0;
  Real t_photons = 0.0;
  Real t_poisson = 0.0;
  Real t_cellBin = 0.0;
  Real t_chemistry = 0.0;
  Real t_super = 0.0;
  Real t_patchBin = 0.0;
  Real t_clear = 0.0;
  Real t_velo = 0.0;
  Real t_diffu = 0.0;
    

  // Transport kernel for particles and photons.
  MPI_Barrier(Chombo_MPI::comm);
  t_advect -= MPI_Wtime();
  this->advect_particles(a_dt);
  t_advect += MPI_Wtime();
  t_diffuse -= MPI_Wtime();
  this->diffuse_particles(a_dt);
  t_diffuse += MPI_Wtime();
  t_isect -= MPI_Wtime();
  this->intersect_particles(a_dt);
  t_isect += MPI_Wtime();
  t_remap -= MPI_Wtime();
  m_ito->remap();
  t_remap += MPI_Wtime();
  t_deposit -= MPI_Wtime();
  m_ito->deposit_particles();
  MPI_Barrier(Chombo_MPI::comm);
  t_deposit += MPI_Wtime();

  // Remove particles that leaked into EBs
#if 1 // Write this in a better format, please!
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->remove_eb_particles();
  }
#endif


  // Move photons
  MPI_Barrier(Chombo_MPI::comm);
  t_photons -= MPI_Wtime();
  this->advance_photons(a_dt);
  MPI_Barrier(Chombo_MPI::comm);
  t_photons += MPI_Wtime();



  // Sort the particles per cell
  MPI_Barrier(Chombo_MPI::comm);
  t_cellBin -= MPI_Wtime();
  m_ito->sort_particles_by_cell();
  this->sort_bulk_photons_by_cell();
  this->sort_source_photons_by_cell();
  MPI_Barrier(Chombo_MPI::comm);
  t_cellBin += MPI_Wtime();

  // Chemistry kernel.
  MPI_Barrier(Chombo_MPI::comm);
  t_chemistry -= MPI_Wtime();
  this->advance_reaction_network(a_dt);
  MPI_Barrier(Chombo_MPI::comm);
  t_chemistry += MPI_Wtime();

  // Make superparticles
  MPI_Barrier(Chombo_MPI::comm);
  t_super -= MPI_Wtime();
  if((m_step+1) % m_merge_interval == 0 && m_merge_interval > 0){
    m_ito->make_superparticles(m_ppc);
  }
  MPI_Barrier(Chombo_MPI::comm);
  t_super += MPI_Wtime();

  // Sort particles per patch
  MPI_Barrier(Chombo_MPI::comm);
  t_patchBin -= MPI_Wtime();
  m_ito->sort_particles_by_patch();
  this->sort_bulk_photons_by_patch();
  this->sort_source_photons_by_patch();
  MPI_Barrier(Chombo_MPI::comm);
  t_patchBin += MPI_Wtime();

  // Compute the electric field, recompute velocities and diffusion coefficients
  MPI_Barrier(Chombo_MPI::comm);
  t_poisson -= MPI_Wtime();
  this->solve_poisson();
  MPI_Barrier(Chombo_MPI::comm);
  t_poisson += MPI_Wtime();

  // Clear other data holders for now. BC comes later
  MPI_Barrier(Chombo_MPI::comm);
  t_clear -= MPI_Wtime();
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->clear(solver_it()->get_eb_particles());
    solver_it()->clear(solver_it()->get_domain_particles());
  }
  MPI_Barrier(Chombo_MPI::comm);
  t_clear += MPI_Wtime();

  // Prepare next step
  MPI_Barrier(Chombo_MPI::comm);
  t_velo -= MPI_Wtime();
  this->compute_ito_velocities();
  MPI_Barrier(Chombo_MPI::comm);
  t_velo += MPI_Wtime();
  t_diffu -= MPI_Wtime();
  this->compute_ito_diffusion();
  MPI_Barrier(Chombo_MPI::comm);
  t_diffu += MPI_Wtime();

  t_total += MPI_Wtime();
  
#if DEBUG
  t_advect    *= 100./t_total;
  t_diffuse   *= 100./t_total;
  t_isect     *= 100./t_total;
  t_remap     *= 100./t_total;
  t_deposit   *= 100./t_total;
  t_photons   *= 100./t_total;
  t_poisson   *= 100./t_total;
  t_chemistry *= 100./t_total;
  t_super     *= 100./t_total;
  t_clear     *= 100./t_total;
  t_velo      *= 100./t_total;
  t_diffu     *= 100./t_total;
  pout() << "total time = " << t_total << "\n"
	 << "advect = " << t_advect << "%\n"
    	 << "diffuse = " << t_diffuse << "%\n"
	 << "isect = " << t_isect << "%\n"
    	 << "remap = " << t_remap << "%\n"
	 << "deposit = " << t_deposit << "%\n"
    	 << "photons = " << t_photons << "%\n"
	 << "poisson = " << t_poisson << "%\n"
	 << "cell bin = " << t_cellBin << "%\n"
    	 << "chemistry = " << t_chemistry << "%\n"
	 << "super = " << t_super << "%\n"
    	 << "patch bin = " << t_patchBin << "%\n"
	 << "clear = " << t_clear << "%\n"
    	 << "velo = " << t_velo << "%\n"
	 << "diffu = " << t_diffu << "%\n"
	 << endl;
#endif


  
  return a_dt;
}

void ito_plasma_godunov::advect_particles(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advect_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    if(solver->is_mobile()){
      
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids()[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // First step
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];
	    p.oldPosition() = p.position();
	    p.position() += 0.5*p.velocity()*a_dt;
	  }

	  // Interpolate velocities
	  solver->interpolate_velocities(lvl, dit());

	  // Second step
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];
	    p.position() = p.oldPosition() + p.velocity()*a_dt;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::diffuse_particles(const Real a_dt){
  CH_TIME("ito_plasma_godunov::diffuse_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::diffuse_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    
    if(solver->is_diffusive()){
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids()[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // Diffusion hop.
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];
	    const RealVect ran = solver->random_gaussian();
	    const RealVect hop = ran*sqrt(2.0*p.diffusion()*a_dt);
	    p.position() += hop;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::intersect_particles(const Real a_dt){
  CH_TIME("ito_plasma_godunov::intersect_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::intersect_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();

    solver->intersect_particles();
  }
}
