/*!
  @file   ito_plasma_godunov.cpp
  @author Robert Marskar
  @date   June 2020
  @brief  Implementation of ito_plasma_godunov
*/

#include "ito_plasma_godunov.H"

#include <ParmParse.H>

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

  pp.get("verbosity",     m_verbosity);
  pp.get("ppc",           m_ppc);
  pp.get("max_cells_hop", m_max_cells_hop);
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

  // Transport kernel for particles and photons. 
  this->advect_particles(a_dt);
  this->diffuse_particles(a_dt);
  this->intersect_particles(a_dt);
  m_ito->remap();
  m_ito->deposit_particles();

  // Move photons
  this->advance_photons(a_dt);

  // Compute the electric field, recompute velocities and diffusion coefficients
  this->solve_poisson();

  // Chemistry kernel.
  this->advance_reaction_network(a_dt);

  // Make superparticles
  //  m_ito->make_superparticles(m_ppc);

  // Prepare next step
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
  
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
