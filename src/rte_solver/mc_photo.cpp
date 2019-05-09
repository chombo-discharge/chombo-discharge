/*!
  @file   mc_photo.cpp
  @brief  Implementation of mc_photo.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "mc_photo.H"
#include "data_ops.H"
#include "units.H"

#include <time.h>
#include <chrono>

#include <EBLevelDataOps.H>
#include <BoxIterator.H>
#include <EBArith.H>
#include <ParmParse.H>
#include <ParticleIO.H>

mc_photo::mc_photo(){
  this->set_verbosity(-1);
  this->set_stationary(false);
  this->set_rng();
}

mc_photo::~mc_photo(){

}

bool mc_photo::advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi){
  data_ops::set_value(a_state, 0.0);

  const RealVect origin  = m_physdom->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();
  const int boxsize      = m_amr->get_max_box_size();
  if(boxsize != m_amr->get_blocking_factor()){
    MayDay::Abort("mc_photo::advance - only constant box sizes are supported for particle methods");
  }

#if 0 // Debug
  if(procID() == 0) std::cout << "advancing mc_photo" << std::endl;
#endif

  EBAMRParticles absorbed_photons;
  m_amr->allocate(absorbed_photons);

  // Generate photons
  this->generate_photons(m_particles, a_source, a_dt);

  const int N = 1;
  for (int i = 0; i < N; i++){
    this->move_and_absorb_photons(absorbed_photons, m_particles, a_dt/N);
  }

  // Deposit absorbed photons onto mesh
#if 1 // Actual code
  this->deposit_photons(a_state, absorbed_photons);
#else // Debug
  this->deposit_photons(a_state, m_particles);
#endif

  data_ops::floor(a_state, 0.0);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    EBLevelDataOps::setCoveredVal(*a_state[lvl], 0.0);
  }

  return true;
}

void mc_photo::set_rng(){
  CH_TIME("mc_photo::set_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::set_rng" << endl;
  }

  // Seed the RNG
  ParmParse pp("mc_photo");
  pp.get("seed", m_seed);
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);

  m_udist01 = new uniform_real_distribution<Real>( 0.0, 1.0);
  m_udist11 = new uniform_real_distribution<Real>(-1.0, 1.0);
}
  
void mc_photo::allocate_internals(){
  CH_TIME("mc_photo::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  const int buffer = 0;
  const int ncomp  = 1;
  m_amr->allocate(m_state,  m_phase, ncomp); // This is the deposited 
  m_amr->allocate(m_source, m_phase, ncomp);
  
  m_amr->allocate(m_particles);
  m_amr->allocate(m_pvr, buffer);
}
  
void mc_photo::cache_state(){
  CH_TIME("mc_photo::cache_state");
  if(m_verbosity > 5){
    pout() << m_name + "::cache_state" << endl;
  }
}

void mc_photo::deallocate_internals(){
  m_amr->deallocate(m_state);
  m_amr->deallocate(m_source);
}

void mc_photo::regrid(const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("mc_photo::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  this->allocate_internals();
}

void mc_photo::compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state){
  CH_TIME("mc_photo::compute_boundary_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_boundary_flux" << endl;
  }
  data_ops::set_value(a_ebflux, 0.0);
}

void mc_photo::compute_domain_flux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_state){
  CH_TIME("mc_photo::compute_domain_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_domain_flux" << endl;
  }
  data_ops::set_value(a_domainflux, 0.0);
}

void mc_photo::compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state){
  MayDay::Abort("mc_photo::compute_flux - I don't think this should ever be called.");
}

void mc_photo::compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state){
  MayDay::Abort("mc_photo::compute_density - I don't think this should ever be called.");
}

void mc_photo::write_plot_file(){
  CH_TIME("mc_photo::write_plot_file");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_file" << endl;
  }
}

int mc_photo::query_ghost() const {
  return 3;
}

RealVect mc_photo::random_direction(){

  Real u1 = 2.0;
  Real u2 = 2.0;
  Real a  = u1*u1 + u2*u2;
  while(a >= 1.0 || a < 1.E-10){
    u1 = (*m_udist11)(*m_rng);
    u2 = (*m_udist11)(*m_rng);
    a  = u1*u1 + u2*u2;
  }

#if CH_SPACEDIM == 2
  RealVect ret = RealVect((u1*u1 - u2*u2), 2*u1*u2)/a;
#else
  const Real b = 2*sqrt(1.0 - a);
  RealVect ret = RealVect(b*u1, b*u2, 1-2*a);
#endif

  return ret;
}

void mc_photo::generate_photons(EBAMRParticles& a_particles, const EBAMRCellData& a_source, const Real a_dt){
  CH_TIME("mc_photo::generate_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::generate_photons" << endl;
  }



  const RealVect origin  = m_physdom->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    const bool has_coar          = lvl > 0;

    if(has_coar) { // If there is a coarser level, remove particles from the overlapping region and put them onto this level
      const int ref_ratio = m_amr->get_ref_rat()[lvl-1];
      collectValidParticles(a_particles[lvl]->outcast(),
      			    *a_particles[lvl-1],
      			    m_pvr[lvl]->mask(),
      			    dx*RealVect::Unit,
      			    ref_ratio,
			    false,
			    origin);
      a_particles[lvl]->outcast().clear(); // Delete particles generated on the coarser level and regenerate them on this level
    }

    // Create new particles on this level using fine-resolution data
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = (*a_source[lvl])[dit()].getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = IntVectSet(box);

      FArrayBox& source = (*a_source[lvl])[dit()].getFArrayBox();

      // Generate new particles in this box
      List<Particle> particles;
      for (VoFIterator vofit(IntVectSet(box), ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const IntVect iv = vof.gridIndex();
	const RealVect pos = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin);
	const Real kappa = ebisbox.volFrac(vof);
	
	const Real mean = source(iv,0)*kappa*vol*a_dt;
	std::poisson_distribution<int> dist(mean);
	const int num_particles = dist(*m_rng);
	for (int i = 0; i < num_particles; i++){
	  const RealVect dir = random_direction();
	  particles.append(Particle(1.0, pos, dir*units::s_c0));
	}
      }

      // Add new particles to data holder
      (*a_particles[lvl])[dit()].addItemsDestructive(particles);
    }
  }
}

void mc_photo::deposit_photons(EBAMRCellData& a_state, const EBAMRParticles& a_particles){
  CH_TIME("mc_photo::generate_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::generate_photons" << endl;
  }

  const RealVect origin  = m_physdom->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      MeshInterp interp(box, dx*RealVect::Unit, origin);
      InterpType type = InterpType::CIC;
      interp.deposit((*a_particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), type);
    }
  }
}

void mc_photo::move_and_absorb_photons(EBAMRParticles& a_absorbed, EBAMRParticles& a_original, const Real a_dt){
  CH_TIME("mc_photo::move_and_absorb_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::move_and_absorb_photons" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const Real dx = m_amr->get_dx()[lvl];

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<Particle>& absorbed  = (*a_absorbed[lvl])[dit()].listItems();
      List<Particle>& particles = (*a_original[lvl])[dit()].listItems();
      
      for (ListIterator<Particle> lit(particles); lit.ok(); ++lit){
	Particle& particle = lit();

	const RealVect vel     = particle.velocity();
	const RealVect old_pos = particle.position();
	const Real kappa       = m_photon_group->get_kappa(old_pos);
	const Real v           = vel.vectorLength();
	const Real max_dx      = Min(0.1/kappa, v*a_dt);
	const int nsteps       = max_dx <= v*a_dt ? ceil(v*a_dt/max_dx) : v*a_dt;
	const Real dtstep      = a_dt/nsteps;

	for (int istep = 0; istep < nsteps; istep++){
	  // We need to determine the maximum allowed step that this photon can move. It is either one grid cell or 0.1*kappa
	  //	  const Real kappa  = 0.5*(m_photon_group->get_kappa(old_pos) + m_photon_group->get_kappa(new_pos));
	  const RealVect dx = vel*dtstep;
	  const Real p      = dx.vectorLength()*kappa;
	  const Real r      = (*m_udist01)(*m_rng);
	  const bool absorb = r < p;

	  if(absorb) {
	    const Real r = (*m_udist01)(*m_rng);
	    particle.position() += r*dx;
	    absorbed.transfer(lit);
	    break;
	  }
	  else{
	    particle.position() += dx;
	  }
	}
      }
    }

    a_original[lvl]->gatherOutcast();
    a_original[lvl]->remapOutcast();

    a_absorbed[lvl]->gatherOutcast();
    a_absorbed[lvl]->remapOutcast();
  }
}
