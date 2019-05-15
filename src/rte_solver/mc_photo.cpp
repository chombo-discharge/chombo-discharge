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

#include <EBAlias.H>
#include <EBLevelDataOps.H>
#include <BoxIterator.H>
#include <EBArith.H>
#include <ParmParse.H>
#include <ParticleIO.H>

mc_photo::mc_photo(){
  this->set_verbosity(-1);
  this->set_stationary(false);
  this->set_rng();
  this->set_pseudophotons();
  this->set_deposition_type();
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

  EBAMRPhotons absorbed_photons;
  m_amr->allocate(absorbed_photons);

  // Generate photons
  this->generate_photons(m_photons, a_source, a_dt);                // Generate photons
  this->move_and_absorb_photons(absorbed_photons, m_photons, a_dt); // Move photons
  this->deposit_photons(a_state, m_photons);                        // Compute photoionization profile

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
  
  m_amr->allocate(m_photons);
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
#if CH_SPACEDIM == 2
  return random_direction2D();
#else
  return random_direction3D();
#endif
}

#if CH_SPACEDIM == 2
RealVect mc_photo::random_direction2D(){
  const Real EPS = 1.E-8;
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = (*m_udist11)(*m_rng);
    x2 = (*m_udist11)(*m_rng);
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
}
#endif

#if CH_SPACEDIM==3
RealVect mc_photo::random_direction3D(){
  const Real EPS = 1.E-8;
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = (*m_udist11)(*m_rng);
    x2 = (*m_udist11)(*m_rng);
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
}
#endif

void mc_photo::generate_photons(EBAMRPhotons& a_particles, const EBAMRCellData& a_source, const Real a_dt){
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
      List<photon> particles;
      for (VoFIterator vofit(IntVectSet(box), ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const IntVect iv    = vof.gridIndex();
	const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin);
	const Real kappa    = ebisbox.volFrac(vof);

	// RNG 
	const Real mean = source(iv,0)*kappa*vol*a_dt;
	std::poisson_distribution<int> dist(mean);

	const int num_phys_photons = dist(*m_rng);
	if(num_phys_photons > 0){

	  const int num_photons = (num_phys_photons <= m_max_photons) ? num_phys_photons : m_max_photons;
	  const Real weight      = (1.0*num_phys_photons)/num_photons;

	  // Generate computational photons
	  for (int i = 0; i < num_photons; i++){
	    const RealVect dir = random_direction();
	    particles.append(photon(pos, dir*units::s_c0, m_photon_group->get_kappa(pos), weight));
	  }
	}
      }

      // Add new particles to data holder
      (*a_particles[lvl])[dit()].addItemsDestructive(particles);
    }
  }
}

void mc_photo::deposit_photons(EBAMRCellData& a_state, const EBAMRPhotons& a_particles){
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
      interp.deposit((*a_particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), m_deposition);
    }

    // Add in the contribution from the ghost cells
    const RefCountedPtr<Copier>& reversecopier = m_amr->get_reverse_copier(m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_state[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);
  }
}

void mc_photo::set_deposition_type(){
  CH_TIME("mc_photo::setdeposition_type");
  if(m_verbosity > 5){
    pout() << m_name + "::setdeposition_type" << endl;
  }

  std::string str;
  ParmParse pp("mc_photo");
  pp.get("deposition", str);

  if(str == "ngp"){
    m_deposition = InterpType::NGP;
  }
  else if(str == "cic"){
    m_deposition = InterpType::CIC;
  }
  else if(str == "tsc"){
    m_deposition = InterpType::TSC;
  }
  else if(str == "w4"){
    m_deposition = InterpType::W4;
  }
  else{
    MayDay::Abort("mc_photo::set_deposition_type - unknown interpolant requested");
  }
}

void mc_photo::set_pseudophotons(){
  CH_TIME("mc_photo::set_pseudophotons");
  if(m_verbosity > 5){
    pout() << m_name + "::set_pseudophotons" << endl;
  }

  ParmParse pp("mc_photo");
  pp.get("max_photons", m_max_photons);
  
  if(m_max_photons <= 0){ // = -1 => no restriction
    m_max_photons = 99999999;
  }
}

void mc_photo::move_and_absorb_photons(EBAMRPhotons& a_absorbed, EBAMRPhotons& a_original, const Real a_dt){
  CH_TIME("mc_photo::move_and_absorb_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::move_and_absorb_photons" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  const RealVect origin  = m_physdom->get_prob_lo();


  // Advance over levels
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<photon>& absorbed  = (*a_absorbed[lvl])[dit()].listItems();
      List<photon>& particles = (*a_original[lvl])[dit()].listItems();
      
      for (ListIterator<photon> lit(particles); lit.ok(); ++lit){
	photon& particle = lit();

	// Position and direction
	RealVect& pos = particle.position();
	const RealVect n = particle.velocity()/units::s_c0;
	
	// Draw the absorption position for this photon
	std::exponential_distribution<Real> dist(m_photon_group->get_kappa(pos));
	const RealVect absorbed_pos = dist(*m_rng)*n;

	// Check if we should absorb
	bool absorb = false;
	const int nsteps = 100;
	const RealVect dx    = absorbed_pos/nsteps;
	for (int istep = 0; istep < nsteps; istep++){
	  pos += dx;

	  const Vector<dielectric>& dielectrics = m_compgeom->get_dielectrics();
	  const RefCountedPtr<BaseIF>& impfunc = dielectrics[0].get_function();
	  absorb = impfunc->value(pos) > 0.0;

	  if(absorb) {
	    const Real r = (*m_udist01)(*m_rng);
	    particle.position() += r*dx;
	    absorbed.transfer(lit);
	    break;
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
