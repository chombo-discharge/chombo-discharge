/*!
  @file   mc_photo.cpp
  @brief  Implementation of mc_photo.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "mc_photo.H"
#include "data_ops.H"
#include "units.H"
#include "poly.H"

#include <time.h>
#include <chrono>

#include <PolyGeom.H>
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
  this->set_bisect_step();
  this->set_domain_bc();
  this->set_random_kappa();
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

  int num_photons, num_outcast;

  // Generate photons
  this->clear(m_photons);
  this->generate_photons(m_photons, a_source, a_dt);                // Generate photons
  this->move_and_absorb_photons(absorbed_photons, m_photons, a_dt); // Move photons
  this->remap_photons(m_photons);                                   // Remap photons
  this->remap_photons(absorbed_photons);                            // Remap photons
  this->deposit_photons(a_state, m_photons);                        // Compute photoionization profile
  
  return true;
}

void mc_photo::clear(EBAMRPhotons& a_photons){
  CH_TIME("mc_photo::clear");
  if(m_verbosity > 5){
    pout() << m_name + "::clear" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    a_photons[lvl]->clear();
  }
}

void mc_photo::set_rng(){
  CH_TIME("mc_photo::set_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::set_rng" << endl;
  }

  // Seed the RNG
  ParmParse pp("mc_photo");
  pp.get("seed", m_seed);
  pp.get("poiss_exp_swap", m_poiss_exp_swap);
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
  m_amr->allocate(m_state,  m_phase, ncomp); 
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

int mc_photo::random_poisson(const Real& a_mean){
  if(a_mean < m_poiss_exp_swap){
    std::poisson_distribution<int> pdist(a_mean);
    return pdist(*m_rng);
  }
  else {
    std::normal_distribution<Real> ndist(a_mean, sqrt(a_mean));
    return (int) round(ndist(*m_rng));
  }
}

int mc_photo::domainbc_map(const int a_dir, const Side::LoHiSide a_side) {
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}

Real mc_photo::random_exponential(const Real a_mean){
  std::exponential_distribution<Real> dist(a_mean);
  return dist(*m_rng);
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

	const Real mean = source(iv,0)*vol*a_dt;//*kappa;
	const int num_phys_photons = random_poisson(mean);
	if(num_phys_photons > 0){

	  const int num_photons = (num_phys_photons <= m_max_photons) ? num_phys_photons : m_max_photons;
	  const Real weight      = (1.0*num_phys_photons)/num_photons;

	  // Generate computational photons
	  for (int i = 0; i < num_photons; i++){
	    const RealVect dir = random_direction();
	    if(m_random_kappa){
	      particles.append(photon(pos, dir*units::s_c0, m_photon_group->get_random_kappa(), weight));
	    }
	    else{
	      particles.append(photon(pos, dir*units::s_c0, m_photon_group->get_kappa(pos), weight));
	    }
	  }
	}
      }

      // Add new particles to data holder
      (*a_particles[lvl])[dit()].addItemsDestructive(particles);
    }
  }

  // Count number of photons
  m_num_photons = this->count_photons(m_photons);
}

void mc_photo::deposit_photons(EBAMRCellData& a_state, const EBAMRPhotons& a_particles){
  CH_TIME("mc_photo::deposit_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_photons" << endl;
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

void mc_photo::set_random_kappa(){
  CH_TIME("mc_photo::set_random_kappa");
  if(m_verbosity > 5){
    pout() << m_name + "::set_random_kappa" << endl;
  }

  std::string str;
  ParmParse pp("mc_photo");
  pp.get("random_kappa", str);

  m_random_kappa = (str == "true") ? true : false;
}

void mc_photo::set_deposition_type(){
  CH_TIME("mc_photo::set_deposition_type");
  if(m_verbosity > 5){
    pout() << m_name + "::set_deposition_type" << endl;
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

void mc_photo::set_domain_bc(){
  CH_TIME("mc_photo::set_domain_bc");
  if(m_verbosity > 5){
    pout() << m_name + "::set_domain_bc" << endl;
  }

  m_domainbc.resize(2*SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
      const int idx = domainbc_map(dir, side);

      ParmParse pp("mc_photo");
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

      std::string sidestr = (side == Side::Lo) ? "_low" : "_high";
      std::string bc_string = "bc_" + str_dir + sidestr;
      std::string type;
      pp.get(bc_string.c_str(), type);
      if(type == "outflow"){
	m_domainbc[idx] = wallbc::outflow;
      }
      else if(type == "symmetry"){
	m_domainbc[idx] = wallbc::symmetry;
      }
      else if(type == "wall"){
	m_domainbc[idx] = wallbc::wall;
      }
      else {
	std::string error = "mc_photo::set_domain_bc - unsupported boundary condition requested: " + bc_string;
	MayDay::Abort(error.c_str());
      }
    }
  }
}

void mc_photo::set_bisect_step(){
  CH_TIME("mc_photo::set_bisect_step");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bisect_step" << endl;
  }

  ParmParse pp("mc_photo");
  pp.get("bisect_step", m_bisect_step);
}

void mc_photo::move_and_absorb_photons(EBAMRPhotons& a_absorbed, EBAMRPhotons& a_original, const Real a_dt){
  CH_TIME("mc_photo::move_and_absorb_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::move_and_absorb_photons" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  const RealVect origin  = m_physdom->get_prob_lo();
  const RealVect prob_hi = m_physdom->get_prob_hi();
  const RefCountedPtr<BaseIF>& impfunc = m_compgeom->get_gas_if();

  // Advance over levels
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<photon>& absorbed  = (*a_absorbed[lvl])[dit()].listItems();
      List<photon>& particles = (*a_original[lvl])[dit()].listItems();
      
      for (ListIterator<photon> lit(particles); lit.ok(); ++lit){
	photon& particle = lit();

	RealVect& pos               = particle.position();
	const RealVect unit_v       = particle.velocity()/units::s_c0;
	const RealVect absorbed_pos = pos + unit_v*random_exponential(particle.kappa());
	const RealVect path         = absorbed_pos - pos;
	const Real path_len         = path.vectorLength();

	// Check if absorbed_pos - pos is smaller than distance to any object
	const Real dist = impfunc->value(pos);

	bool check_eb = false; // Flag for intersection test with eb
	bool check_bc = false; // Flag for intersection test with domain

	// See if we should check for different types of boundary intersections
	if(impfunc->value(pos) < path_len){
	  check_eb = true;
	}
	for (int dir = 0; dir < SpaceDim; dir++){
	  if(absorbed_pos[dir] < origin[dir] || absorbed_pos[dir] > prob_hi[dir]){
	    check_bc = true;
	  }
	}
						  
	// The remaining code does the boundary intersection tests
	if(!check_eb && !check_bc){ // Test necessary
	  pos = absorbed_pos;
	}
	else{

	  // Determine if the particle contacts one of the domain walls, and after how long it contacts
	  bool contact_bc = false;
	  Real contact_s = 2.0;
	  int contact_dir;
	  Side::LoHiSide contact_side;
	  RealVect contact_point;
	  if(check_bc){ // This test does a line-plane intersection test on each side and selects the shortest intersection
	    for (int dir = 0; dir < SpaceDim; dir++){
	      for (SideIterator sit; sit.ok(); ++sit){
		const Side::LoHiSide side = sit();

		const RealVect p0    = (side == Side::Lo) ? origin : prob_hi; // A point on the domain
		const RealVect n0    = sign(side)*RealVect(BASISV(dir));
		const Real norm_path = PolyGeom::dot(n0, path);

		if(norm_path > 0.0){
		  const Real s = PolyGeom::dot(p0-pos, n0)/norm_path;
		  if(s >= 0.0 && s <= 1.0){
		    if(s < contact_s){
		      contact_bc    = true;
		      contact_side  = side;
		      contact_dir   = dir;
		      contact_s     = s;
		      contact_point = pos + contact_s*path;
		    }
		  }
		}
	      }
	    }
	  }
	  else{
	    contact_bc = false;
	  }
	  
	  // Now check the where we contact the embedded boundary
	  if(check_eb){
	    const int nsteps  = ceil(path_len/m_bisect_step);
	    const RealVect dx = (absorbed_pos - pos)/nsteps;

	    // Check each interval
	    bool absorb = false;
	    for (int istep = 0; istep < nsteps; istep++){
	      const Real fa = impfunc->value(pos);
	      const Real fb = impfunc->value(pos+dx);
	      const bool absorb = fa*fb <= 0.0;

	      if(absorb){ 
		// We happen to know that f(pos+dx) > 0.0 and f(pos) < 0.0 so we must now compute the precise location
		// where the absorption happens. For that we use a Brent root finder on the interval [pos, pos+dx].
		const RealVect xcol = poly::brent_root_finder(impfunc, pos, pos+dx);
		pos = xcol;
		absorbed.transfer(lit);
		break;
	      }
	      else{ // Move to next interval
		pos += dx;
	      }
	    }
	  }
	  else if(contact_bc){
	    const int idx = domainbc_map(contact_dir, contact_side);

	    // Check the boundary condition type
	    if(m_domainbc[idx] == wallbc::symmetry){ // Need to bounce back
	      pos += contact_s*path; // Gives point on wall

	      // Wall normal vector

		
	      // Modify velocity vector
	      RealVect& v  = particle.velocity();
	      const RealVect n0 = sign(contact_side)*RealVect(BASISV(contact_dir));

	      v = v - 2.0*PolyGeom::dot(v, n0)*n0; // New direction
	      pos += (1-contact_s)*path_len*v/v.vectorLength();
	    }
	    else if(m_domainbc[idx] == wallbc::wall){
	      pos += contact_s*path;
	      absorbed.transfer(lit);
	    }
	    else if(m_domainbc[idx] == wallbc::outflow){ // Do nothing
	    }
	  }
	}
      }
    }
  }
}

void mc_photo::remap_photons(EBAMRPhotons& a_photons){
  CH_TIME("mc_photo::remap_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::remap_photons" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  const RealVect origin  = m_physdom->get_prob_lo();

  // Upsweep, build outcast lists on each level and transfer from coarser levels to current levels PVR
  for (int lvl = 0; lvl <= finest_level; lvl++){

    // Build outcast list on this level
    a_photons[lvl]->outcast().clear();
    a_photons[lvl]->gatherOutcast();
    a_photons[lvl]->remapOutcast();

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    // Transfer particles from coarser level to this levels PVR
    if(has_coar){
      collectValidParticles(a_photons[lvl]->outcast(),
			    *a_photons[lvl-1],
			    m_pvr[lvl]->mask(),
			    m_amr->get_dx()[lvl]*RealVect::Unit,
			    m_amr->get_ref_rat()[lvl-1],
			    false,
			    origin);
    }
    a_photons[lvl]->gatherOutcast();
    a_photons[lvl]->remapOutcast();
  }

  // Downsweep
  for (int lvl = finest_level; lvl >= 0; lvl--){
    const bool has_coar = lvl > 0;

    if(has_coar){ // Add the current levels outcast list to the coarser level, then rebin the coarser level
      List<photon>& this_outcast = a_photons[lvl]->outcast();
      List<photon>& coar_outcast = a_photons[lvl-1]->outcast();
      coar_outcast.catenate(this_outcast);
      //      this_outcast.clear();

      a_photons[lvl]->gatherOutcast();
      a_photons[lvl]->remapOutcast();

      a_photons[lvl-1]->gatherOutcast();
      a_photons[lvl-1]->remapOutcast();
    }

    a_photons[lvl]->outcast().clear(); // Done with this level, anything that did not get transferred is lost
  }
}

int mc_photo::count_photons(const EBAMRPhotons& a_photons) const {
  CH_TIME("mc_photo::count_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::count_photons" << endl;
  }

  int num_photons = 0;

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    num_photons += a_photons[lvl]->numValid();
  }

  return num_photons;
}

int mc_photo::count_outcast(const EBAMRPhotons& a_photons) const {
  CH_TIME("mc_photo::count_outcast");
  if(m_verbosity > 5){
    pout() << m_name + "::count_outcast" << endl;
  }

  int num_outcast = 0;

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    num_outcast += a_photons[lvl]->numOutcast();
  }

  return num_outcast;
}
