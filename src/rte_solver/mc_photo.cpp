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
  m_name       = "mc_photo";
  m_class_name = "mc_photo";
}

mc_photo::~mc_photo(){

}

bool mc_photo::advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi){
  data_ops::set_value(a_state, 0.0);

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();
  const int boxsize      = m_amr->get_max_box_size();
  if(boxsize != m_amr->get_blocking_factor()){
    MayDay::Abort("mc_photo::advance - only constant box sizes are supported for particle methods");
  }

  EBAMRPhotons absorbed_photons;
  m_amr->allocate(absorbed_photons);

  int num_photons, num_outcast;

  // Generate photons
  //  Real t1 = MPI_Wtime();
  this->clear(m_photons);                                           // Clear internal data
  //  Real t2 = MPI_Wtime();
  this->generate_photons(m_photons, a_source, a_dt);                // Generate photons
  //  Real t3 = MPI_Wtime();
  this->move_and_absorb_photons(absorbed_photons, m_photons, a_dt); // Move photons
  //  Real t4 = MPI_Wtime();
  this->remap_photons(m_photons);                                   // Remap photons
  //  Real t5 = MPI_Wtime();
  this->deposit_photons(a_state, m_photons, m_deposition);          // Compute photoionization profile
  //  Real t6 = MPI_Wtime();

  // Real T = t6-t1;
  // pout() << "total time = " << T << endl
  // 	 << "clear = " << 100.*(t2-t1)/T << endl
  //   	 << "gener = " << 100.*(t3-t2)/T << endl
  // 	 << "move  = " << 100.*(t4-t3)/T << endl
  //   	 << "remap = " << 100.*(t5-t4)/T << endl
  // 	 << "depos = " << 100.*(t6-t5)/T << endl;
  
  return true;
}

void mc_photo::parse_options(){
  CH_TIME("mc_photo::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }
  
  parse_rng();
  parse_pseudophotons();
  parse_photogen();
  parse_source_type();
  parse_deposition();
  parse_bisect_step();
  parse_domain_bc();
  parse_random_kappa();
  parse_pvr_buffer();
  parse_plot_vars();
}

void mc_photo::parse_rng(){
  CH_TIME("mc_photo::parse_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_rng" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_class_name.c_str());
  pp.get("seed", m_seed);
  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);

  m_udist01 = new uniform_real_distribution<Real>( 0.0, 1.0);
  m_udist11 = new uniform_real_distribution<Real>(-1.0, 1.0);
}

void mc_photo::parse_pseudophotons(){
  CH_TIME("mc_photo::parse_pseudophotons");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_pseudophotons" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());
  
  pp.get("max_photons", m_max_photons);
  if(m_max_photons <= 0){ // = -1 => no restriction
    m_max_photons = 99999999;
  }
}

void mc_photo::parse_photogen(){
  CH_TIME("mc_photo::parse_photogen");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_photogen" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("photon_generation", str);

  if(str == "deterministic"){
    m_photogen = photon_generation::deterministic;
  }
  else if(str == "stochastic"){
    m_photogen = photon_generation::stochastic;
  }
  else{
    MayDay::Abort("mc_photo::set_photon_generation - unknown photon generation type requested");
  }
}

void mc_photo::parse_source_type(){
  CH_TIME("mc_photo::parse_source_type");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_source_type" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("source_type", str);

  if(str == "number"){
    m_src_type = source_type::number;
  }
  else if(str == "volume"){
    m_src_type = source_type::per_vol;
  }
  else if(str == "volume_rate"){
    m_src_type = source_type::per_vols;
  }
  else if(str == "rate"){
    m_src_type = source_type::per_s;
  }
  else{
    MayDay::Abort("mc_photo::set_source_type - unknown source type requested");
  }
}

void mc_photo::parse_deposition(){
  CH_TIME("mc_photo::parse_deposition");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_deposition" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("deposition", str);

  m_deposit_numbers = false;
  if(str == "num"){
    m_deposition = DepositionType::NGP;
    m_deposit_numbers = true;
  }
  else if(str == "ngp"){
    m_deposition = DepositionType::NGP;
  }
  else if(str == "cic"){
    m_deposition = DepositionType::CIC;
  }
  else if(str == "tsc"){
    m_deposition = DepositionType::TSC;
  }
  else if(str == "w4"){
    m_deposition = DepositionType::W4;
  }
  else{
    MayDay::Abort("mc_photo::set_deposition_type - unknown interpolant requested");
  }

  pp.get("plot_deposition", str);
  m_plot_numbers = false;
  if(str == "num"){
    m_plot_deposition = DepositionType::NGP;
    m_plot_numbers = true;
  }
  else if(str == "ngp"){
    m_plot_deposition = DepositionType::NGP;
  }
  else if(str == "cic"){
    m_plot_deposition = DepositionType::CIC;
  }
  else if(str == "tsc"){
    m_plot_deposition = DepositionType::TSC;
  }
  else if(str == "w4"){
    m_plot_deposition = DepositionType::W4;
  }
  else{
    MayDay::Abort("mc_photo::set_deposition_type - unknown interpolant requested");
  }
}

void mc_photo::parse_bisect_step(){
  CH_TIME("mc_photo::parse_bisect_step");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_bisect_step" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());
  pp.get("bisect_step", m_bisect_step);
}

void mc_photo::parse_domain_bc(){
  CH_TIME("mc_photo::parse_domain_bc");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_domain_bc" << endl;
  }
  
  m_domainbc.resize(2*SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
      const int idx = domainbc_map(dir, side);

      ParmParse pp(m_class_name.c_str());
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

void mc_photo::parse_random_kappa(){
  CH_TIME("mc_photo::parse_random_kappa");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_random_kappa" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("random_kappa", str);

  m_random_kappa = (str == "true") ? true : false; 
}

void mc_photo::parse_pvr_buffer(){
  CH_TIME("mc_photo::parse_pvr_buffer");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_pvr_buffer" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());
  pp.get("pvr_buffer", m_pvr_buffer);
}

void mc_photo::parse_plot_vars(){
  CH_TIME("mc_photo::parse_plot_vars");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_plot_vars" << endl;
  }

  m_plot_phi = false;
  m_plot_src = false;

  ParmParse pp(m_class_name.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi") m_plot_phi = true;
    else if(str[i] == "src") m_plot_src = true;
  }
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
  
void mc_photo::allocate_internals(){
  CH_TIME("mc_photo::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  const int buffer = 0;
  const int ncomp  = 1;
  m_amr->allocate(m_state,  m_phase, ncomp); 
  m_amr->allocate(m_source, m_phase, ncomp);
  m_amr->allocate(m_crse,   m_phase, ncomp);
		  
  m_amr->allocate(m_photons);
  m_amr->allocate(m_pvr, m_pvr_buffer);
}
  
void mc_photo::cache_state(){
  CH_TIME("mc_photo::cache_state");
  if(m_verbosity > 5){
    pout() << m_name + "::cache_state" << endl;
  }

  // Allocate cache
  m_amr->allocate(m_photocache);

#if 0 // Debug, count number of photons
  long long num_precopy = 0;
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    num_precopy += m_photons[lvl]->numParticles();
  }
#endif

  // Copy photons onto cache
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    collectValidParticles(m_photocache[lvl]->outcast(),
			  *m_photons[lvl],
			  m_pvr[lvl]->mask(),
			  m_amr->get_dx()[lvl]*RealVect::Unit,
			  1,
			  false, 
			  m_amr->get_prob_lo());
    m_photocache[lvl]->remapOutcast();

    // Clear old photons
    m_photons[lvl]->clear();
  }

#if 0 // Debug, count number of photons
  long long num_postcopy = 0;
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    num_postcopy += m_photocache[lvl]->numParticles();
  }

  if(procID() == 0){
    std::cout << "precopy = " << num_precopy << "\t postcopy = " << num_postcopy << std::endl;
  }
#endif
}

void mc_photo::deallocate_internals(){
  m_amr->deallocate(m_state);
  m_amr->deallocate(m_source);
}

void mc_photo::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("mc_photo::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }


  // This reallocates all internal stuff. After this, the only object with any knowledge of the past
  // is m_photocache, which has all the photons in the old data holders.
  this->allocate_internals();


  // Here are the steps:
  // 
  // 1. We are regridding a_lmin and above, so we need to move all photons onto level a_lmin-1. There's
  //    probably a better way to do this
  // 
  // 2. Levels 0 through a_lmin-1 did not change. Copy photons back into those levels. 
  // 
  // 3. Levels a_lmin-1 through a_new_finest_level changed. Remap photons into the correct boxes. 
  //
  // 4. Redeposit photons
  
  const int base_level  = Max(0, a_lmin-1);
  const RealVect origin = m_amr->get_prob_lo();

#if 0 // Debug, count number of photons
  long long num_pretransfer = 0;
  for (int lvl = 0; lvl <= a_old_finest_level; lvl++){
    num_pretransfer += m_photocache[lvl]->numParticles();
  }
#endif
  

  // 1. Move all photons (that are in the cache) onto the coarsest grid level
  List<photon>& base_cache = m_photocache[0]->outcast();
  base_cache.clear();
  for (int lvl = base_level; lvl <= a_old_finest_level; lvl++){
    for (DataIterator dit = m_photocache[lvl]->dataIterator(); dit.ok(); ++dit){
      for (ListIterator<photon> li((*m_photocache[lvl])[dit()].listItems()); li.ok(); ){
	m_photocache[base_level]->outcast().transfer(li);
      }
    }
  }
    

#if 0 // Debug, count number of photons
  long long num_posttransfer = 0;
  for (int lvl = 0; lvl <= a_old_finest_level; lvl++){
    num_posttransfer += m_photocache[lvl]->numParticles();
  }

  if(procID() == 0){
    std::cout << "pre = " << num_pretransfer << "\t post = " << num_posttransfer << std::endl;
  }
#endif

  // 2. Levels 0 through base_level did not change, so copy photons back onto those levels. 
  for (int lvl = 0; lvl <= base_level; lvl++){
    collectValidParticles(*m_photons[lvl],
			  *m_photocache[lvl],
			  m_pvr[lvl]->mask(),
			  m_amr->get_dx()[lvl]*RealVect::Unit,
			  1,
			  false, 
			  origin);

    m_photons[lvl]->gatherOutcast();
    m_photons[lvl]->remapOutcast();
  }

  // 3. Levels above base_level may have changed, so we need to move the particles into the correct boxes now
  for (int lvl = base_level; lvl < a_new_finest_level; lvl++){
    collectValidParticles(m_photons[lvl+1]->outcast(),
			  *m_photons[lvl],
			  m_pvr[lvl+1]->mask(),
			  m_amr->get_dx()[lvl+1]*RealVect::Unit,
			  m_amr->get_ref_rat()[lvl],
			  false,
			  origin);
    m_photons[lvl+1]->remapOutcast();
  }

  // 4. Redeposit photons
  deposit_photons(m_state, m_photons, m_deposition);
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

  MayDay::Abort("mc_photo::write_plot_file - not implemented");
}

void mc_photo::write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("mc_photo::write_checkpoint_level");
  if(m_verbosity > 5){
    pout() << m_name + "::write_checkpoint_level" << endl;
  }

  // Write state vector
  write(a_handle, *m_state[a_level], m_name);

  // Write particles. Must be implemented.
}

void mc_photo::read_checkpoint_level(HDF5Handle& a_handle, const int a_level){
  CH_TIME("mc_photo::read_checkpoint_level");
  if(m_verbosity > 5){
    pout() << m_name + "::read_checkpoint_level" << endl;
  }

  // Read state vector
  read<EBCellFAB>(a_handle, *m_state[a_level], m_name, m_amr->get_grids()[a_level], Interval(0,0), false);

  // Read particles. Should be implemented
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

  const RealVect origin  = m_amr->get_prob_lo();
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
      a_particles[lvl]->outcast().clear(); // Delete particles generated on the coarser level. We'll generate them on this level.
    }

    // Create new particles on this level using fine-resolution data
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = (*a_source[lvl])[dit()].getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = IntVectSet(box);

      const EBCellFAB& source = (*a_source[lvl])[dit()];
      const FArrayBox& srcFAB = source.getFArrayBox();

      Real sum = srcFAB.sum(0);
      
      // Generate new particles in this box
      List<photon> particles;
      if(sum > 0){

	// Regular cells
	for (BoxIterator bit(box); bit.ok(); ++bit){
	  const IntVect iv   = bit();
	  const RealVect pos = origin + (RealVect(iv) + 0.5*RealVect::Unit)*dx;

	  const int num_phys_photons = draw_photons(srcFAB(iv,0), vol, a_dt);
	  
	  if(num_phys_photons > 0){
	    const int num_photons = (num_phys_photons <= m_max_photons) ? num_phys_photons : m_max_photons;
	    const Real weight      = (1.0*num_phys_photons)/num_photons;

	    // Generate computational photons 
	    for (int i = 0; i < num_photons; i++){
	      const RealVect dir = random_direction();
	      if(m_random_kappa){
		particles.append(photon(pos, dir*units::s_c0, m_rte_species->get_random_kappa(), weight));
	      }
	      else{
		particles.append(photon(pos, dir*units::s_c0, m_rte_species->get_kappa(pos), weight));
	      }
	    }
	  }
	}

	// Irregular cells
	for (VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin);
	  const Real kappa    = ebisbox.volFrac(vof);

	  const int num_phys_photons = draw_photons(source(vof,0), vol, a_dt);
	  
	  if(num_phys_photons > 0){
	    const int num_photons = (num_phys_photons <= m_max_photons) ? num_phys_photons : m_max_photons;
	    const Real weight      = (1.0*num_phys_photons)/num_photons;

	    // Generate computational photons 
	    for (int i = 0; i < num_photons; i++){
	      const RealVect dir = random_direction();
	      if(m_random_kappa){
		particles.append(photon(pos, dir*units::s_c0, m_rte_species->get_random_kappa(), weight));
	      }
	      else{
		particles.append(photon(pos, dir*units::s_c0, m_rte_species->get_kappa(pos), weight));
	      }
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

int mc_photo::draw_photons(const Real a_source, const Real a_volume, const Real a_dt){
  CH_TIME("mc_photo::draw_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::draw_photons" << endl;
  }

  int num_photons = 0;

  // Check if we need any type of source term normalization
  Real factor;
  if(m_src_type == source_type::number){
    factor = 1.0;
  }
  else if(m_src_type == source_type::per_vol){
    factor = a_volume;
  }
  else if(m_src_type == source_type::per_vols){
    factor = a_volume*a_dt;
  }
  else if(m_src_type == source_type::per_s){
    factor = a_dt;
  }

  // Draw a number of photons with the desired algorithm
  if(m_photogen == photon_generation::stochastic){
    const Real mean = a_source*factor;
    num_photons = random_poisson(mean);
  }
  else if(m_photogen == photon_generation::deterministic){
    num_photons = round(a_source*factor);
  }
  else{
    MayDay::Abort("mc::draw_photons - unknown generation requested. Aborting...");
  }

  return num_photons;
}

void mc_photo::deposit_photons(EBAMRCellData& a_state, const EBAMRPhotons& a_particles, const DepositionType::Which& a_deposition){
  CH_TIME("mc_photo::deposit_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_photons" << endl;
  }
  
  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  data_ops::set_value(a_state, 0.0);
  data_ops::set_value(m_crse,  0.0);

  DepositionType::Which deposition = a_deposition;

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. If we have a coarser level whose cloud hangs into this level, interpolate the coarser level here first
    if(has_coar){
      RefCountedPtr<EBPWLFineInterp>& interp = m_amr->get_eb_pwl_interp(m_phase)[lvl];
      interp->interpolate(*a_state[lvl], *m_crse[lvl-1], interv);
    }
    
    // 2. Deposit this levels particles and exchange ghost cells
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      EBParticleInterp interp(box, dx*RealVect::Unit, origin);
      interp.deposit((*a_particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), deposition);
    }

    const RefCountedPtr<Copier>& reversecopier = m_amr->get_reverse_copier(m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_state[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 3. If we have a finer level, copy contributions from this level to the temporary holder that is used for
    //    interpolation of "hanging clouds"
    if(has_fine){
      a_state[lvl]->copyTo(*m_crse[lvl]);
    }

    // 4. Add in contributions from photons on finer levels
#if 0
    if(has_fine){
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box box          = dbl.get(dit());
	EBParticleInterp interp(box, dx*RealVect::Unit, origin);
	interp.deposit((*m_joint_photons[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), deposition);
      }
    }
#endif
  }

  // In this case we should deposit numbers and not densities
  if(m_deposit_numbers && a_deposition == DepositionType::NGP){
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      data_ops::scale(*a_state[lvl], pow(m_amr->get_dx()[lvl], SpaceDim));
    }
  }
  else{
    data_ops::kappa_scale(a_state);
  }
}

void mc_photo::move_and_absorb_photons(EBAMRPhotons& a_absorbed, EBAMRPhotons& a_original, const Real a_dt){
  CH_TIME("mc_photo::move_and_absorb_photons");
  if(m_verbosity > 5){
    pout() << m_name + "::move_and_absorb_photons" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  const RealVect origin  = m_amr->get_prob_lo();
  const RealVect prob_hi = m_amr->get_prob_hi();
  const RefCountedPtr<BaseIF>& impfunc = m_compgeom->get_gas_if();

  // Advance over levels
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<photon>& deposit_surf = (*a_absorbed[lvl])[dit()].listItems();
      List<photon>& deposit_bulk = (*a_original[lvl])[dit()].listItems();

      // Copy a list we can iterate through and then clear the list to be filled
      List<photon>  all_photons = deposit_bulk;

      // Clear these, they will be filled.
      deposit_bulk.clear();
      deposit_surf.clear();

      // Iterate through everything
      for (ListIterator<photon> lit(all_photons); lit.ok(); ++lit){
	photon& particle = lit();

	// Draw a randomly absorbed position
	RealVect& pos               = particle.position();
	const RealVect unit_v       = particle.velocity()/units::s_c0;
	const RealVect absorbed_pos = pos + unit_v*random_exponential(particle.kappa());
	const RealVect path         = absorbed_pos - pos;
	const Real path_len         = path.vectorLength();


	// See if we should check for different types of boundary intersections. These are basically
	// cheap initial tests that allow us to skip computations for some particles. 
	bool check_eb = false; // Flag for intersection test with eb
	bool check_bc = false; // Flag for intersection test with domain
	if(impfunc->value(pos) < path_len){ // This tests if we should check for EB intersections
	  check_eb = true;
	}
	for (int dir = 0; dir < SpaceDim; dir++){ // This tests if the photon ended up outside the domain
	  if(absorbed_pos[dir] < origin[dir] || absorbed_pos[dir] > prob_hi[dir]){
	    check_bc = true; // We are guaranteed that the photon crossed a domain boundary
	  }
	}
						  
	// These two loops do the boundary intersection tests. 
	if(!check_eb && !check_bc){ // Impossible for this photon to strike a boundary. Yay!
	  pos = absorbed_pos;
	  deposit_bulk.add(particle);
	}
	else{
	  // Here are the contact points for the EB and domain stuff.
	  // We have x(s) = x0 + s*(x1-x0), s=[0,1].
	  // Take the smallest s to be the intersection points
	  Real s_eb = 2.0;
	  Real s_bc = 2.0;
	  bool contact_bc = false;
	  bool contact_eb = false;
	  int contact_dir;
	  Side::LoHiSide contact_side;
	  
	  // 1. Determine if the particle contacts one of the domain walls, and after how long it contacts
	  if(check_bc){ // This test does a line-plane intersection test on each side and direction
	    for (int dir = 0; dir < SpaceDim; dir++){
	      for (SideIterator sit; sit.ok(); ++sit){
		const Side::LoHiSide side = sit();

		const RealVect p0    = (side == Side::Lo) ? origin : prob_hi; // A point on the domain
		const RealVect n0    = sign(side)*RealVect(BASISV(dir));      // Domain side normal vector
		const Real norm_path = PolyGeom::dot(n0, path);               // Used for test

		if(norm_path > 0.0){ // Moves towards wall
		  const Real s = PolyGeom::dot(p0-pos, n0)/norm_path; // Intersection test
		  if(s >= 0.0 && s <= 1.0){
		    if(s < s_bc){
		      contact_bc    = true;
		      contact_side  = side;
		      contact_dir   = dir;
		      s_bc          = s;
		    }
		  }
		}
	      }
	    }
	  }
	  
	  // 2. Now check the where we contact the embedded boundary. We do a step-wise increment and
	  //    check where the photon cross the boundary. If it does, find the position and flag the photon. 
	  if(check_eb){
	    const int nsteps  = ceil(path_len/m_bisect_step);
	    const RealVect dx = (absorbed_pos - pos)/nsteps;
	    RealVect cur_pos  = pos;
	    
	    // Check each interval
	    for (int istep = 0; istep < nsteps; istep++){
	      const Real fa = impfunc->value(cur_pos);
	      const Real fb = impfunc->value(cur_pos+dx);

	      if(fa*fb <= 0.0){ 
		// We happen to know that f(pos+dx) > 0.0 and f(pos) < 0.0 so we must now compute the precise location
		// where the photon crossed the EB. For that we use a Brent root finder on the interval [pos, pos+dx].
		const RealVect xcol = poly::brent_root_finder(impfunc, cur_pos, cur_pos+dx);
		s_eb = (xcol - pos).vectorLength()/path.vectorLength();
		contact_eb = true;
		break;
	      }
	      else{ // Move to next interval
		cur_pos += dx;
	      }
	    }
	  }

	  // 3. If we don't contact either boundary, absorb the particle in the bulk. Otherwise, absorb in on a surface
	  //    somewhere
	  if(!contact_eb && !contact_bc){
	    pos = absorbed_pos;
	    deposit_bulk.add(particle);
	  }
	  else{
	    pos += Min(s_eb, s_bc)*path; // Move the photon to the absorption point. This is either on the boundary
	    deposit_surf.add(particle);

#if 0 // This code needs to be rewritten, I don't think it does 
	    // Domain boundaries may be bounce-back. 
	    if(contact_bc && s_bc < s_eb){ // Photon was absorbed on the domain 
	      const int idx = domainbc_map(contact_dir, contact_side);

	      // Check the boundary condition type
	      if(m_domainbc[idx] == wallbc::symmetry){ // Need to bounce back

		// Modify velocity vector
		RealVect& v  = particle.velocity();
		const RealVect n0 = sign(contact_side)*RealVect(BASISV(contact_dir));

		v = v - 2.0*PolyGeom::dot(v, n0)*n0; // New direction
		pos += (1-s_bc)*path_len*v/v.vectorLength();
	      }
	      else if(m_domainbc[idx] == wallbc::wall){
		absorbed.transfer(lit);
	      }
	      else if(m_domainbc[idx] == wallbc::outflow){ // Do nothing
	      }
	    }
	    else if(contact_eb && s_eb < s_bc) { // Photon was absorbed on the EB
	      absorbed.transfer(lit);
	    }
#endif
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
  const RealVect origin  = m_amr->get_prob_lo();

  // 1. Gather everything on the coarsest level
  List<photon>& coarsest_outcast = a_photons[0]->outcast();
  a_photons[0]->gatherOutcast();
  for (int lvl = 1; lvl <= finest_level; lvl++){
    a_photons[lvl]->outcast().clear();
    a_photons[lvl]->gatherOutcast();
    coarsest_outcast.catenate(a_photons[lvl]->outcast());
  }

  // Remap coarsest level
  a_photons[0]->remapOutcast();
  coarsest_outcast.clear();
  
  for (int lvl = 1; lvl <= finest_level; lvl++){

    // 1. Collect coarser level particles into this levels PVR 
    collectValidParticles(a_photons[lvl]->outcast(),
			  *a_photons[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_amr->get_dx()[lvl]*RealVect::Unit,
			  m_amr->get_ref_rat()[lvl-1],
			  false, 
			  origin);
    a_photons[lvl]->remapOutcast();

    // 2. There may be particles that remained on this levels DBL but may not belong to the PVR. Move those 
    //    particles one level down and remap that level once more
    if(m_pvr_buffer > 0){
      collectValidParticles(a_photons[lvl-1]->outcast(),
			    *a_photons[lvl],
			    m_pvr[lvl]->mask(),
			    m_amr->get_dx()[lvl]*RealVect::Unit,
			    1,
			    true, 
			    origin);

      a_photons[lvl-1]->remapOutcast();
    }

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

void mc_photo::binmapPhotons(std::map<IntVect, joint_photon, CompIntVect>& a_mip,
			     const List<photon>&                           a_photons,
			     const RealVect&                               a_meshSpacing,
			     const RealVect&                               a_origin){

  for (ListIterator<photon> li(a_photons); li; ++li){
    const photon& p = li();
    const IntVect bin = locateBin(p.position(), a_meshSpacing, a_origin);

    std::map<IntVect, joint_photon, CompIntVect>::iterator it = a_mip.find(bin);

    if(it == a_mip.end()){ // Create a new joint_photon
      a_mip[bin] = joint_photon(p.mass(), p.position(), 1);
    }
    else{ // Increment to the one that is already here
      it->second.add_photon(&p);
    }
  }
}

void mc_photo::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("mc_photo::write_plot_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_data" << endl;
  }

  if(m_plot_phi) {
    const int ncomp = 1;
    EBAMRCellData dummy;
    m_amr->allocate(dummy, m_phase, ncomp);
    deposit_photons(dummy, m_photons, m_plot_deposition);
    write_data(a_output, a_comp, dummy,  true);
  }
  if(m_plot_src) write_data(a_output, a_comp, m_source, false);
}
