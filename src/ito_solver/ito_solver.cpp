/*!
  @file   ito_solver.cpp
  @brief  Declaration of an abstract class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
*/

#include "simple_ito_particle.H"
#include "ito_solver.H"
#include "data_ops.H"
#include "EBParticleInterp.H"
#include "units.H"
#include "EBGhostCloud.H"
#include "ito_layout.H"


#include "particle_ops.H"

#include <EBArith.H>
#include <ParmParse.H>
#include <EBAlias.H>
#include <BaseEBCellFactory.H>
#include <ParticleIO.H>

#include <chrono>

#define ITO_DEBUG 0

ito_solver::ito_solver(){
  m_name       = "ito_solver";
  m_class_name = "ito_solver";

  m_mobility_interp = mobility_interp::mobility;
}

ito_solver::~ito_solver(){

}

std::string ito_solver::get_name(){
  return m_name;
}

const std::string ito_solver::get_realm() const{
  return m_realm;
}

void ito_solver::set_realm(const std::string a_realm){
  m_realm = a_realm;

  // This is for later, in case we want to move averaging/interpolating onto a different realm. 
  m_fluid_realm = m_realm;
}

RefCountedPtr<ito_species>& ito_solver::get_species(){
  return m_species;
}

void ito_solver::parse_options(){
  CH_TIME("ito_solver::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }

  this->parse_superparticles();
  this->parse_rng();
  this->parse_plot_vars();
  this->parse_deposition();
  this->parse_bisect_step();
  this->parse_pvr_buffer();
  this->parse_diffusion_hop();
  this->parse_redistribution();
  this->parse_conservation();
  this->parse_checkpointing();
}

void ito_solver::parse_superparticles(){
  CH_TIME("ito_solver::parse_superparticles");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_superparticles" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_class_name.c_str());
  pp.get("kd_direction", m_kd_direction);

  m_kd_direction = min(m_kd_direction, SpaceDim-1);
}

void ito_solver::parse_rng(){
  CH_TIME("ito_solver::parse_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_rng" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_class_name.c_str());
  pp.get("seed",       m_seed_rng);
  pp.get("normal_max", m_normal_max);
  if(m_seed_rng < 0) { // Random seed if input < 0
    m_seed_rng = std::chrono::system_clock::now().time_since_epoch().count();
  }
  
  m_rng     = RAN::mt19937_64(m_seed_rng);
  m_udist01 = RAN::uniform_real_distribution<Real>( 0.0, 1.0);
  m_udist11 = RAN::uniform_real_distribution<Real>(-1.0, 1.0);
  m_gauss01 = RAN::normal_distribution<Real>(0.0, 1.0);
  m_udist0d = RAN::uniform_int_distribution<int>(0, SpaceDim-1);
}

void ito_solver::parse_plot_vars(){
  CH_TIME("mc_photo::parse_plot_vars");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_plot_vars" << endl;
  }

  m_plot_phi = false;
  m_plot_vel = false;
  m_plot_dco = false;
  m_plot_particles        = false;
  m_plot_eb_particles     = false;
  m_plot_domain_particles = false;
  m_plot_source_particles = false;
  m_plot_energy_density   = false;
  m_plot_average_energy   = false;

  ParmParse pp(m_class_name.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi")            m_plot_phi              = true;
    else if(str[i] == "vel")            m_plot_vel              = true;
    else if(str[i] == "dco")            m_plot_dco              = true;
    else if(str[i] == "part")           m_plot_particles        = true;
    else if(str[i] == "eb_part")        m_plot_eb_particles     = true;
    else if(str[i] == "dom_part")       m_plot_domain_particles = true;
    else if(str[i] == "src_part")       m_plot_source_particles = true;
    else if(str[i] == "energy_density") m_plot_energy_density   = true;
    else if(str[i] == "average_energy") m_plot_average_energy   = true;
  }
}

void ito_solver::parse_deposition(){
  CH_TIME("ito_solver::parse_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_rng" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  std::string str;

  // Deposition for particle-mesh operations
  pp.get("deposition", str);
  if(str == "ngp"){
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
    MayDay::Abort("ito_solver::parse_deposition - unknown interpolant requested");
  }

  // Deposition for plotting only
  pp.get("plot_deposition", str);

  if(str == "ngp"){
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
    MayDay::Abort("ito_solver::parse_deposition - unknown interpolant requested");
  }

  // Mobility interpolation.
  pp.get("mobility_interp",str);

  if(str == "mobility"){
    m_mobility_interp = mobility_interp::mobility;
  }
  else if(str == "velocity"){
    m_mobility_interp = mobility_interp::velocity;
  }
  else{
    MayDay::Abort("ito_solver::parse_deposition - unknown interpolation method for mobility");
  }

  pp.get("irr_ngp_deposition", m_irreg_ngp_deposition);
  pp.get("irr_ngp_interp",     m_irreg_ngp_interpolation);
}

void ito_solver::parse_bisect_step(){
  CH_TIME("ito_solver::parse_bisect_step");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_bisect_step" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  pp.get("bisect_step", m_bisect_step);
}

void ito_solver::parse_pvr_buffer(){
  CH_TIME("ito_solver::parse_pvr_buffer");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_pvr_buffer" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  pp.get("pvr_buffer",  m_pvr_buffer);
  pp.get("halo_buffer", m_halo_buffer);

  std::string str;
  pp.get("halo_deposition", str);
  if(str == "ngp"){
    m_ngp_halo = true;
  }
  else if(str == "native"){
    m_ngp_halo = false;
  }
  else{
    MayDay::Abort("ito_solver::parse_pvr_buffer - unknown argument to 'halo_deposition'");
  }
}

void ito_solver::parse_diffusion_hop(){
  CH_TIME("ito_solver::parse_diffusion_hop");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_diffusion_hop" << endl;
  }

  ParmParse pp(m_class_name.c_str());

  pp.get("max_diffusion_hop", m_max_diffusion_hop);
}

void ito_solver::parse_redistribution(){
  CH_TIME("ito_solver::parse_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_redistribution" << endl;
  }

  ParmParse pp(m_class_name.c_str());

  pp.get("redistribute", m_redistribute);
}

void ito_solver::parse_conservation(){
  CH_TIME("ito_solver::parse_conservation");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_conservation" << endl;
  }

  ParmParse pp(m_class_name.c_str());

  pp.get("blend_conservation", m_blend_conservation);
}

void ito_solver::parse_checkpointing(){
  CH_TIME("ito_solver::parse_checkpointing");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_checkpointing" << endl;
  }

  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("checkpointing", str);
  pp.get("ppc_restart", m_ppc_restart);
  if(str == "particles"){
    m_checkpointing = which_checkpoint::particles;
  }
  else if(str == "numbers"){
    m_checkpointing = which_checkpoint::numbers;
  }
  else{
    MayDay::Abort("ito_solver::parse_checkpointing - unknown checkpointing method requested");
  }
}

Vector<std::string> ito_solver::get_plotvar_names() const {
  CH_TIME("ito_solver::get_plotvar_names");
  if(m_verbosity > 5){
    pout() << m_name + "::get_plotvar_names" << endl;
  }

  Vector<std::string> names(0);
  if(m_plot_phi) {
    names.push_back(m_name + " phi");
  }
  if(m_plot_dco && m_diffusive){
    names.push_back(m_name + " diffusion_coefficient");
  }
  if(m_plot_vel && m_mobile){
    names.push_back("x-Velocity " + m_name);
    names.push_back("y-Velocity " + m_name);
    if(SpaceDim == 3){
      names.push_back("z-Velocity " + m_name);
    }
  }
  if(m_plot_particles)         names.push_back(m_name + " particles");
  if(m_plot_eb_particles)      names.push_back(m_name + " eb_particles");
  if(m_plot_domain_particles)  names.push_back(m_name + " domain_particles");
  if(m_plot_source_particles)  names.push_back(m_name + " source_particles");
  if(m_plot_covered_particles) names.push_back(m_name + " covered_particles");
  if(m_plot_energy_density)    names.push_back(m_name + " energy * phi");
  if(m_plot_average_energy)    names.push_back(m_name + " average_energy");

  return names;
}

int ito_solver::get_num_plotvars() const {
  CH_TIME("ito_solver::get_num_plotvars");
  if(m_verbosity > 5){
    pout() << m_name + "::get_num_plotvars" << endl;
  }

  int num_plotvars = 0;
  
  if(m_plot_phi)                num_plotvars += 1;
  if(m_plot_dco && m_diffusive) num_plotvars += 1;
  if(m_plot_vel && m_mobile)    num_plotvars += SpaceDim;
  if(m_plot_particles)          num_plotvars = num_plotvars + 1;
  if(m_plot_eb_particles)       num_plotvars = num_plotvars + 1;
  if(m_plot_domain_particles)   num_plotvars = num_plotvars + 1;
  if(m_plot_source_particles)   num_plotvars = num_plotvars + 1;
  if(m_plot_covered_particles)  num_plotvars = num_plotvars + 1;
  if(m_plot_energy_density)     num_plotvars += 1;
  if(m_plot_average_energy)     num_plotvars += 1;

  return num_plotvars;
}

int ito_solver::get_pvr_buffer() const {
  CH_TIME("ito_solver::get_pvr_buffer");
  if(m_verbosity > 5){
    pout() << m_name + "::get_pvr_buffer" << endl;
  }

  return m_pvr_buffer;
}

int ito_solver::get_halo_buffer() const {
  CH_TIME("ito_solver::get_halo_buffer");
  if(m_verbosity > 5){
    pout() << m_name + "::get_halo_buffer" << endl;
  }

  return m_halo_buffer;
}

void ito_solver::set_pvr_buffer(const int a_buffer) {
  CH_TIME("ito_solver::set_pvr_buffer");
  if(m_verbosity > 5){
    pout() << m_name + "::set_pvr_buffer" << endl;
  }

  m_pvr_buffer = a_buffer;
}

void ito_solver::set_halo_buffer(const int a_buffer)  {
  CH_TIME("ito_solver::set_halo_buffer");
  if(m_verbosity > 5){
    pout() << m_name + "::set_halo_buffer" << endl;
  }

  m_halo_buffer = a_buffer;
}


size_t ito_solver::get_num_particles(const which_container a_container, const bool a_local) const{
  CH_TIME("ito_solver::get_num_particles(string, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::get_num_particles(string, bool)" << endl;
  }

  size_t N;

  const particle_container<ito_particle>& particles = m_particle_containers.at(a_container);

  if(a_local){
    N = particles.get_num_valid_local();
  }
  else{
    N = particles.get_num_valid_global();
  }

  return N;
}

void ito_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("ito_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  
  m_compgeom = a_compgeom;
}

void ito_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("ito_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }

  m_amr = a_amr;
}

void ito_solver::register_operators(){
  CH_TIME("ito_solver::register_operators");
  if(m_verbosity > 5){
    pout() << m_name + "::register_operators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("cdr_solver::register_operators - need to set amr_mesh!");
  }
  else{
    m_amr->register_operator(s_eb_coar_ave,     m_realm, m_phase);
    m_amr->register_operator(s_eb_fill_patch,   m_realm, m_phase);
    m_amr->register_operator(s_eb_mg_interp,    m_realm, m_phase);
    m_amr->register_operator(s_eb_copier,       m_realm, m_phase);
    m_amr->register_operator(s_eb_ghostcloud,   m_realm, m_phase);
    m_amr->register_operator(s_eb_noncons_div,  m_realm, m_phase);
    
    if(m_redistribute)    m_amr->register_operator(s_eb_redist,  m_realm, m_phase);
    if(m_halo_buffer > 0) m_amr->register_mask(s_particle_halo, m_halo_buffer, m_realm);
  }
}

void ito_solver::set_phase(phase::which_phase a_phase){
  CH_TIME("ito_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void ito_solver::set_verbosity(const int a_verbosity){
  CH_TIME("ito_solver::set_verbosity");
  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
}

void ito_solver::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("ito_solver::set_time");
  if(m_verbosity > 5){
    pout() << m_name + "::set_time" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void ito_solver::initial_data(){
  CH_TIME("ito_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }

  particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  particles.clear_particles();

  // Add particles, remove the ones that are inside the EB, and then depsit
  particles.add_particles(m_species->get_initial_particles());
  this->remove_covered_particles(particles, EB_representation::implicit_function, 0.0);
  this->deposit_particles(m_state, particles, m_deposition);
}

void ito_solver::compute_loads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level){
  CH_TIME("ito_solver::compute_loads");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_loads" << endl;
  }

  const particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  a_loads.resize(a_dbl.size(),0);
  
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit){
    const int numPart = particles[a_level][dit()].numItems();
    a_loads[dit().intCode()] = numPart;
  }

  // Gather loads globally
#ifdef CH_MPI
  int count = a_loads.size();
  Vector<long int> tmp(count);
  MPI_Allreduce(&(a_loads[0]),&(tmp[0]), count, MPI_LONG, MPI_SUM, Chombo_MPI::comm);
  a_loads = tmp;
#endif
}

void ito_solver::remove_covered_particles(const EB_representation a_representation, const Real a_tol){
  CH_TIME("ito_solver::remove_covered_particles(EB_representation, tolerance)");
  if(m_verbosity > 5){
    pout() << m_name + "::remove_covered_particles(EB_representation, tolerance)" << endl;
  }

  this->remove_covered_particles(which_container::bulk, a_representation, a_tol);
}


void ito_solver::remove_covered_particles(const which_container a_container, const EB_representation a_representation, const Real a_tol){
  CH_TIME("ito_solver::remove_covered_particles(container, EB_representation, tolerance)");
  if(m_verbosity > 5){
    pout() << m_name + "::remove_covered_particles(container, EB_representation, tolerance)" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(a_container);
  
  this->remove_covered_particles(particles, a_representation, a_tol);
}

void ito_solver::remove_covered_particles(particle_container<ito_particle>& a_particles, const EB_representation a_representation, const Real a_tol){
  CH_TIME("ito_solver::remove_covered_particles(particles, EB_representation)");
  if(m_verbosity > 5){
    pout() << m_name + "::remove_covered_particles(particles, EB_representation)" << endl;
  }

  switch(a_representation){
  case EB_representation::implicit_function:
    this->remove_covered_particles_if(a_particles, a_tol);
    break;
  case EB_representation::discrete:
    this->remove_covered_particles_discrete(a_particles);
    break;
  case EB_representation::voxel:
    this->remove_covered_particles_voxels(a_particles);
    break;
  default:
    MayDay::Abort("ito_solver::remove_covered_particles - unsupported EB representation requested");
  }
}

void ito_solver::remove_covered_particles_if(particle_container<ito_particle>& a_particles, const Real a_tol){
  CH_TIME("ito_solver::remove_covered_particles_if(particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::remove_covered_particles_if(particles)" << endl;
  }

  const RefCountedPtr<BaseIF>& func = (m_phase == phase::gas) ? m_compgeom->get_gas_if() : m_compgeom->get_sol_if();

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    const Real tol               = a_tol*dx;

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      // Check if particles are outside the implicit function. 
      for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	ito_particle& p = lit();

	const Real f = func->value(p.position());
	if(f > tol) particles.remove(lit);
      }
    }
  }
}


void ito_solver::transfer_covered_particles_if(particle_container<ito_particle>& a_src, particle_container<ito_particle>& a_dst, const Real a_tol){
  CH_TIME("ito_solver::transfer_covered_particles_if(container, container, tolerance)");
  if(m_verbosity > 5){
    pout() << m_name + "::transfer_covered_particles_if(container, container, tolerance)" << endl;
  }

  const RefCountedPtr<BaseIF>& func = (m_phase == phase::gas) ? m_compgeom->get_gas_if() : m_compgeom->get_sol_if();

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    const Real tol               = a_tol*dx;

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ito_particle>& src = a_src[lvl][dit()].listItems();
      List<ito_particle>& dst = a_dst[lvl][dit()].listItems();

      // Check if particles are outside the implicit function. 
      for (ListIterator<ito_particle> lit(src); lit.ok(); ++lit){
	ito_particle& p = lit();

	const Real f = func->value(p.position());
	if(f > tol) {
	  dst.transfer(lit);
	}
      }
    }
  }
}

void ito_solver::remove_covered_particles_discrete(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::remove_covered_particles_discrete(particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::remove_covered_particles_discrete(particles)" << endl;
  }

  const RealVect prob_lo = m_amr->get_prob_lo();

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      if(ebisbox.isAllCovered()){ // Box is all covered, remove everything. 
	particles.clear();
      }
      else if(ebisbox.isAllRegular()){ // Do nothing
      }
      else{
	for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	  ito_particle& p = lit();

	  
	  const RealVect rv  = (p.position() - prob_lo)/dx;
	  const IntVect iv   = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));


	  if(ebisbox.isCovered(iv)){
	    particles.remove(lit);
	  }
	  else if(ebisbox.isIrregular(iv)){
	    const VolIndex vof = VolIndex(iv, 0);
	    const RealVect xc  = ebisbox.bndryCentroid(vof);
	    const RealVect nc  = ebisbox.normal(vof);

	    const Real proj    = PolyGeom::dot(p.position() - xc, nc);

	    if(proj < 0.0){
	      particles.remove(lit);
	    }
	  }
	}
      }
    }
  }
}

void ito_solver::remove_covered_particles_voxels(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::remove_covered_particles_voxels(particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::remove_covered_particles_voxels(particles)" << endl;
  }

  const RealVect prob_lo = m_amr->get_prob_lo();

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      if(ebisbox.isAllCovered()){ // Box is all covered, remove everything. 
	particles.clear();
      }
      else if(ebisbox.isAllRegular()){ // Do nothing
      }
      else{
	for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	  ito_particle& p = lit();
	  
	  const RealVect rv  = (p.position() - prob_lo)/dx;
	  const IntVect iv   = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

	  if(ebisbox.isCovered(iv)){
	    particles.remove(lit);
	  }
	}
      }
    }
  }
}

void ito_solver::transfer_covered_particles(const EB_representation a_representation, const Real a_tol){
  CH_TIME("ito_solver::transfer_covered_particles(EB_representation, tol)");
  if(m_verbosity > 5){
    pout() << m_name + "::transfer_covered_particles(EB_representation, tol)" << endl;
  }

  this->transfer_covered_particles(which_container::bulk, which_container::covered, a_representation, a_tol);
}

void ito_solver::transfer_covered_particles(const which_container a_containerFrom, const which_container a_containerTo, const EB_representation a_representation, const Real a_tol){
  CH_TIME("ito_solver::transfer_covered_particles(EB_representation, string, string, tol)");
  if(m_verbosity > 5){
    pout() << m_name + "::transfer_covered_particles(EB_representation, string, string, tol)" << endl;
  }

  particle_container<ito_particle>& src = this->get_particles(a_containerFrom);
  particle_container<ito_particle>& dst = this->get_particles(a_containerTo);

  this->transfer_covered_particles(src, dst, a_representation, a_tol);
}

void ito_solver::transfer_covered_particles(particle_container<ito_particle>& a_containerFrom,
					    particle_container<ito_particle>& a_containerTo,
					    const EB_representation           a_representation,
					    const Real                        a_tol){
  CH_TIME("ito_solver::transfer_covered_particles(EB_representation, container, container, tol)");
  if(m_verbosity > 5){
    pout() << m_name + "::transfer_covered_particles(EB_representation, container, container, tol)" << endl;
  }

  switch(a_representation){
  case EB_representation::implicit_function:
    this->transfer_covered_particles_if(a_containerFrom, a_containerTo, a_tol);
    break;
  default:
    MayDay::Abort("ito_solver::intersect_particles - unsupported EB representation requested");
  }
}

void ito_solver::intersect_particles(const EB_representation a_representation, const bool a_delete){
  CH_TIME("ito_solver::intersect_particles(EB_representation, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::intersect_particles(EB_representation, bool)" << endl;
  }

  this->intersect_particles(which_container::bulk, which_container::eb, which_container::domain, a_representation, a_delete);
}

void ito_solver::intersect_particles(const which_container       a_particles,
				     const which_container       a_eb_particles,
				     const which_container       a_dom_particles,
				     const EB_representation a_representation,				     
				     const bool              a_delete){
  CH_TIME("ito_solver::intersect_particles(string, string, string, bool, EB_representation)");
  if(m_verbosity > 5){
    pout() << m_name + "::intersect_particles(string, string, string, bool, EB_representation)" << endl;
  }

  particle_container<ito_particle>& particles     = this->get_particles(a_particles);
  particle_container<ito_particle>& eb_particles  = this->get_particles(a_eb_particles);
  particle_container<ito_particle>& dom_particles = this->get_particles(a_dom_particles);

  this->intersect_particles(particles, eb_particles, dom_particles, a_representation, a_delete);
}


void ito_solver::intersect_particles(particle_container<ito_particle>& a_particles,
				     particle_container<ito_particle>& a_eb_particles,
				     particle_container<ito_particle>& a_dom_particles,
				     const EB_representation           a_representation,
				     const bool                        a_delete){
  CH_TIME("ito_solver::intersect_particles(container, container, container, EB_representation, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::intersect_particles(container, container, container, EB_representation, bool)" << endl;
  }

  switch(a_representation){
  case EB_representation::implicit_function:
    this->intersect_particles_if(a_particles, a_eb_particles, a_dom_particles, a_delete);
    break;
  default:
    MayDay::Abort("ito_solver::intersect_particles - unsupported EB representation requested");
  }
}

void ito_solver::intersect_particles_if(particle_container<ito_particle>& a_particles,
					particle_container<ito_particle>& a_eb_particles,
					particle_container<ito_particle>& a_domain_particles,
					const bool                        a_delete){
  CH_TIME("ito_solver::intersect_particles_if(container, container, container, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::intersect_particles_if(container, container, container, bool)" << endl;
  }

  const RealVect prob_lo = m_amr->get_prob_lo();
  const RealVect prob_hi = m_amr->get_prob_hi();
  const Real     SAFETY  = 1.E-6;

  // This is the implicit function used for intersection tests
  const RefCountedPtr<BaseIF>& impfunc = m_amr->get_baseif(m_phase);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    for (DataIterator dit = m_amr->get_grids(m_realm)[lvl]; dit.ok(); ++dit){

      const Real     dx      = m_amr->get_dx()[lvl];
      const EBISBox& ebisbox = m_amr->get_ebisl(m_realm, m_phase)[lvl][dit()];

      if(!ebisbox.isAllRegular()){
	List<ito_particle>& particles    = a_particles[lvl][dit()].listItems();
	List<ito_particle>& ebParticles  = a_eb_particles[lvl][dit()].listItems();
	List<ito_particle>& domParticles = a_domain_particles[lvl][dit()].listItems();

	ebParticles.clear();
	domParticles.clear();

	for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	  ito_particle& p = lit();

	  const RealVect newPos  = p.position();
	  const RealVect oldPos  = p.oldPosition();
	  const RealVect path    = newPos - oldPos;
	  const Real     pathLen = path.vectorLength();

	  // Check if we should check of different types of boundary intersections. These are checp initial tests that allow
	  // us to skip intersection tests for some photons.
	  bool checkEB  = false;
	  bool checkDom = false;

	  if(impfunc->value(oldPos) < pathLen){ // Checks distance to EB. 
	    checkEB = true;
	  }
	  for (int dir = 0; dir < SpaceDim; dir++){
	    if(newPos[dir] < prob_lo[dir] || newPos[dir] > prob_hi[dir]){ // Checks if we crossed a domain boundary. 
	      checkDom = true; 
	    }
	  }

	  // Must do intersection test on at least one of these. 
	  if(checkEB || checkDom){ 
	    Real dom_s = 1.E99;
	    Real eb_s  = 1.E99;

	    bool contact_domain = false;
	    bool contact_eb     = false;
	      
	    if(checkDom) contact_domain = particle_ops::domain_bc_intersection(oldPos, newPos, path, prob_lo, prob_hi, dom_s);
	    if(checkEB)  contact_eb     = particle_ops::eb_bc_intersection(impfunc, oldPos, newPos, pathLen, dx, eb_s);
	  
	    if(contact_eb || contact_domain){ // Particle trajectory crossed something. 
	      if(eb_s < dom_s){ // It was the EB first. 
		p.position() = oldPos + eb_s*path;
		if(a_delete){
		  ebParticles.transfer(lit);
		}
		else{
		  ebParticles.add(lit());
		}
	      }
	      else{ // It was the domain boundary. 
		p.position() = oldPos + Max(0.0,dom_s-SAFETY)*path;
		if(a_delete){
		  domParticles.transfer(lit);
		}
		else{
		  domParticles.add(lit());
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void ito_solver::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("ito_solver::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  // Reallocate mesh data
  m_amr->reallocate(m_state,         m_phase, a_lmin);
  m_amr->reallocate(m_scratch,       m_phase, a_lmin);
  m_amr->reallocate(m_mobility_func, m_phase, a_lmin);
  m_amr->reallocate(m_depositionNC,  m_phase, a_lmin);
  m_amr->reallocate(m_massDiff,      m_phase, a_lmin);
  
  // Only allocate memory if we actually have a mobile solver
  if(m_mobile){
    m_amr->reallocate(m_velo_func, m_phase, a_lmin);
  }
  else{ 
    m_amr->allocate_ptr(m_velo_func);
  }

  // Only allocate memory if we actually a diffusion solver
  if(m_diffusive){
    m_amr->reallocate(m_diffco_cell, m_phase, a_lmin);
  }
  else{
    m_amr->allocate_ptr(m_diffco_cell);
  }

  // Particle data regrids
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids(m_realm);
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const Vector<int>& ref_rat             = m_amr->get_ref_rat();

  for (auto& container : m_particle_containers){
    particle_container<ito_particle>& particles = container.second;
    particles.regrid(grids, domains, dx, ref_rat, a_lmin, a_new_finest_level);
  }
}

void ito_solver::set_species(RefCountedPtr<ito_species> a_species){
  CH_TIME("ito_solver::set_species");
  if(m_verbosity > 5){
    pout() << m_name + "::set_species" << endl;
  }
  
  m_species   = a_species;
  m_name      = a_species->get_name();
  m_diffusive = m_species->is_diffusive();
  m_mobile    = m_species->is_mobile();
}

void ito_solver::allocate_internals(){
  CH_TIME("ito_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }
  
  const int ncomp = 1;

  m_amr->allocate(m_state,         m_realm, m_phase, ncomp);
  m_amr->allocate(m_scratch,       m_realm, m_phase, ncomp);
  m_amr->allocate(m_mobility_func, m_realm, m_phase, ncomp);
  m_amr->allocate(m_depositionNC,  m_realm, m_phase, ncomp);
  m_amr->allocate(m_massDiff,      m_realm, m_phase, ncomp);

  // Only allocate memory for velocity if we actually have a mobile solver
  if(m_mobile){
    m_amr->allocate(m_velo_func, m_realm, m_phase, SpaceDim);
  }
  else{ 
    m_amr->allocate_ptr(m_velo_func);
  }

  // Only allocate memory if we actually have a diffusion solver
  if(m_diffusive){
    m_amr->allocate(m_diffco_cell, m_realm, m_phase, 1);
  }
  else{
    m_amr->allocate_ptr(m_diffco_cell);
  }

  m_particle_containers.emplace(which_container::bulk,    particle_container<ito_particle>());
  m_particle_containers.emplace(which_container::eb,      particle_container<ito_particle>());
  m_particle_containers.emplace(which_container::domain,  particle_container<ito_particle>());
  m_particle_containers.emplace(which_container::source,  particle_container<ito_particle>());
  m_particle_containers.emplace(which_container::covered, particle_container<ito_particle>());
  m_particle_containers.emplace(which_container::scratch, particle_container<ito_particle>());

  for (auto& container : m_particle_containers){
    m_amr->allocate(container.second, m_pvr_buffer, m_realm);
  }
}

void ito_solver::write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ito_solver::write_checkpoint_level");
  if(m_verbosity > 5){
    pout() << m_name + "::write_checkpoint_level" << endl;
  }

  // Write state. 
  write(a_handle, *m_state[a_level], m_name);

  // Write particles.
  if(m_checkpointing == which_checkpoint::particles){
    this->write_checkpoint_level_particles(a_handle, a_level);
  }
  else{ // In this case we need to write the number of physical particles in a grid cell. 
    this->write_checkpoint_level_fluid(a_handle, a_level);
  }
}

void ito_solver::write_checkpoint_level_particles(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ito_solver::write_checkpoint_level_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::write_checkpoint_level_particles" << endl;
  }

  const int halo        = 0;
  const std::string str = m_name + "_particles";

  particle_container<simple_ito_particle> realmParticles;
  m_amr->allocate(realmParticles,  m_pvr_buffer, m_realm);

  const particle_container<ito_particle>& myParticles = this->get_particles(which_container::bulk);

  // Make ito_particle into simple_ito_particle. This saves a shitload of disk space. 
  for (DataIterator dit(m_amr->get_grids(m_realm)[a_level]); dit.ok(); ++dit){
    List<simple_ito_particle>& other_particles = (realmParticles[a_level])[dit()].listItems();

    other_particles.clear();
      
    for (ListIterator<ito_particle> lit(myParticles[a_level][dit()].listItems()); lit.ok(); ++lit){
      other_particles.append(simple_ito_particle(lit().mass(), lit().position(), lit().energy()));
    }
  }

  writeParticlesToHDF(a_handle, realmParticles[a_level], str);
}

void ito_solver::write_checkpoint_level_fluid(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ito_solver::write_checkpoint_level_fluid");
  if(m_verbosity > 5){
    pout() << m_name + "::write_checkpoint_level_fluid" << endl;
  }

  const std::string str = m_name + "_particles";

  const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[a_level];
  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_level];
    
  const int comp     = 0;
  const int ncomp    = 1;
  const RealVect dx  = m_amr->get_dx()[a_level]*RealVect::Unit;
  const RealVect plo = m_amr->get_prob_lo();
  const IntVect ghos = IntVect::Zero;

  // Relevant container
  const particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  // Make something that can hold the particle numbers (stored as a Real)
  EBCellFactory fact(ebisl);
  LevelData<EBCellFAB> particleNumbers(dbl, ncomp, ghos, fact);
  EBLevelDataOps::setVal(particleNumbers, 0.0);

  // Now go through the grid and add the number of particles in each cell
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    BaseFab<Real>& pNum = particleNumbers[dit()].getSingleValuedFAB(); // No multivalued cells please. 

    // Get cell particles
    BinFab<ito_particle> pCel(dbl.get(dit()), dx, plo);
    pCel.addItems(particles[a_level][dit()].listItems());
				 
    for (BoxIterator bit(dbl.get(dit())); bit.ok(); ++bit){
      pNum(bit(), comp) = 0.0;
      for (ListIterator<ito_particle> lit(pCel(bit(), comp)); lit.ok(); ++lit){
	pNum(bit()) += lit().mass();
      }
    }
  }

  write(a_handle, particleNumbers, str);
}

void ito_solver::read_checkpoint_level(HDF5Handle& a_handle, const int a_level){
  CH_TIME("ito_solver::read_checkpoint_level");
  if(m_verbosity > 5){
    pout() << m_name + "::read_checkpoint_level" << endl;
  }

  // Read state vector
  read<EBCellFAB>(a_handle, *m_state[a_level], m_name, m_amr->get_grids(m_realm)[a_level], Interval(0,0), false);

  // Read particles.
  const std::string str = m_name + "_particles";
  if(m_checkpointing == which_checkpoint::particles){

    // Get simple_ito_particles from data file
    Vector<RefCountedPtr<ParticleData<simple_ito_particle> > > simpleParticles;
    m_amr->allocate(simpleParticles, m_realm);
    readParticlesFromHDF(a_handle, *simpleParticles[a_level], str);

    particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

    // Make simple_ito_particles into ito_ particles
    for (DataIterator dit(m_amr->get_grids(m_realm)[a_level]); dit.ok(); ++dit){
      List<ito_particle>& particlesDit = particles[a_level][dit()].listItems();
      
      for (ListIterator<simple_ito_particle> lit((*simpleParticles[a_level])[dit()].listItems()); lit.ok(); ++lit){
	particlesDit.append(ito_particle(lit().mass(), lit().position(), RealVect::Zero, 0.0, 0.0, lit().energy()));
      }
    }
  }
  else{

    // Read particles per cell
    EBAMRCellData ppc;
    m_amr->allocate(ppc, m_realm, m_phase, 1, 0);
    read<EBCellFAB>(a_handle, *ppc[a_level], str, m_amr->get_grids(m_realm)[a_level], Interval(0,0), false);
    
    // Restart particles
    this->restart_particles(*ppc[a_level], a_level);
  }
    
}

void ito_solver::restart_particles(LevelData<EBCellFAB>& a_num_particles, const int a_level){
  CH_TIME("ito_solver::restart_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::restart_particles" << endl;
  }

  const Real dx          = m_amr->get_dx()[a_level];
  const RealVect prob_lo = m_amr->get_prob_lo();

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_level];

  particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box           = dbl.get(dit());
    const EBISBox& ebisbox   = m_amr->get_ebisl(m_realm, m_phase)[a_level][dit()];
    const BaseFab<Real>& ppc = a_num_particles[dit()].getSingleValuedFAB();

    // Clear just to be safe
    List<ito_particle>& myParticles = particles[a_level][dit()].listItems();
    myParticles.clear();

    // Do regular cells
    for (BoxIterator bit(box); bit.ok(); ++bit){
      if(ebisbox.isRegular(bit())){

	// Compute weights and remainder
	const unsigned long long numParticles = (unsigned long long) llround(ppc(bit()));
	unsigned long long w, N, r;
	data_ops::compute_particle_weights(w, N, r, numParticles, m_ppc_restart);

	const RealVect minLo = -0.5*RealVect::Unit;
	const RealVect minHi =  0.5*RealVect::Unit;
	const RealVect norma = RealVect::Zero;
	const RealVect centr = RealVect::Zero;
	const RealVect pos   = prob_lo + (RealVect(bit()) + 0.5*RealVect::Unit)*dx;
	const Real kappa     = 1.0;


	// Now add N partices. If r > 0 we add another one with weight w + r
	for (unsigned long long i = 0; i < N; i++){
	  const RealVect particlePosition = this->random_position(pos, minLo, minHi, centr, norma, dx, kappa);
	  const Real particleWeight       = (Real) w;
	  
	  myParticles.add(ito_particle(particleWeight, particlePosition));
	}

	if(r > 0){
	  const RealVect particlePosition = this->random_position(pos, minLo, minHi, centr, norma, dx, kappa);
	  const Real particleWeight       = (Real) (w+r);

	  myParticles.add(ito_particle(particleWeight, particlePosition));
	}
      }
    }

    // Do irregular cells
    VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_level])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect iv    = vof.gridIndex();
      const RealVect cent = ebisbox.bndryCentroid(vof);
      const RealVect norm = ebisbox.normal(vof);
      const RealVect pos  = prob_lo + dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const Real kappa    = ebisbox.volFrac(vof);


      // Compute a small box that encloses the cut-cell volume
      RealVect minLo = -0.5*RealVect::Unit;
      RealVect minHi =  0.5*RealVect::Unit;
      if(kappa < 1.0){
	data_ops::compute_min_valid_box(minLo, minHi, norm, cent);
      }

      // Compute weights and remainder
      const unsigned long long numParticles = (unsigned long long) llround(ppc(iv));
      unsigned long long w, N, r;
      data_ops::compute_particle_weights(w, N, r, numParticles, m_ppc_restart);

      // Now add N partices. If r > 0 we add another one with weight w + r
      for (unsigned long long i = 0; i < N; i++){
	const RealVect particlePosition = this->random_position(pos, minLo, minHi, cent, norm, dx, kappa);
	const Real particleWeight       = (Real) w;
	  
	myParticles.add(ito_particle(particleWeight, particlePosition));
      }

      if(r > 0){
	const RealVect particlePosition = this->random_position(pos, minLo, minHi, cent, norm, dx, kappa);
	const Real particleWeight       = (Real) (w+r);

	myParticles.add(ito_particle(particleWeight, particlePosition));
      }
    }
  }
}

void ito_solver::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("ito_solver::write_plot_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_data" << endl;
  }

  // Write phi
  if(m_plot_phi){
    const Interval src(0, 0);
    const Interval dst(a_comp, a_comp);

    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      if(m_realm == a_output.get_realm()){
	m_state[lvl]->localCopyTo(src, *a_output[lvl], dst);
      }
      else{
	m_state[lvl]->copyTo(src, *a_output[lvl], dst);
      }
    }
    //data_ops::set_covered_value(a_output, a_comp, 0.0);
    a_comp++;
  }

  // Plot diffusion coefficient
  if(m_plot_dco && m_diffusive){
    const Interval src(0, 0);
    const Interval dst(a_comp, a_comp);
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      if(m_realm == a_output.get_realm()){
	m_diffco_cell[lvl]->localCopyTo(src, *a_output[lvl], dst);
      }
      else{
	m_diffco_cell[lvl]->copyTo(src, *a_output[lvl], dst);
      }
    }
    data_ops::set_covered_value(a_output, a_comp, 0.0);
    a_comp++;
  }

  // Write velocities
  if(m_plot_vel && m_mobile){
    const int ncomp = SpaceDim;
    const Interval src(0, ncomp-1);
    const Interval dst(a_comp, a_comp + ncomp-1);

    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      if(m_realm == a_output.get_realm()){
	m_velo_func[lvl]->localCopyTo(src, *a_output[lvl], dst);
      }
      else{
	m_velo_func[lvl]->copyTo(src, *a_output[lvl], dst);
      }
    }

    for (int c = 0; c < SpaceDim; c++){
      data_ops::set_covered_value(a_output, a_comp + c, 0.0);
    }

    a_comp += ncomp;
  }

  if(m_plot_particles){
    this->deposit_particles(m_scratch, m_particle_containers.at(which_container::bulk), m_plot_deposition);
    this->write_data(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_eb_particles){
    this->deposit_particles(m_scratch, m_particle_containers.at(which_container::eb), m_plot_deposition);
    this->write_data(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_domain_particles){
    this->deposit_particles(m_scratch, m_particle_containers.at(which_container::domain), m_plot_deposition);
    this->write_data(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_source_particles){
    this->deposit_particles(m_scratch, m_particle_containers.at(which_container::source), m_plot_deposition);
    this->write_data(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_energy_density){
    this->deposit_energy_density(m_scratch, m_particle_containers.at(which_container::bulk), m_plot_deposition);
    this->write_data(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_average_energy){
    this->sort_particles_by_cell(which_container::bulk);
    this->compute_average_energy(m_scratch, m_particle_containers.at(which_container::bulk));
    this->write_data(a_output, a_comp, m_scratch,  false);
    this->sort_particles_by_patch(which_container::bulk);
  }
}

void ito_solver::write_data(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp){
  CH_TIME("ito_solver::write_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_data" << endl;
  }

  const int comp = 0;
  const int ncomp = a_data[0]->nComp();

  const Interval src_interv(0, ncomp-1);
  const Interval dst_interv(a_comp, a_comp + ncomp - 1);

  // Copy data onto scratch
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, ncomp);
  data_ops::copy(scratch, a_data);

  // Interp if we should
  if(a_interp){
    m_amr->interpolate_to_centroids(scratch, m_realm, phase::gas);
  }

  m_amr->average_down(scratch, m_realm, m_phase);
  m_amr->interp_ghost(scratch, m_realm, m_phase);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    if(m_realm == a_output.get_realm()){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else{
      scratch[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }

  data_ops::set_covered_value(a_output, a_comp, 0.0);

  a_comp += ncomp;
}

void ito_solver::set_mass_to_conductivity(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::unset_mass_to_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::unset_mass_to_conductivity" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	ito_particle& p = lit();
	p.tmp()   = p.mass();
	p.mass() *= p.mobility();
      }
    }
  }
}

void ito_solver::unset_mass_to_conductivity(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::unset_mass_to_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::unset_mass_to_conductivity" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	ito_particle& p = lit();
	p.mass() = p.tmp();
      }
    }
  }
}

void ito_solver::set_mass_to_diffusivity(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::unset_mass_to_diffusivity");
  if(m_verbosity > 5){
    pout() << m_name + "::unset_mass_to_diffusivity" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	ito_particle& p = lit();
	p.tmp()   = p.mass();
	p.mass() *= p.diffusion();
      }
    }
  }
}

void ito_solver::unset_mass_to_diffusivity(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::unset_mass_to_diffusivity");
  if(m_verbosity > 5){
    pout() << m_name + "::unset_mass_to_diffusivity" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	ito_particle& p = lit();
	p.mass() = p.tmp();
      }
    }
  }
}

void ito_solver::set_mass_to_energy(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::set_mass_to_energy");
  if(m_verbosity > 5){
    pout() << m_name + "::set_mass_to_energy" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	ito_particle& p = lit();
	p.tmp()   = p.mass();
	p.mass() *= p.energy();
      }
    }
  }
}

void ito_solver::unset_mass_to_energy(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::unset_mass_to_energy");
  if(m_verbosity > 5){
    pout() << m_name + "::unset_mass_to_energy" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ito_particle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	ito_particle& p = lit();
	p.mass() = p.tmp();
      }
    }
  }
}

void ito_solver::deposit_conductivity(){
  CH_TIME("ito_solver::deposit_conductivity()");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_conductivity()" << endl;
  }

  this->deposit_conductivity(m_state, m_particle_containers.at(which_container::bulk));
}

void ito_solver::deposit_conductivity(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::deposit_conductivity(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_conductivity(state, particles)" << endl;
  }

  this->deposit_conductivity(a_state, a_particles, m_deposition);
}

void ito_solver::deposit_conductivity(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles, const DepositionType::Which a_deposition){
  CH_TIME("ito_solver::deposit_conductivity(state, particles, deposition_type)");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_conductivity(state, particles, deposition_type)" << endl;
  }

  this->set_mass_to_conductivity(a_particles);                 // Make mass = mass*mu
  this->deposit_particles(a_state, a_particles, a_deposition); // Deposit mass*mu
  this->unset_mass_to_conductivity(a_particles);               // Make mass = mass/mu
}

void ito_solver::deposit_diffusivity(){
  CH_TIME("ito_solver::deposit_diffusivity()");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_diffusivity()" << endl;
  }

  this->deposit_diffusivity(m_state, m_particle_containers.at(which_container::bulk));
}

void ito_solver::deposit_diffusivity(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::deposit_diffusivity(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_diffusivity(state, particles)" << endl;
  }

  this->deposit_diffusivity(a_state, a_particles, m_deposition);
}

void ito_solver::deposit_diffusivity(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles, const DepositionType::Which a_deposition){
  CH_TIME("ito_solver::deposit_diffusivity(state, particles, deposition_type)");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_diffusivity(state, particles, deposition_type)" << endl;
  }

  this->set_mass_to_diffusivity(a_particles);                                  // Make mass = mass*D
  this->deposit_particles(a_state, a_particles, a_deposition); // Deposit mass*D
  this->unset_mass_to_diffusivity(a_particles);                                // Make mass = mass/D
}

void ito_solver::deposit_energy_density(){
  CH_TIME("ito_solver::deposit_energy_density()");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_energy_density()" << endl;
  }

  this->deposit_energy_density(m_state, m_particle_containers.at(which_container::bulk));
}

void ito_solver::deposit_energy_density(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::deposit_energy_density(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_energy_density(state, particles)" << endl;
  }

  this->deposit_energy_density(a_state, a_particles, m_deposition);
}

void ito_solver::deposit_energy_density(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles, const DepositionType::Which a_deposition){
  CH_TIME("ito_solver::deposit_energy_density(state, particles, deposition_type)");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_energy_density(state, particles, deposition_type)" << endl;
  }

  this->set_mass_to_energy(a_particles);                       // Make mass = mass*E
  this->deposit_particles(a_state, a_particles, a_deposition); // Deposit mass*E
  this->unset_mass_to_energy(a_particles);                     // Make mass = mass/E
}


void ito_solver::compute_average_mobility(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::compute_average_mobility(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_average_mobility(state, particles)" << endl;
  }

  constexpr int comp = 0;

  data_ops::set_value(a_state, 0.0);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const BinFab<ito_particle>& particles = a_particles.get_cell_particles(lvl, dit());

      BaseFab<Real>& state = (*a_state[lvl])[dit()].getSingleValuedFAB();

      for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit){
	const IntVect iv = bit();

	const List<ito_particle>& listParticles = particles(iv, comp);

	if(listParticles.length() > 0){
	  Real M  = 0.0;
	  Real MU = 0.0;
	  for (ListIterator<ito_particle> lit(listParticles); lit.ok(); ++lit){
	    const ito_particle& p = lit();
	    const Real m  = p.mass();
	    const Real mu = p.mobility();

	    M  += m;
	    MU += m*mu;
	  }

	  state(iv, comp) = MU/M;
	}
      }
    }
  }
}

void ito_solver::compute_average_diffusion(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::compute_average_diffusion(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_average_diffusion(state, particles)" << endl;
  }

  constexpr int comp = 0;

  data_ops::set_value(a_state, 0.0);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const BinFab<ito_particle>& particles = a_particles.get_cell_particles(lvl, dit());

      BaseFab<Real>& state = (*a_state[lvl])[dit()].getSingleValuedFAB();

      for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit){
	const IntVect iv = bit();

	const List<ito_particle>& listParticles = particles(iv, comp);

	if(listParticles.length() > 0){
	  Real M = 0.0;
	  Real D = 0.0;
	  for (ListIterator<ito_particle> lit(listParticles); lit.ok(); ++lit){
	    const ito_particle& p = lit();
	    const Real m = p.mass();
	    const Real d = p.diffusion();

	    M += m;
	    D += m*d;
	  }

	  state(iv, comp) = D/M;
	}
      }
    }
  }
}

void ito_solver::compute_average_energy(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::compute_average_energy(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_average_energy(state, particles)" << endl;
  }

  constexpr int comp = 0;

  data_ops::set_value(a_state, 0.0);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const BinFab<ito_particle>& particles = a_particles.get_cell_particles(lvl, dit());

      BaseFab<Real>& state = (*a_state[lvl])[dit()].getSingleValuedFAB();

      for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit){
	const IntVect iv = bit();

	const List<ito_particle>& listParticles = particles(iv, comp);

	if(listParticles.length() > 0){
	  Real M = 0.0;
	  Real E = 0.0;
	  for (ListIterator<ito_particle> lit(listParticles); lit.ok(); ++lit){
	    const ito_particle& p = lit();
	    const Real m = p.mass();
	    const Real e = p.energy();

	    M += m;
	    E += m*e;
	  }

	  state(iv, comp) = E/M;
	}
      }
    }
  }
}

void ito_solver::deposit_particles(){
  CH_TIME("ito_solver::deposit_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_particles" << endl;
  }

  this->deposit_particles(which_container::bulk);
}

void ito_solver::deposit_particles(const which_container a_container){
  CH_TIME("ito_solver::deposit_particles(container)");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_particles(container)" << endl;
  }

  this->deposit_particles(m_state, m_particle_containers.at(a_container), m_deposition);
}

void ito_solver::deposit_nonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC){
  CH_TIME("ito_solver::deposit_nonConservative");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_nonConservative" << endl;
  }

  const std::string cur_realm = a_depositionNC.get_realm();

  if(m_blend_conservation){
    irreg_amr_stencil<noncons_div>& stencils = m_amr->get_noncons_div_stencils(cur_realm, m_phase);
    stencils.apply(a_depositionNC, a_depositionKappaC);
  }
  else{
    data_ops::set_value(a_depositionNC, 0.0);
  }
}

void ito_solver::deposit_hybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_mass_diff, const EBAMRIVData& a_depositionNC){
  CH_TIME("ito_solver::deposit_hybrid");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_hybrid" << endl;
  }

  const std::string cur_realm = a_depositionH.get_realm();

  const int comp  = 0;
  const int ncomp = 1;


  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(cur_realm)[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(cur_realm, m_phase)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& divH               = (*a_depositionH[lvl])[dit()];  // On input, this contains kappa*depositionWeights
      BaseIVFAB<Real>& deltaM       = (*a_mass_diff[lvl])[dit()];
      const BaseIVFAB<Real>& divNC  = (*a_depositionNC[lvl])[dit()]; 
      const EBISBox& ebisbox        = ebisl[dit()];

      VoFIterator& vofit = (*m_amr->get_vofit(cur_realm, m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa    = ebisbox.volFrac(vof);
	const Real dc       = divH(vof, comp);
	const Real dnc      = divNC(vof, comp);

	// Note that if dc - kappa*dnc can be negative, i.e. we may end up STEALING mass
	// from other cells. This is why there is a flag m_blend_conservation which always
	// gives positive definite results. 
	divH(vof, comp)   = dc + (1.0-kappa)*dnc;        // On output, contains hybrid divergence
	deltaM(vof, comp) = (1-kappa)*(dc - kappa*dnc);  // Remember, dc already scaled by kappa.
      }
    }
  }
}


void ito_solver::increment_redist(const EBAMRIVData& a_mass_diff){
  CH_TIME("ito_solver::increment_redist");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist" << endl;
  }

  const std::string cur_realm = a_mass_diff.get_realm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(cur_realm)[lvl];
    
    EBLevelRedist& level_redist = *(m_amr->get_level_redist(cur_realm, m_phase)[lvl]);
    level_redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      level_redist.increment((*a_mass_diff[lvl])[dit()], dit(), interv);
    }
  }
}

void ito_solver::level_redistribution(EBAMRCellData& a_state){
  CH_TIME("ito_solver::level_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::level_redistribution" << endl;
  }

  const std::string cur_realm = a_state.get_realm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& level_redist = *(m_amr->get_level_redist(cur_realm, m_phase)[lvl]);
    level_redist.redistribute(*a_state[lvl], interv);
    level_redist.setToZero();
  }
}

void ito_solver::coarse_fine_increment(const EBAMRIVData& a_mass_diff){
  CH_TIME("ito_solver::coarse_fine_increment");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_increment" << endl;
  }

  const std::string cur_realm = a_mass_diff.get_realm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(0,0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(cur_realm)[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(cur_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(cur_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(cur_realm, m_phase)[lvl];

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < 0;

    if(has_coar){
      fine2coar_redist->setToZero();

    }
    if(has_fine){
      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      if(has_coar){
	fine2coar_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
      }

      if(has_fine){
	coar2fine_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
	coar2coar_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
      }
    }
  }
}


void ito_solver::coarse_fine_redistribution(EBAMRCellData& a_state){
  CH_TIME("ito_solver::coarse_fine_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_redistribution" << endl;
  }

  const std::string cur_realm = a_state.get_realm();

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);


  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->get_dx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(cur_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(cur_realm, m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(cur_realm, m_phase)[lvl];
    if(has_coar){
      fine2coar_redist->redistribute(*a_state[lvl-1], interv);
      fine2coar_redist->setToZero();
    }

    if(has_fine){
      coar2fine_redist->redistribute(*a_state[lvl+1], interv);
      coar2coar_redist->redistribute(*a_state[lvl],   interv);

      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }
  }
}

void ito_solver::deposit_weights(EBAMRCellData& a_state, const particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::deposit_weights");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_weights" << endl;
  }

  this->deposit_particles(a_state, m_particle_containers.at(which_container::bulk), DepositionType::NGP);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    data_ops::scale(*a_state[lvl], pow(m_amr->get_dx()[lvl], SpaceDim));
  }
}

void ito_solver::add_particles(ListBox<ito_particle>& a_part, const int a_lvl, const DataIndex a_dit, const bool a_destructive){
  CH_TIME("ito_solver::add_particles(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::add_particles(lvl, dit)" << endl;
  }

  particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  ListBox<ito_particle>& my_particles = particles[a_lvl][a_dit];

  if(a_destructive){
    my_particles.addItemsDestructive(a_part.listItems());
  }
  else{
    my_particles.addItems(a_part.listItems());
  }
}

bool ito_solver::is_mobile() const{
  CH_TIME("ito_solver::is_mobile");
  
  return m_mobile;
}
  

bool ito_solver::is_diffusive() const{
  CH_TIME("ito_solver::is_diffusive");
  
  return m_diffusive;
}

void ito_solver::pre_regrid(const int a_base, const int a_old_finest_level){
  CH_TIME("ito_solver::pre_regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::pre_regrid" << endl;
  }

  for (auto& container : m_particle_containers){
    particle_container<ito_particle>& particles = container.second;

    particles.pre_regrid(a_base);
  }
}

particle_container<ito_particle>& ito_solver::get_particles(const which_container a_container){
  return m_particle_containers.at(a_container);
}

const particle_container<ito_particle>& ito_solver::get_particles(const which_container a_container) const {
  return m_particle_containers.at(a_container);
}

EBAMRCellData& ito_solver::get_state(){
  CH_TIME("ito_solver::get_state");
  if(m_verbosity > 5){
    pout() << m_name + "::get_state" << endl;
  }

  return m_state;
}

EBAMRCellData& ito_solver::get_velo_func(){
  CH_TIME("ito_solver::get_velo_func");
  if(m_verbosity > 5){
    pout() << m_name + "::get_velo_func" << endl;
  }

  return m_velo_func;
}

EBAMRCellData& ito_solver::get_diffco_func(){
  CH_TIME("ito_solver::get_diffco_func");
  if(m_verbosity > 5){
    pout() << m_name + "::get_diffco_func" << endl;
  }

  return m_diffco_cell;
}

EBAMRCellData& ito_solver::get_scratch(){
  CH_TIME("ito_solver::get_scratch");
  if(m_verbosity > 5){
    pout() << m_name + "::get_scratch" << endl;
  }

  return m_scratch;
}

EBAMRCellData& ito_solver::get_mobility_func(){
  CH_TIME("ito_solver::get_mobility_func");
  if(m_verbosity > 5){
    pout() << m_name + "::get_mobility_func" << endl;
  }

  return m_mobility_func;
}

void ito_solver::set_diffco_func(const Real a_diffco){
  CH_TIME("ito_solver::set_diffco_func");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco_func" << endl;
  }

  data_ops::set_value(m_diffco_cell, a_diffco);
}

void ito_solver::set_velocity_func(const RealVect a_vel){
  CH_TIME("ito_solver::set_velocity_func");
  if(m_verbosity > 5){
    pout() << m_name + "::set_velocity_func" << endl;
  }

  for (int comp = 0; comp < SpaceDim; comp++){
    data_ops::set_value(m_velo_func, a_vel[comp], comp);
  }
}

void ito_solver::set_mobility(const Real a_mobility){
  CH_TIME("ito_solver::set_mobility");
  if(m_verbosity > 5){
    pout() << m_name + "::set_mobility" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(which_container::bulk);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ito_particle>& particlesDit = particles[lvl][dit()].listItems();

      for (ListIterator<ito_particle> lit(particlesDit); lit.ok(); ++lit){
	lit().mobility() = a_mobility;
      }
    }
  }
}

void ito_solver::interpolate_velocities(){
  CH_TIME("ito_solver::interpolate_velocities");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_velocities" << endl;
  }

  if(m_mobile){
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	this->interpolate_velocities(lvl, dit());
      }
    }
  }
}

void ito_solver::interpolate_velocities(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ito_solver::interpolate_velocities");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_velocities" << endl;
  }

  particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  if(m_mobile){
    const EBCellFAB& velo_func = (*m_velo_func[a_lvl])[a_dit];
    const EBISBox& ebisbox     = velo_func.getEBISBox();
    const FArrayBox& vel_fab   = velo_func.getFArrayBox();
    const RealVect dx          = m_amr->get_dx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->get_prob_lo();
    const Box box              = m_amr->get_grids(m_realm)[a_lvl][a_dit];

    List<ito_particle>& particleList = particles[a_lvl][a_dit].listItems();

    // This interpolates the velocity function on to the particle velocities
    EBParticleInterp meshInterp(box, ebisbox, dx, origin, m_irreg_ngp_interpolation);
    meshInterp.interpolateVelocity(particleList, vel_fab, m_deposition);

    // Go through the particles and set their velocities to velo_func*mobility
    for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
      ito_particle& p = lit();
      p.velocity() *= p.mobility();
    }
  }
}

void ito_solver::interpolate_mobilities(){
  CH_TIME("ito_solver::interpolate_mobilities()");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_mobilities()" << endl;
  }

  if(m_mobile){

    if(m_mobility_interp == mobility_interp::velocity){
      data_ops::vector_length(m_scratch, m_velo_func); // Compute |E| (or whatever other function you've decided to provide).
      m_amr->average_down(m_scratch, m_realm, m_phase);
      m_amr->interp_ghost(m_scratch, m_realm, m_phase);
    }
      
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
      
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	this->interpolate_mobilities(lvl, dit());
      }
    }
  }
}

void ito_solver::interpolate_mobilities(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ito_solver::interpolate_mobilities(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_mobilities(lvl, dit)" << endl;
  }


  if(m_mobility_interp == mobility_interp::mobility){
    this->interpolate_mobilities_mu(a_lvl, a_dit);
  }
  else if (m_mobility_interp == mobility_interp::velocity){
    this->interpolate_mobilities_vel(a_lvl, a_dit);
  }

}

void ito_solver::interpolate_mobilities_mu(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ito_solver::interpolate_mobilities_mu(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_mobilities_mu(lvl, dit)" << endl;
  }

  particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  if(m_mobile){
    const EBCellFAB& mob_func  = (*m_mobility_func[a_lvl])[a_dit];
    const EBISBox& ebisbox     = mob_func.getEBISBox();
    const FArrayBox& mob_fab   = mob_func.getFArrayBox();
    const RealVect dx          = m_amr->get_dx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->get_prob_lo();
    const Box box              = m_amr->get_grids(m_realm)[a_lvl][a_dit];

    List<ito_particle>& particleList = particles[a_lvl][a_dit].listItems();
    EBParticleInterp meshInterp(box, ebisbox, dx, origin, m_irreg_ngp_interpolation);
    
    meshInterp.interpolateMobility(particleList, mob_fab, m_deposition);
  }
}

void ito_solver::interpolate_mobilities_vel(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ito_solver::interpolate_mobilities_vel(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_mobilities_vel(lvl, dit)" << endl;
  }

  particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  if(m_mobile){
    const EBCellFAB& mob_func  = (*m_mobility_func[a_lvl])[a_dit];
    const EBISBox& ebisbox     = mob_func.getEBISBox();
    const FArrayBox& mob_fab   = mob_func.getFArrayBox();
    const RealVect dx          = m_amr->get_dx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->get_prob_lo();
    const Box box              = m_amr->get_grids(m_realm)[a_lvl][a_dit];

    FArrayBox& scratch = (*m_scratch[a_lvl])[a_dit].getFArrayBox();

    List<ito_particle>& particleList = particles[a_lvl][a_dit].listItems();
    EBParticleInterp meshInterp(box, ebisbox, dx, origin, m_irreg_ngp_interpolation);
    
    // First, interpolate |E| to the particle position, it will be stored on m_tmp. 
    meshInterp.interpolateMobility(particleList, scratch, m_deposition);

    for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
      ito_particle& p = lit();
      
      p.tmp() = p.mobility();
    }

    // This interpolates mu*|E| to the particle position and stores it on the mobility. After that, we compute mu_p = (mu*E)/E
    scratch *= mob_fab;
    meshInterp.interpolateMobility(particleList, scratch, m_deposition);
    for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
      ito_particle& p = lit();

      p.mobility() *= 1./p.tmp();
    }
  }
}

void ito_solver::update_mobilities(){
  CH_TIME("ito_solver::update_mobilities");
  if(m_verbosity > 5){
    pout() << m_name + "::update_mobilities" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      this->update_mobilities(lvl, dit());
    }
  }
}

void ito_solver::update_mobilities(const int a_level, const DataIndex a_dit){
  CH_TIME("ito_solver::update_mobilities(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::update_mobilities(lvl, dit)" << endl;
  }

  particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  if(m_mobile){
    List<ito_particle>& particleList = particles[a_level][a_dit].listItems();

    for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
      ito_particle& p = lit();
      
      p.mobility() = m_species->mobility(p.energy());
    }
  }
}

void ito_solver::interpolate_diffusion(){
  CH_TIME("ito_solver::interpolate_diffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_diffusion" << endl;
  }

  if(m_diffusive){
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	this->interpolate_diffusion(lvl, dit());
      }
    }
  }
}

void ito_solver::interpolate_diffusion(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ito_solver::interpolate_diffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_diffusion" << endl;
  }

  particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  if(m_diffusive){
    const EBCellFAB& dco_cell  = (*m_diffco_cell[a_lvl])[a_dit];
    const EBISBox& ebisbox     = dco_cell.getEBISBox();
    const FArrayBox& dco_fab   = dco_cell.getFArrayBox();
    const RealVect dx          = m_amr->get_dx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->get_prob_lo();
    const Box box              = m_amr->get_grids(m_realm)[a_lvl][a_dit];

    List<ito_particle>& particleList = particles[a_lvl][a_dit].listItems();

    EBParticleInterp meshInterp(box, ebisbox,dx, origin, m_irreg_ngp_interpolation);
    meshInterp.interpolateDiffusion(particleList, dco_fab, m_deposition);
  }
}

void ito_solver::update_diffusion(){
  CH_TIME("ito_solver::update_diffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::update_diffusion" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      this->update_diffusion(lvl, dit());
    }
  }
}

void ito_solver::update_diffusion(const int a_level, const DataIndex a_dit){
  CH_TIME("ito_solver::update_diffusion(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::update_diffusion(lvl, dit)" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(which_container::bulk);

  if(m_diffusive){
    List<ito_particle>& particleList = particles[a_level][a_dit].listItems();

    for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
      ito_particle& p = lit();
      
      p.mobility() = m_species->diffusion(p.energy());
    }
  }
}



Real ito_solver::compute_dt() const {
  CH_TIME("ito_solver::compute_dt(allAMRlevels)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_dt(allAMRlevels)" << endl;
  }

  Real dt = 1.E99;
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real levelDt = this->compute_dt(lvl);
    dt = Min(dt, levelDt);
  }

  return dt;
}

Real ito_solver::compute_dt(const int a_lvl) const{
  CH_TIME("ito_solver::compute_dt(lvl)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_dt(lvl)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_lvl];
  const Real dx = m_amr->get_dx()[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->compute_dt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
    if(result != MPI_SUCCESS){
      MayDay::Error("ito_solver::compute_dt(lvl) - communication error on norm");
    }
    dt = tmp;
#endif  

  return dt;
}
  

Real ito_solver::compute_dt(const int a_lvl, const DataIndex a_dit, const Real a_dx) const{
  CH_TIME("ito_solver::compute_dt(lvl, dit, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_dt(lvl, dit, dx)" << endl;
  }

  Real dt = 1.E99;

  const particle_container<ito_particle>& particles = this->get_particles(which_container::bulk);

  const List<ito_particle>& particleList = particles[a_lvl][a_dit].listItems();
  
  ListIterator<ito_particle> lit(particleList);
  if(m_mobile && !m_diffusive){
    for (lit.rewind(); lit.ok(); ++lit){
      const ito_particle& p = particleList[lit];
      const RealVect& v = p.velocity();

      const int maxDir  = v.maxDir(true);
      const Real vMax   = Abs(v[maxDir]);
      const Real thisDt = (vMax > 0.0) ? a_dx/vMax : 1.E99;

      dt = Min(dt, thisDt);
    }
  }
  else if(!m_mobile && m_diffusive){
    for (lit.rewind(); lit.ok(); ++lit){
      const ito_particle& p = particleList[lit];

      const Real D      = p.diffusion();
      const Real thisDt = (D > 0.0) ? a_dx*a_dx/(2.0*D) : 1.E99;

      dt = Min(dt, thisDt);
    }
  }
  else if(m_mobile && m_diffusive){
    for (lit.rewind(); lit.ok(); ++lit){
      const ito_particle& p = particleList[lit];
      
      const RealVect& v = p.velocity();
      const int maxDir = v.maxDir(true);
      const Real vMax  = Abs(v[maxDir]);
      const Real D     = p.diffusion();

      const Real dtAdvect = (vMax > 0.0) ? a_dx/vMax         : 1.E99;
      const Real dtDiffus = (D > 0.0)    ? a_dx*a_dx/(2.0*D) : 1.E99;

      const Real thisDt = 1./(1./dtAdvect + 1./dtDiffus);
      
      dt = Min(dt, thisDt);
    }
  }

  return dt;
}

Real ito_solver::compute_min_dt(const Real a_maxCellsToMove) const {
  CH_TIME("ito_solver::compute_min_dt(allAMRlevels, maxCellsToMove)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_min_dt(allAMRlevels, maxCellsToMove)" << endl;
  }

  Real dt = 1.E99;
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real levelDt = this->compute_min_dt(a_maxCellsToMove, lvl);
    dt = Min(dt, levelDt);
  }

  return dt;
}

Real ito_solver::compute_min_dt(const Real a_maxCellsToMove, const int a_lvl) const{
  CH_TIME("ito_solver::compute_min_dt(maxCellsToMove, lvl)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_min_dt(maxCellsToMove, lvl)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_lvl];
  const Real dx = m_amr->get_dx()[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->compute_min_dt(a_maxCellsToMove, a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
    if(result != MPI_SUCCESS){
      MayDay::Error("ito_solver::compute_dt(lvl) - communication error on norm");
    }
    dt = tmp;
#endif  

  return dt;
}

Real ito_solver::compute_min_dt(const Real a_maxCellsToMove, const int a_lvl, const DataIndex a_dit, const Real a_dx) const{
  CH_TIME("ito_solver::compute_min_dt(maxCellsToMove, lvl, dit, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_min_dt(maxCellsToMove, lvl, dit, dx)" << endl;
  }

  Real dt = 1.E99;

  const particle_container<ito_particle>& particles = this->get_particles(which_container::bulk);

  const Real dMax  = a_maxCellsToMove*a_dx;
  const Real dMax2 = dMax*dMax;
  const Real W0    = m_normal_max;
  const Real W02   = m_normal_max*m_normal_max;

  const List<ito_particle>& particleList = particles[a_lvl][a_dit].listItems();
  
  ListIterator<ito_particle> lit(particleList);
  if(m_mobile && !m_diffusive){
    for (lit.rewind(); lit; ++lit){
      const ito_particle& p = particleList[lit];
      const RealVect& v = p.velocity();

      const int maxDir = v.maxDir(true);
      const Real thisDt = dMax/Abs(v[maxDir]);

      dt = Min(dt, thisDt);
    }
  }
  else if(!m_mobile && m_diffusive){
    for (lit.rewind(); lit; ++lit){
      const ito_particle& p = particleList[lit];
      
      const Real thisDt = dMax2/(2.0*p.diffusion()*W02);
      dt = Min(dt, thisDt);
    }
  }
  else if(m_mobile && m_diffusive){
    for (lit.rewind(); lit; ++lit){
      const ito_particle& p = particleList[lit];
      
      const RealVect& v = p.velocity();
      const int maxDir = v.maxDir(true);
      const Real vMax  = Abs(v[maxDir]);
      const Real D     = p.diffusion();

      
      if(vMax > 0.0){
	const Real a = vMax;
	const Real b = W0*sqrt(2.0*D);
	const Real c = dMax;
	
	const Real A = a*a;
	const Real B = -(b*b + 2*a*c);
	const Real C = c*c;

	const Real thisDt = (-B - sqrt(B*B - 4.*A*C))/(2.*A);
	dt = Min(dt, thisDt);
      }
      else{
	if(D > 0.0){
	  const Real thisDt = dMax2/(2.0*D*W02);
	  dt = Min(dt, thisDt);
	}
      }
    }
  }

  return dt;

}

Real ito_solver::compute_min_drift_dt(const Real a_maxCellsToMove) const{
  CH_TIME("ito_solver::compute_min_drift_dt(allAMRlevels, maxCellsToMove)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_min_drift_dt(allAMRlevels, maxCellsToMove)" << endl;
  }

  Vector<Real> dt = this->compute_drift_dt(a_maxCellsToMove);

  Real minDt = dt[0];
  for (int lvl = 0; lvl < dt.size(); lvl++){
    minDt = Min(minDt, dt[lvl]);
  }

  return minDt;
}

Vector<Real> ito_solver::compute_drift_dt(const Real a_maxCellsToMove) const {
  CH_TIME("ito_solver::compute_drift_dt(amr, maxCellsToMove)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(amr, maxCellsToMove)" << endl;
  }

  Vector<Real> dt = this->compute_drift_dt();
  for (int lvl = 0; lvl < dt.size(); lvl++){
    dt[lvl] = dt[lvl]*a_maxCellsToMove;
  }

  return dt;

}

Real ito_solver::compute_advective_dt() const {
  CH_TIME("ito_solver::compute_advective_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_advective_dt()" << endl;
  }
  
  Real minDt = 1.E99;
  const Vector<Real> levelDts = this->compute_drift_dt();

  for (int lvl = 0; lvl < levelDts.size(); lvl++){
    minDt = Min(levelDts[lvl], minDt);
  }

  return minDt;
}

Vector<Real> ito_solver::compute_drift_dt() const {
  CH_TIME("ito_solver::compute_drift_dt(amr)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(amr)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  Vector<Real> dt(1 + finest_level, 1.2345E67);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    dt[lvl] = this->compute_drift_dt(lvl);
  }

  return dt;
}

Real ito_solver::compute_drift_dt(const int a_lvl) const {
  CH_TIME("ito_solver::compute_drift_dt(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(level)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_lvl];
  const RealVect dx = m_amr->get_dx()[a_lvl]*RealVect::Unit;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->compute_drift_dt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
    if(result != MPI_SUCCESS){
      MayDay::Error("ito_solver::compute_drift_level(lvl) - communication error on norm");
    }
    dt = tmp;
#endif  

  return dt;
}

Real ito_solver::compute_drift_dt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const{
  CH_TIME("ito_solver::compute_drift_dt(level, dataindex, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(level, dataindex, dx)" << endl;
  }

  const particle_container<ito_particle>& particles = m_particle_containers.at(which_container::bulk);

  constexpr Real safety = 1.E-10;

  const List<ito_particle>& particleList = particles[a_lvl][a_dit].listItems();

  Real dt = 1.E99;

  if(m_mobile){
    for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
      const ito_particle& p = particleList[lit];
      const RealVect& v     = p.velocity();
      const int maxDir      = v.maxDir(true);
      const Real thisDt     = a_dx[maxDir]/(safety + Abs(v[maxDir]));

      dt = Min(dt, thisDt);
    }
  }

  return dt;
}

Real ito_solver::compute_min_diffusion_dt(const Real a_maxCellsToHop) const{
  CH_TIME("ito_solver::compute_min_diffusion_dt(min, maxCellsToHop)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_min_diffusion_dt(min, maxCellsToHop)" << endl;
  }

  Vector<Real> dt = this->compute_diffusion_dt(a_maxCellsToHop);
  Real minDt = dt[0];
  for (int lvl = 0; lvl < dt.size(); lvl++){
    minDt = Min(minDt, dt[lvl]);
  }

  return minDt;
}

Vector<Real> ito_solver::compute_diffusion_dt(const Real a_maxCellsToHop) const{
  CH_TIME("ito_solver::compute_min_diffusion_dt(maxCellsToHop)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_min_diffusion_dt(maxCellsToHop)" << endl;
  }

  const Real factor  = a_maxCellsToHop/m_normal_max;
  const Real factor2 = factor*factor;
  
  Vector<Real> dt = this->compute_diffusion_dt();
  for (int lvl = 0; lvl < dt.size(); lvl++){
    dt[lvl] = dt[lvl]*factor2;
  }

  return dt;
}

Real ito_solver::compute_diffusive_dt() const {
  CH_TIME("ito_solver::compute_diffusive_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusive_dt" << endl;
  }

  Real minDt = 1.E99;
  const Vector<Real> levelDts = this->compute_diffusion_dt();

  for (int lvl = 0; lvl < levelDts.size(); lvl++){
    minDt = Min(levelDts[lvl], minDt);
  }

  return minDt;
}

Vector<Real> ito_solver::compute_diffusion_dt() const{
  CH_TIME("ito_solver::compute_diffusion_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_dt" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
    
  Vector<Real> dt(1 + finest_level, 1.2345E6);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    dt[lvl] = this->compute_diffusion_dt(lvl);
  }

  return dt;
}

Real ito_solver::compute_diffusion_dt(const int a_lvl) const{
  CH_TIME("ito_solver::compute_diffusion_dt(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_dt(level)" << endl;
  }

  
  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_lvl];
  const RealVect dx = m_amr->get_dx()[a_lvl]*RealVect::Unit;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->compute_diffusion_dt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("ito_solver::compute_diffusion_level(lvl) - communication error on norm");
  }
  dt = tmp;
#endif

  return dt;
}

Real ito_solver::compute_diffusion_dt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const{
  CH_TIME("ito_solver::compute_diffusion_dt(level, dataindex, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_dt(level, dataindex, dx)" << endl;
  }

  const particle_container<ito_particle>& particles = this->get_particles(which_container::bulk);
  
  const List<ito_particle>& particleList = particles[a_lvl][a_dit].listItems();

  Real dt = 1.E99;

  if(m_diffusive){
    for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
      const ito_particle& p = particleList[lit];
    
      const Real thisDt = a_dx[0]*a_dx[0]/(2.0*p.diffusion());
    
      dt = Min(dt, thisDt);
    }
  }

  return dt;
}

void ito_solver::remap(){
  CH_TIME("ito_solver::remap");
  if(m_verbosity > 5){
    pout() << m_name + "::remap" << endl;
  }

  this->remap(which_container::bulk);

}

void ito_solver::remap(const which_container a_container){
  CH_TIME("ito_solver::remap(container)");
  if(m_verbosity > 5){
    pout() << m_name + "::remap(container)" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(a_container);
  
  particles.remap();
}

DepositionType::Which ito_solver::get_deposition() const {
  return m_deposition;
}

phase::which_phase ito_solver::get_phase() const{
  return m_phase;
}

void ito_solver::sort_particles_by_cell(const which_container a_container){
  CH_TIME("ito_solver::sort_particles_by_cell(container)");
  if(m_verbosity > 5){
    pout() << m_name + "::sort_particles_by_cell(container)" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(a_container);

  particles.sort_particles_by_cell();
}


void ito_solver::sort_particles_by_patch(const which_container a_container){
  CH_TIME("ito_solver::sort_particles_by_patch(container)");
  if(m_verbosity > 5){
    pout() << m_name + "::sort_particles_by_patch(container)" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(a_container);

  particles.sort_particles_by_patch();
}

void ito_solver::make_superparticles(const which_container a_container, const int a_particlesPerPatch){
  CH_TIME("ito_solver::make_superparticles(int)");
  if(m_verbosity > 5){
    pout() << m_name + "::make_superparticles(int)" << endl;
  }
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    this->make_superparticles(a_container, a_particlesPerPatch, lvl);
  }

}

void ito_solver::make_superparticles(const which_container a_container, const int a_particlesPerPatch, const int a_level){
  CH_TIME("ito_solver::make_superparticles(int, level)");
  if(m_verbosity > 5){
    pout() << m_name + "::make_superparticles(int, level)" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(which_container::bulk);

  const DisjointBoxLayout& dbl = particles.get_grids()[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    this->make_superparticles(a_container, a_particlesPerPatch, a_level, dit());
  }
}

void ito_solver::make_superparticles(const which_container a_container, const int a_particlesPerCell, const int a_level, const DataIndex a_dit){
  CH_TIME("ito_solver::make_superparticles(int, level, patch)");
  if(m_verbosity > 5){
    pout() << m_name + "::make_superparticles(int, level, patch)" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(a_container);
  
  const int comp = 0;
  //  const Box box  = m_amr->get_grids(m_realm)[a_level].get(a_dit);

  const Box box  = particles.get_grids()[a_level][a_dit];

  // This are the particles in the box we're currently looking at. 
  BinFab<ito_particle>& cellParticles = particles.get_cell_particles(a_level, a_dit);

  // Iterate over particles
  for (BoxIterator bit(box); bit.ok(); ++bit){
    const IntVect iv = bit();
  
    List<ito_particle>& particles = cellParticles(iv, comp);

    if(particles.length() > 0){
      this->bvh_merge(particles, a_particlesPerCell);
    }
  }
}

void ito_solver::bvh_merge(List<ito_particle>& a_particles, const int a_particlesPerCell){
  CH_TIME("ito_solver::bvh_merge");

#if ITO_DEBUG
  Real mass_before = 0.0;
  Real energy_before = 0.0;
  
  Real mass_after = 0.0;
  Real energy_after = 0.0;
  RealVect center_before = RealVect::Zero;
  RealVect center_after = RealVect::Zero;
#endif
  
  std::vector<point_mass> pointMasses(a_particles.length());
  int i = 0;
  Real mass = 0.0;
  for (ListIterator<ito_particle> lit(a_particles); lit.ok(); ++lit, i++){
    const ito_particle& p = lit();
    pointMasses[i].define(p.position(), p.mass(), p.energy()); // Note, p.energy() is average energy and not total energy. 
    mass += p.mass();
    //    i++;

#if ITO_DEBUG
    mass_before   += p.mass();
    center_before += p.mass()*p.position();
    energy_before += p.mass()*p.energy();
#endif
  }

#if ITO_DEBUG
  center_before *= 1./mass_before;
  if(mass_before   < 1.0)            MayDay::Abort("ito_solver::bvh_merge - bad initial mass!");
  if(mass_before   != mass_before)   MayDay::Abort("ito_solver::bvh_merge - initial mass is NaN");
  if(energy_before != energy_before) MayDay::Abort("ito_solver::bvh_merge - initial energy is NaN");
#endif
  
  // 2. Build the BVH tree and get the leaves of the tree
  const int dir = (m_kd_direction < 0) ? m_udist0d(m_rng) : m_kd_direction;
  m_tree.define(pointMasses, mass);
  m_tree.build_tree(dir, a_particlesPerCell);
  const std::vector<std::shared_ptr<bvh_node<point_mass> > >& leaves = m_tree.get_leaves();

  // 3. Clear particles in this cell and add new ones.
  a_particles.clear();
  for (int i = 0; i < leaves.size(); i++){
    point_mass pointMass(leaves[i]->get_data());
    ito_particle p(pointMass.mass(), pointMass.pos(), RealVect::Zero, 0.0, 0.0, pointMass.energy());
    a_particles.add(p);

#if ITO_DEBUG
    mass_after   += p.mass();
    energy_after += p.mass()*p.energy();
    center_after += p.mass()*p.position();
#endif
  }

#if ITO_DEBUG
  center_after *= 1./mass_after;

  constexpr Real BVH_ERR_TOL = 1.E-6;
  const bool break_energy = Abs(energy_before/energy_after - 1.0) > BVH_ERR_TOL;
  const bool break_center = (center_before/center_after - RealVect::Unit).vectorLength() > BVH_ERR_TOL;
  const bool break_mass   = Abs(mass_before/mass_after - 1.0) > BVH_ERR_TOL;
  if(break_mass)   pout() << "ito_solver::bvh_merge failed. Mass before = "   << mass_before   << "\t Mass after = "   << mass_after << endl;
  if(break_center) pout() << "ito_solver::bvh_merge failed. center before = " << center_before << "\t center after = " << center_after << endl;
  if(break_energy) pout() << "ito_solver::bvh_merge failed. Energy before = " << energy_before << "\t Energy after = " << energy_after << endl;
#endif
}

void ito_solver::clear(const which_container a_container){
  CH_TIME("ito_solver::clear(string)");
  if(m_verbosity > 5){
    pout() << m_name + "::clear(string)" << endl;
  }

  particle_container<ito_particle>& particles = this->get_particles(a_container);

  this->clear(particles);
}

void ito_solver::clear(particle_container<ito_particle>& a_particles){
  CH_TIME("ito_solver::clear(particle_container)");
  if(m_verbosity > 5){
    pout() << m_name + "::clear(particle_container)" << endl;
  }

  this->clear(a_particles.get_particles());
}

void ito_solver::clear(AMRParticles<ito_particle>& a_particles){
  CH_TIME("ito_solver::clear");
  if(m_verbosity > 5){
    pout() << m_name + "::clear" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    a_particles[lvl]->clear();
  }
}

RealVect ito_solver::random_position(const RealVect a_pos,
				     const RealVect a_lo,
				     const RealVect a_hi,
				     const RealVect a_bndryCentroid,
				     const RealVect a_bndryNormal,
				     const Real     a_dx,
				     const Real     a_kappa) {

  RealVect pos;
  if(a_kappa < 1.0){ // Rejection sampling. 
    pos = this->random_position(a_lo, a_hi, a_bndryCentroid, a_bndryNormal);
  }
  else{ // Regular cell. Get a position. 
    pos = this->random_position(a_lo, a_hi);
  }

  pos = a_pos + pos*a_dx;

  return pos;
}

RealVect ito_solver::random_position(const RealVect a_lo,
				     const RealVect a_hi,
				     const RealVect a_bndryCentroid,
				     const RealVect a_bndryNormal) {
  RealVect pos = this->random_position(a_lo, a_hi);
  bool valid   = PolyGeom::dot(pos-a_bndryCentroid, a_bndryNormal) >= 0.0;

  while(!valid){
    pos    = this->random_position(a_lo, a_hi);
    valid = PolyGeom::dot(pos-a_bndryCentroid, a_bndryNormal) >= 0.0;
  }

  return pos;
}

RealVect ito_solver::random_position(const RealVect a_lo, const RealVect a_hi) {

  RealVect pos = RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++){
    pos[dir] = a_lo[dir] + 0.5*(1.0 + m_udist11(m_rng))*(a_hi[dir] - a_lo[dir]);
  }

  return pos;
}
