/*!
  @file   air6_mc8.H
  @brief  6-species (3/3 charged/excited) and 8-photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air6_mc8.H"
#include "air6_mc8_species.H"
#include "data_ops.H"
#include "units.H" 

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include <ParmParse.H>
#include <PolyGeom.H>

std::string air6_mc8::s_bolsig_mobility = "E/N (Td)	Mobility *N (1/m/V/s)";
std::string air6_mc8::s_bolsig_diffco   = "E/N (Td)	Diffusion coefficient *N (1/m/s)";
std::string air6_mc8::s_bolsig_alpha    = "E/N (Td)	Total ionization freq. /N (m3/s)";
std::string air6_mc8::s_bolsig_eta      = "E/N (Td)	Total attachment freq. /N (m3/s)";


air6_mc8::air6_mc8() {

  instantiate_species();      

  parse_transport_file(); 
  parse_transport();
  parse_gas_params();
  parse_electron_mobility();
  parse_electron_diffco();
  parse_alpha();
  parse_eta();
  parse_photoi();
  parse_chemistry();
  parse_see();
  parse_domain_bc();

  init_rng();                 // Initialize random number generators
  
  parse_initial_particles();  // Parse initial particles

}

air6_mc8::~air6_mc8() {

}

void air6_mc8::read_file_entries(lookup_table& a_table, const std::string a_string){
  Real x, y;
  bool read_line = false;
  std::ifstream infile(m_transport_file);
  std::string line;

  while (std::getline(infile, line)){

    // Right trim string
    line.erase(line.find_last_not_of(" \n\r\t")+1);

    if(line == a_string){ // Begin reading
      read_line = true;
    }
    else if(line == "" & read_line){ // Stop reading
      read_line = false;
    }

    if(read_line){
      std::istringstream iss(line);
      if (!(iss >> x >> y)) {
    	continue;
      }
      a_table.add_entry(x, y);
    }
  }
  infile.close();
}

void air6_mc8::parse_transport_file(){
  ParmParse pp("air6_mc8");
  pp.get("transport_file",  m_transport_file);
  pp.get("uniform_tables",  m_uniform_entries);
  std::ifstream infile(m_transport_file);
  if(!infile.good()){
    MayDay::Abort("air6_mc8::parse_transport_file - could not find transport data");
  }
  else{
    infile.close();
  }
}

void air6_mc8::parse_transport(){
  ParmParse pp("air6_mc8");

  std::string str;

  pp.get("mobile_electrons", str);    m_mobile_electrons    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_diffusive_electrons = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);      m_diffusive_ions      = (str == "true") ? true : false;
  pp.get("mobile_ions", str);         m_mobile_ions         = (str == "true") ? true : false;
  
  pp.get("ion_mobility", m_ion_mobility);

  m_ion_diffusion = m_ion_mobility*(units::s_kb*m_T)/units::s_Qe;
}

void air6_mc8::parse_gas_params(){
  ParmParse pp("air6_mc8");

  // Pressure form input script
  pp.get("pressure",    m_p);
  pp.get("temperature", m_T);
  pp.get("frac_N2",     m_N2frac);
  pp.get("frac_O2",     m_O2frac);

  m_p  *= units::s_atm2pascal;
  m_N = m_p*units::s_Na/(m_T*units::s_R);
}

void air6_mc8::parse_electron_mobility(){
  ParmParse pp("air6_mc8");

  read_file_entries(m_e_mobility, air6_mc8::s_bolsig_mobility);
  m_e_mobility.scale_x(m_N*units::s_Td);
  m_e_mobility.scale_y(1./m_N); 
  m_e_mobility.make_uniform(m_uniform_entries);
}

void air6_mc8::parse_electron_diffco(){
  ParmParse pp("air6_mc8");
  
  read_file_entries(m_e_diffco, air6_mc8::s_bolsig_diffco);
  m_e_diffco.scale_x(m_N*units::s_Td);
  m_e_diffco.scale_y(1./m_N); 
  m_e_diffco.make_uniform(m_uniform_entries);
}

void air6_mc8::parse_alpha(){
  ParmParse pp("air6_mc8");
  read_file_entries(m_e_alpha, air6_mc8::s_bolsig_alpha);
  m_e_alpha.scale_x(m_N*units::s_Td);
  m_e_alpha.scale_y(m_N); 
  m_e_alpha.make_uniform(m_uniform_entries);
}

void air6_mc8::parse_eta(){
  ParmParse pp("air6_mc8");
  read_file_entries(m_e_eta, air6_mc8::s_bolsig_eta);
  m_e_eta.scale_x(m_N*units::s_Td);
  m_e_eta.scale_y(m_N); 
  m_e_eta.make_uniform(m_uniform_entries);
}

void air6_mc8::parse_photoi(){

  ParmParse pp("air6_mc8");

  pp.get("c4v0_exc_eff", m_c4v0_exc_eff);
  pp.get("c4v1_exc_eff", m_c4v1_exc_eff);
  pp.get("b1v1_exc_eff", m_b1v1_exc_eff);

  pp.get("c4v0_X1v0_photoi_eff", m_c4v0_X1v0_photoi_eff);
  pp.get("c4v0_X1v1_photoi_eff", m_c4v0_X1v1_photoi_eff);
  pp.get("c4v1_X1v0_photoi_eff", m_c4v1_X1v0_photoi_eff);
  pp.get("c4v1_X1v1_photoi_eff", m_c4v1_X1v1_photoi_eff);
  pp.get("c4v1_X1v2_photoi_eff", m_c4v1_X1v2_photoi_eff);
  pp.get("c4v1_X1v3_photoi_eff", m_c4v1_X1v3_photoi_eff);
  pp.get("b1v1_X1v0_photoi_eff", m_b1v1_X1v0_photoi_eff);
  pp.get("b1v1_X1v1_photoi_eff", m_b1v1_X1v1_photoi_eff);

  pp.get("c4v0_X1v0_tau_rad", m_c4v0_X1v0_tau_r);
  pp.get("c4v0_X1v1_tau_rad", m_c4v0_X1v1_tau_r);
  pp.get("c4v1_X1v0_tau_rad", m_c4v1_X1v0_tau_r);
  pp.get("c4v1_X1v1_tau_rad", m_c4v1_X1v1_tau_r);
  pp.get("c4v1_X1v2_tau_rad", m_c4v1_X1v2_tau_r);
  pp.get("c4v1_X1v3_tau_rad", m_c4v1_X1v3_tau_r);
  pp.get("b1v1_X1v0_tau_rad", m_b1v1_X1v0_tau_r);
  pp.get("b1v1_X1v1_tau_rad", m_b1v1_X1v1_tau_r);

  pp.get("c4v0_tau_predissoc", m_c4v0_tau_p);
  pp.get("c4v1_tau_predissoc", m_c4v1_tau_p);
  pp.get("b1v1_tau_predissoc", m_b1v1_tau_p);

  pp.get("quenching_pressure", m_pq);

  m_pq *= units::s_atm2pascal;

  // Set all quenching lifetimes to radiative lifetime x p/pq
  m_c4v0_X1v0_tau_q = m_c4v0_X1v0_tau_r*m_p/m_pq;
  m_c4v0_X1v1_tau_q = m_c4v0_X1v1_tau_r*m_p/m_pq;
  m_c4v1_X1v0_tau_q = m_c4v1_X1v0_tau_r*m_p/m_pq;
  m_c4v1_X1v1_tau_q = m_c4v1_X1v1_tau_r*m_p/m_pq;
  m_c4v1_X1v2_tau_q = m_c4v1_X1v2_tau_r*m_p/m_pq;
  m_c4v1_X1v3_tau_q = m_c4v1_X1v3_tau_r*m_p/m_pq;
  m_b1v1_X1v0_tau_q = m_b1v1_X1v0_tau_r*m_p/m_pq;
  m_b1v1_X1v1_tau_q = m_b1v1_X1v1_tau_r*m_p/m_pq;
}

void air6_mc8::parse_chemistry(){
  ParmParse pp("air6_mc8");
  std::string str;

  pp.get("chemistry", str);
  if(str == "ssa"){
    m_scomp = source_comp::ssa;
  }
  else if(str == "tau"){
    m_scomp = source_comp::tau;
  }
  else if(str == "rre"){
    m_scomp = source_comp::rre;
  }
  else{
    MayDay::Abort("air6_mc8::parse_reaction_settings - stop!");
  }

  pp.get("poisson_exp_swap", m_poiss_exp_swap);
}

void air6_mc8::parse_see(){
  ParmParse pp("air6_mc8");
  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void air6_mc8::init_rng(){
  ParmParse pp("air6_mc8");
  pp.get("rng_seed", m_rng_seed);
  
  if(m_rng_seed < 0) {
    m_rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  m_udist11 = new std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_rng     = new std::mt19937_64(m_rng_seed);
}

void air6_mc8::instantiate_species(){
  m_num_species = 6;
  m_num_photons = 8;

  m_elec_idx = 0;
  m_plus_idx = 1;
  m_minu_idx = 2;
  m_c4v0_idx = 3;
  m_c4v1_idx = 4;
  m_b1v1_idx = 5;

  m_c4v0_X1v0_idx = 0;
  m_c4v0_X1v1_idx = 1;
  m_c4v1_X1v0_idx = 2;
  m_c4v1_X1v1_idx = 3;
  m_c4v1_X1v2_idx = 4;
  m_c4v1_X1v3_idx = 5;
  m_b1v1_X1v0_idx = 6;
  m_b1v1_X1v1_idx = 7;

  m_species.resize(m_num_species);
  m_species[m_elec_idx]  = RefCountedPtr<species>      (new air6_mc8::electron());
  m_species[m_plus_idx]  = RefCountedPtr<species>      (new air6_mc8::M_plus());
  m_species[m_minu_idx]  = RefCountedPtr<species>      (new air6_mc8::M_minus());
  m_species[m_c4v0_idx]  = RefCountedPtr<species>      (new air6_mc8::N2_c4v0());
  m_species[m_c4v1_idx]  = RefCountedPtr<species>      (new air6_mc8::N2_c4v1());
  m_species[m_b1v1_idx]  = RefCountedPtr<species>      (new air6_mc8::N2_b1v1());

  m_photons.resize(m_num_photons);
  m_photons[m_c4v0_X1v0_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v0_X1v0());
  m_photons[m_c4v0_X1v1_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v0_X1v1());
  m_photons[m_c4v1_X1v0_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v1_X1v0());
  m_photons[m_c4v1_X1v1_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v1_X1v1());
  m_photons[m_c4v1_X1v2_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v1_X1v2());
  m_photons[m_c4v1_X1v3_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v1_X1v3());
  m_photons[m_b1v1_X1v0_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_b1v1_X1v0());
  m_photons[m_b1v1_X1v1_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_b1v1_X1v1());
}

void air6_mc8::parse_initial_particles(){



  // Get some parameters from the input script
  Vector<Real> vec;
  Real weight;
  Real uniform_pairs, uniform_excited;
  Real rad_pairs, rad_excited;
  Real gaussian_pairs, gaussian_excited;
  RealVect center_pairs, center_excited;

  ParmParse pp("air6_mc8");
  pp.get("particle_weight",            weight);
  pp.get("uniform_pairs",              uniform_pairs);
  pp.get("uniform_excited",            uniform_excited);
  pp.get("gaussian_pairs",             gaussian_pairs);
  pp.get("gaussian_excited",           gaussian_excited);
  pp.get("gaussian_pairs_radius",      rad_pairs);
  pp.get("gaussian_excited_radius",    rad_excited);
  pp.getarr("gaussian_pairs_center",   vec, 0, SpaceDim); center_pairs   = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("gaussian_excited_center", vec, 0, SpaceDim); center_excited = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  
  List<Particle> electron_ion_pairs;
  List<Particle> excited_molecules;


  // Add various types of particles
  //  add_uniform_particles(electron_ion_pairs,  round(uniform_pairs),    weight);
  //  add_uniform_particles(excited_molecules,   round(uniform_excited),  weight);
  add_gaussian_particles(electron_ion_pairs, round(gaussian_pairs),   weight, rad_pairs,   center_pairs);
  add_gaussian_particles(excited_molecules,  round(gaussian_excited), weight, rad_excited, center_excited);
  
  // Set initial particles
  m_species[m_elec_idx]->get_initial_particles() = electron_ion_pairs;
  m_species[m_plus_idx]->get_initial_particles() = electron_ion_pairs;
  m_species[m_c4v0_idx]->get_initial_particles() = excited_molecules;
  m_species[m_c4v1_idx]->get_initial_particles() = excited_molecules;
  m_species[m_b1v1_idx]->get_initial_particles() = excited_molecules;

  // Set the deposition scheme
  std::string str;
  InterpType deposition;
  pp.get("particle_deposition", str);
  if(str == "cic"){
    deposition = InterpType::CIC;
  }
  else if(str == "ngp"){
    deposition = InterpType::NGP;
  }
  else if(str == "tsc"){
    deposition = InterpType::TSC;
  }
  else{
    MayDay::Abort("air6_mc8::parse_initial_particles - unknown deposition type requested");
  }
  for (int i = 0; i < m_num_species; i++){
    m_species[i]->get_deposition() = deposition;
  }
}

void air6_mc8::add_uniform_particles(List<Particle>& a_particles, const int a_num, const Real a_weight){

  // Get Lo/Hi sides of domain
  RealVect lo, hi;
  Vector<Real> vec(SpaceDim);
  ParmParse pp("physical_domain");
  pp.getarr("lo_corner", vec, 0, SpaceDim); lo = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("hi_corner", vec, 0, SpaceDim); hi = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  auto rngX = std::uniform_real_distribution<Real>(lo[0], hi[0]);
  auto rngY = std::uniform_real_distribution<Real>(lo[1], hi[1]);
#if CH_SPACEDIM==3
  auto rngZ = std::uniform_real_distribution<Real>(lo[2], hi[2]);
#endif

  for (int i = 0; i < a_num; i++){
    const Real x = rngX(*m_rng);
    const Real y = rngX(*m_rng);
#if CH_SPACEDIM==3
    const Real z = rngZ(*m_rng);
#endif
    RealVect pos = RealVect(D_DECL(x, y, z));
    a_particles.add(Particle(a_weight, pos));
  }
}

void air6_mc8::add_gaussian_particles(List<Particle>& a_particles,
				      const int       a_num,
				      const Real      a_weight,
				      const Real      a_rad,
				      const RealVect  a_center){
  m_gauss = std::normal_distribution<Real>(0., a_rad);

  for (int i = 0; i < a_num; i++){
    RealVect pos = a_center + random_gaussian();
    a_particles.add(Particle(a_weight, pos));
  }
}

Real air6_mc8::ran01() const{
  return (*m_udist01)(*m_rng);
}

Real air6_mc8::ran11() const{
  return (*m_udist11)(*m_rng);
}

RealVect air6_mc8::random_gaussian(){

  const Real rad = m_gauss(*m_rng);
  return rad*random_direction();
}

RealVect air6_mc8::random_direction(){
#if CH_SPACEDIM == 2
  return random_direction2D();
#else
  return random_direction3D();
#endif
}

#if CH_SPACEDIM == 2
RealVect air6_mc8::random_direction2D(){
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
RealVect air6_mc8::random_direction3D(){
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

void air6_mc8::parse_domain_bc(){

  ParmParse pp("air6_mc8");
  std::string str;

  m_wallbc.resize(2*SpaceDim, 0); 
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
	
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

      // Check for wall BCs
      if(side == Side::Lo){
	std::string type;
	std::string bc_string = "domain_bc_" + str_dir + "_lo";
	if(pp.contains(bc_string.c_str())){
	  pp.get(bc_string.c_str(), type);
	  const int idx = 2*dir;
	  if(type == "wall"){
	    m_wallbc[idx] = 1;
	  }
	}
      }
      else if(side == Side::Hi){
	std::string type;
	std::string bc_string = "domain_bc_" + str_dir + "_hi";
	if(pp.contains(bc_string.c_str())){
	  pp.get(bc_string.c_str(), type);
	  const int idx = 2*dir + 1;
	  if(type == "wall"){
	    m_wallbc[idx] = 1;
	  }
	}
      }
    }
  }
}

void air6_mc8::advance_reaction_network(Vector<Real>&          a_particle_sources,
					Vector<Real>&          a_photon_sources,
					const Vector<Real>     a_particle_densities,
					const Vector<RealVect> a_particle_gradients,
					const Vector<Real>     a_photon_densities,
					const RealVect         a_E,
					const RealVect         a_pos,
					const Real             a_dx,
					const Real             a_dt,
					const Real             a_time,
					const Real             a_kappa) const{

  
  // The various forms of the chemistry step is given in the routines below
  if(m_scomp == source_comp::ssa){ // SSA algorithm
    network_ssa(a_particle_sources, a_photon_sources, a_particle_densities, a_particle_gradients, a_photon_densities,
		a_E, a_pos, a_dx, a_dt, a_time, a_kappa);
  }
  else if(m_scomp == source_comp::tau){ // Tau leaping
    network_tau(a_particle_sources, a_photon_sources, a_particle_densities, a_particle_gradients, a_photon_densities,
		a_E, a_pos, a_dx, a_dt, a_time, a_kappa);
  }
  else if(m_scomp == source_comp::rre){ // Reaction rate equation
    network_rre(a_particle_sources, a_photon_sources, a_particle_densities, a_particle_gradients, a_photon_densities,
		a_E, a_pos, a_dx, a_dt, a_time, a_kappa);
  }

  return;
}


Vector<Real> air6_mc8::compute_cdr_diffusion_coefficients(const Real         a_time,
							  const RealVect     a_pos,
							  const RealVect     a_E,
							  const Vector<Real> a_cdr_densities) const {

  Vector<Real> dco(m_num_species, 0.0);
  dco[m_elec_idx] = m_e_diffco.get_entry(a_E.vectorLength());
  dco[m_plus_idx] = m_ion_diffusion;
  dco[m_minu_idx] = m_ion_diffusion;
  
  return dco;

}
  
Vector<RealVect> air6_mc8::compute_cdr_velocities(const Real         a_time,
						  const RealVect     a_pos,
						  const RealVect     a_E,
						  const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_num_species, RealVect::Zero);

  vel[m_elec_idx] = -a_E*m_e_mobility.get_entry(a_E.vectorLength());
  vel[m_plus_idx] =  a_E*m_ion_mobility;
  vel[m_minu_idx] = -a_E*m_ion_mobility;
  
  return vel;
}
  
Vector<Real> air6_mc8::compute_cdr_domain_fluxes(const Real           a_time,
						 const RealVect       a_pos,
						 const int            a_dir,
						 const Side::LoHiSide a_side,
						 const RealVect       a_E,
						 const Vector<Real>   a_cdr_densities,
						 const Vector<Real>   a_cdr_velocities,
						 const Vector<Real>   a_cdr_gradients,
						 const Vector<Real>   a_rte_fluxes,
						 const Vector<Real>   a_extrap_cdr_fluxes) const{
  Vector<Real> fluxes(m_num_species, 0.0);

  int idx, sgn;
  if(a_side == Side::Lo){
    sgn = -1;
    idx = 2*a_dir;
  }
  else{
    sgn = 1;
    idx = 2*a_dir + 1;
  }

  if(m_wallbc[idx] == 0){ // Inflow/outflow
    for (int i = 0; i < fluxes.size(); i++){
      fluxes[i] = a_extrap_cdr_fluxes[i];
    }
  }
  else if(m_wallbc[idx] == 1){ // wall
    for (int i = 0; i < fluxes.size(); i++){
      fluxes[i] = 0.0;
    }
  }
  else{
    MayDay::Abort("morrow_jiang::compute_cdr_domain_fluxes - uknown domain bc requested");
  }

  
  return fluxes;
}
  
Vector<Real> air6_mc8::compute_cdr_electrode_fluxes(const Real         a_time,
						    const RealVect     a_pos,
						    const RealVect     a_normal,
						    const RealVect     a_E,
						    const Vector<Real> a_cdr_densities,
						    const Vector<Real> a_cdr_velocities,
						    const Vector<Real> a_cdr_gradients,
						    const Vector<Real> a_rte_fluxes,
						    const Vector<Real> a_extrap_cdr_fluxes) const{
  return compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			    a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
}

Vector<Real> air6_mc8::compute_cdr_dielectric_fluxes(const Real         a_time,
						     const RealVect     a_pos,
						     const RealVect     a_normal,
						     const RealVect     a_E,
						     const Vector<Real> a_cdr_densities,
						     const Vector<Real> a_cdr_velocities,
						     const Vector<Real> a_cdr_gradients,
						     const Vector<Real> a_rte_fluxes,
						     const Vector<Real> a_extrap_cdr_fluxes) const{
  return compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			    a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

Vector<Real> air6_mc8::compute_cdr_fluxes(const Real         a_time,
					  const RealVect     a_pos,
					  const RealVect     a_normal,
					  const RealVect     a_E,
					  const Vector<Real> a_cdr_densities,
					  const Vector<Real> a_cdr_velocities,
					  const Vector<Real> a_cdr_gradients,
					  const Vector<Real> a_rte_fluxes,
					  const Vector<Real> a_extrap_cdr_fluxes,
					  const Real         a_townsend2,
					  const Real         a_quantum_efficiency) const{
  Vector<Real> fluxes(m_num_species, 0.0);

  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;

  // Switch for setting drift flux to zero for charge species
  Vector<Real> aj(m_num_species, 0.0);
  for (int i = 0; i < m_num_species; i++){
    if(data_ops::sgn(m_species[i]->get_charge())*PolyGeom::dot(a_E, a_normal) < 0){
      aj[i] = 1.0;
    }
    else {
      aj[i] = 0.0;
    }
  }

  // Drift outflow for now
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  return fluxes;
}

Real air6_mc8::initial_sigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

int air6_mc8::poisson_reaction(const Real a_propensity, const Real a_dt) const{
  int value = 0;
  const Real mean = a_propensity*a_dt;

  if(mean < m_poiss_exp_swap){
    std::poisson_distribution<int> dist(mean);
    value = dist(*m_rng);
  }
  else{
    std::normal_distribution<double> dist(mean, sqrt(mean));
    value = dist(*m_rng);
  }

  return value;
}

int air6_mc8::binomial_trials(const int a_trials, const Real a_p) const {
  std::binomial_distribution<int> dist(a_trials, a_p);
  return dist(*m_rng);
}

void air6_mc8::network_ssa(Vector<Real>&          a_particle_sources,
			   Vector<Real>&          a_photon_sources,
			   const Vector<Real>     a_particle_densities,
			   const Vector<RealVect> a_particle_gradients,
			   const Vector<Real>     a_photon_densities,
			   const RealVect         a_E,
			   const RealVect         a_pos,
			   const Real             a_dx,
			   const Real             a_dt,
			   const Real             a_time,
			   const Real             a_kappa) const{
  MayDay::Abort("air6_mc8::network_ssa - not implemented (yet)");
}

void air6_mc8::network_tau(Vector<Real>&          a_particle_sources,
			   Vector<Real>&          a_photon_sources,
			   const Vector<Real>     a_particle_densities,
			   const Vector<RealVect> a_particle_gradients,
			   const Vector<Real>     a_photon_densities,
			   const RealVect         a_E,
			   const RealVect         a_pos,
			   const Real             a_dx,
			   const Real             a_dt,
			   const Real             a_time,
			   const Real             a_kappa) const{

  // Reactions that we track
  //
  // R1: e + M -> e + e + M+  alpha*Xe
  // R2: e + M -> M-+         eta*Xe
  // R3: e + M -> c4v0        alpha*Xe*exc_eff(c4v0)
  // R4: e + M -> c4v1        alpha*Xe*exc_eff(c4v1)
  // R5: e + M -> b1v1        alpha*Xe*exc_eff(b1v1)


  const Real volume = pow(a_dx, SpaceDim);

  Vector<Real> x(m_num_species, 0.0);
  Vector<int> X(m_num_species, 0);
  Vector<int> Y(m_num_photons, 0);

  const Real thresh = 1.E-2;
  
  for (int i = 0; i < m_num_species; i++){
    X[i] = floor(thresh + a_particle_densities[i]*volume); // Integer particles
    x[i] = a_particle_densities[i] - X[i]*volume;          // "Partial particles"
  }

  for (int i = 0; i < m_num_photons; i++){
    Y[i] = floor(thresh + a_photon_densities[i]);
  }

  // Aux functions
  const Real alpha  = m_e_alpha.get_entry(a_E.vectorLength());
  const Real eta    = m_e_eta.get_entry(a_E.vectorLength());

  const Real x_elec = x[m_elec_idx];
  const Real x_plus = x[m_plus_idx];
  const Real x_minu = x[m_minu_idx];
  const Real x_c4v0 = x[m_c4v0_idx];
  const Real x_c4v1 = x[m_c4v1_idx];
  const Real x_b1v1 = x[m_b1v1_idx];

  // Deposit photons first. These are allowed to react
  for (int i = 0; i < m_num_photons; i++){
    X[m_elec_idx] += Y[i];
    X[m_plus_idx] += Y[i];
  }
  
  const int X_elec = X[m_elec_idx];
  const int X_plus = X[m_plus_idx];
  const int X_minu = X[m_minu_idx];
  const int X_c4v0 = X[m_c4v0_idx];
  const int X_c4v1 = X[m_c4v1_idx];
  const int X_b1v1 = X[m_b1v1_idx];


  // Integer number of each species that are evolved
  int s_elec = 0;
  int s_plus = 0;
  int s_minu = 0;
  int s_c4v0 = 0;
  int s_c4v1 = 0;
  int s_b1v1 = 0;


  // e + M -> e + e + M+
  const Real A1 = X_elec*alpha;
  const Real a1 = x_elec*alpha;
  const int  S1 = poisson_reaction(A1, a_dt);
  s_elec += S1;
  s_plus += S1;

  // e + M -> M-
  const Real A2 = X_elec*eta;
  const Real a2 = x_elec*eta;
  const int  S2 = poisson_reaction(A2, a_dt);
  s_elec -= S2;
  s_minu += S2;

  // e + M -> e + M(c4v0)
  const Real A3 = A1*m_c4v0_exc_eff;
  const Real a3 = a1*m_c4v0_exc_eff;
  const int  S3 = poisson_reaction(A3, a_dt);
  s_c4v0 += S3;

  // e + M -> e + M(c4v1)
  const Real A4 = A1*m_c4v1_exc_eff;
  const Real a4 = A1*m_c4v1_exc_eff;
  const int  S4 = poisson_reaction(A4, a_dt);
  s_c4v1 += S4;

  // e + M -> e + M(b1v1)
  const Real A5 = A1*m_b1v1_exc_eff;
  const Real a5 = a1*m_b1v1_exc_eff;
  const int  S5 = poisson_reaction(A5, a_dt);
  s_b1v1 += S5;

  // M(c4v0) -> M through predissociation
  const Real A6 = X_c4v0/m_c4v0_tau_p;
  const Real a6 = x_c4v0/m_c4v0_tau_p;
  const int  S6 = poisson_reaction(A6, a_dt);
  s_c4v0 -= S6;

  // M(c4v1) -> M through predissociation
  const Real A7 = X_c4v1/m_c4v1_tau_p;
  const Real a7 = x_c4v1/m_c4v1_tau_p;
  const int  S7 = poisson_reaction(A7, a_dt);
  s_c4v1 -= S7;

  // M(b1v1) -> M through predissociation
  const Real A8 = X_b1v1/m_b1v1_tau_p;
  const Real a8 = x_b1v1/m_b1v1_tau_p;
  const int  S8 = poisson_reaction(A8, a_dt);
  s_b1v1 -= S8;

  // M(c4v0) -> M + y(c4v0->X1v0)
  const Real A9 = X_c4v0/m_c4v0_X1v0_tau_r;
  const Real a9 = x_c4v0/m_c4v0_X1v0_tau_r;
  const int  S9 = poisson_reaction(A9, a_dt);
  const int SS9 = binomial_trials(S9, m_c4v0_X1v0_photoi_eff);
  s_c4v0 -= S9;
  a_photon_sources[m_c4v0_X1v0_idx] = SS9;

  // M(c4v0) -> M + through quenching
  const Real A10 = X_c4v0/m_c4v0_X1v0_tau_q;
  const Real a10 = x_c4v0/m_c4v0_X1v0_tau_q;
  const int  S10 = poisson_reaction(A10, a_dt);
  s_c4v0 -= S10;

  // M(c4v0) -> M + y(c4v0->X1v1)
  const Real A11 = X_c4v0/m_c4v0_X1v1_tau_r;
  const Real a11 = x_c4v0/m_c4v0_X1v1_tau_r;
  const int  S11 = poisson_reaction(A11, a_dt);
  const int SS11 = binomial_trials(S11, m_c4v0_X1v1_photoi_eff);
  s_c4v0 -= S11;
  a_photon_sources[m_c4v0_X1v1_idx] = 1.0*SS11;

  // M(c4v0) -> M through quenching
  const Real A12 = X_c4v0/m_c4v0_X1v1_tau_q;
  const Real a12 = x_c4v0/m_c4v0_X1v1_tau_q;
  const int  S12 = poisson_reaction(A12, a_dt);
  s_c4v0 -= S12;
  
  // M(c4v1) -> M + y(c4v1->X1v0)
  const Real A13 = X_c4v1/m_c4v1_X1v0_tau_r;
  const Real a13 = x_c4v1/m_c4v1_X1v0_tau_r;
  const int  S13 = poisson_reaction(A13, a_dt);
  const int SS13 = binomial_trials(S13, m_c4v1_X1v0_photoi_eff);
  s_c4v1 -= S13;
  a_photon_sources[m_c4v1_X1v0_idx] = 1.0*SS13;

  // M(c4v0) -> M + through quenching
  const Real A14 = X_c4v1/m_c4v1_X1v0_tau_q;
  const Real a14 = x_c4v1/m_c4v1_X1v0_tau_q;
  const int  S14 = poisson_reaction(A14, a_dt);
  s_c4v1 -= S14;

  // M(c4v1) -> M + y(c4v1->X1v1)
  const Real A15 = X_c4v1/m_c4v1_X1v1_tau_r;
  const Real a15 = x_c4v1/m_c4v1_X1v1_tau_r;
  const int  S15 = poisson_reaction(A15, a_dt);
  const int SS15 = binomial_trials(S15, m_c4v1_X1v1_photoi_eff);
  s_c4v1 -= S15;
  a_photon_sources[m_c4v1_X1v1_idx] = 1.0*SS15;

  // M(c4v0) -> M + through quenching
  const Real A16 = X_c4v1/m_c4v1_X1v1_tau_q;
  const Real a16 = x_c4v1/m_c4v1_X1v1_tau_q;
  const int  S16 = poisson_reaction(A16, a_dt);
  s_c4v1 -= S16;

  // M(c4v1) -> M + y(c4v1->X1v2)
  const Real A17 = X_c4v1/m_c4v1_X1v2_tau_r;
  const Real a17 = x_c4v1/m_c4v1_X1v2_tau_r;
  const int  S17 = poisson_reaction(A17, a_dt);
  const int SS17 = binomial_trials(S17, m_c4v1_X1v2_photoi_eff);
  s_c4v1 -= S17;
  a_photon_sources[m_c4v1_X1v2_idx] = 1.0*SS17;

  // M(c4v1) -> M + through quenching
  const Real A18 = X_c4v1/m_c4v1_X1v2_tau_q;
  const Real a18 = x_c4v1/m_c4v1_X1v2_tau_q;
  const int  S18 = poisson_reaction(A18, a_dt);
  s_c4v1 -= S18;

  // M(c4v1) -> M + y(c4v1->X1v3)
  const Real A19 = X_c4v1/m_c4v1_X1v3_tau_r;
  const Real a19 = x_c4v1/m_c4v1_X1v3_tau_r;
  const int  S19 = poisson_reaction(A19, a_dt);
  const int SS19 = binomial_trials(S19, m_c4v1_X1v3_photoi_eff);
  s_c4v1 -= S19;
  a_photon_sources[m_c4v1_X1v3_idx] = 1.0*SS19;

  // M(c4v1) -> M + through quenching
  const Real A20 = X_c4v1/m_c4v1_X1v3_tau_q;
  const Real a20 = x_c4v1/m_c4v1_X1v3_tau_q;
  const int  S20 = poisson_reaction(A20, a_dt);
  s_c4v1 -= S20;

  // M(b1v1) -> M + y(b1v1->X1v0)
  const Real A21 = X_b1v1/m_b1v1_X1v0_tau_r;
  const Real a21 = x_b1v1/m_b1v1_X1v0_tau_r;
  const int  S21 = poisson_reaction(A21, a_dt);
  const int SS21 = binomial_trials(S21, m_b1v1_X1v0_photoi_eff);
  s_b1v1 -= S21;
  a_photon_sources[m_b1v1_X1v0_idx] = 1.0*SS21;

  // M(c4v1) -> M + through quenching
  const Real A22 = X_b1v1/m_b1v1_X1v0_tau_q;
  const Real a22 = x_b1v1/m_b1v1_X1v0_tau_q;
  const int  S22 = poisson_reaction(A22, a_dt);
  s_b1v1 -= S22;
  
  // M(b1v1) -> M + y(b1v1->X1v1)
  const Real A23 = X_b1v1/m_b1v1_X1v1_tau_r;
  const Real a23 = x_b1v1/m_b1v1_X1v1_tau_r;
  const int  S23 = poisson_reaction(A23, a_dt);
  const int SS23 = binomial_trials(S23, m_b1v1_X1v1_photoi_eff);
  s_b1v1 -= S23;
  a_photon_sources[m_b1v1_X1v1_idx] = 1.0*SS23;

  // M(c4v1) -> M + through quenching
  const Real A24 = X_b1v1/m_b1v1_X1v1_tau_q;
  const Real a24 = x_b1v1/m_b1v1_X1v1_tau_q;
  const int  S24 = poisson_reaction(A24, a_dt);
  s_b1v1 -= S24;


  // Normalize sources
  const Real factor = 1./(volume*a_dt);
  a_particle_sources[m_elec_idx] = s_elec*factor;
  a_particle_sources[m_plus_idx] = s_plus*factor;
  a_particle_sources[m_minu_idx] = s_minu*factor;
  a_particle_sources[m_c4v0_idx] = s_c4v0*factor;
  a_particle_sources[m_c4v1_idx] = s_c4v1*factor;
  a_particle_sources[m_b1v1_idx] = s_b1v1*factor;

  return;
  // "Missing reactions" based on RRE approximation
#if 1 // Original code
  a_particle_sources[m_elec_idx] += a1 - a2;
  a_particle_sources[m_plus_idx] += a1;
  a_particle_sources[m_minu_idx] += a2;
#else // Debug code
  a_particle_sources[m_elec_idx] = (alpha-eta)*a_particle_densities[m_elec_idx];
  a_particle_sources[m_plus_idx] = alpha*a_particle_densities[m_elec_idx];
  a_particle_sources[m_minu_idx] = eta*a_particle_densities[m_elec_idx];
#endif
  
}


void air6_mc8::network_rre(Vector<Real>&          a_particle_sources,
			   Vector<Real>&          a_photon_sources,
			   const Vector<Real>     a_particle_densities,
			   const Vector<RealVect> a_particle_gradients,
			   const Vector<Real>     a_photon_densities,
			   const RealVect         a_E,
			   const RealVect         a_pos,
			   const Real             a_dx,
			   const Real             a_dt,
			   const Real             a_time,
			   const Real             a_kappa) const{
  MayDay::Abort("air6_mc8::network_rre - not implemented (yet)");
}
