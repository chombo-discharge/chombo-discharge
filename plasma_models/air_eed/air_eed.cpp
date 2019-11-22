/*!
  @file   air_eed.H
  @brief  6-species (3/3 charged/excited) and 8-photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air_eed.H"
#include "air_eed_species.H"
#include "data_ops.H"
#include "units.H" 

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <random>

#include <ParmParse.H>
#include <PolyGeom.H>

// Bulk stuff
std::string air_eed::s_bolsig_mobility  = "Energy (eV)	Mobility *N (1/m/V/s)";
std::string air_eed::s_bolsig_diffco    = "Energy (eV)	Diffusion coefficient *N (1/m/s)";
std::string air_eed::s_bolsig_alpha     = "Energy (eV)	Townsend ioniz. coef. alpha/N (m2)";
std::string air_eed::s_bolsig_eta       = "Energy (eV)	Townsend attach. coef. eta/N (m2)";
std::string air_eed::s_bolsig_elastic   = "Energy (eV)	Elastic power loss /N (eV m3/s)";
std::string air_eed::s_bolsig_inelastic = "Energy (eV)	Inelastic power loss /N (eV m3/s)";
std::string air_eed::s_bolsig_lfa       = "Energy (eV) 	Electric field / N (Td)";


air_eed::air_eed() {
  instantiate_species();
  parse_transport_file();
  parse_transport();
  parse_gas_params();
  parse_electron_mobility();
  parse_electron_diffco();
  parse_alpha();
  parse_eta();
  parse_losses();
  parse_lfa();
  parse_photoi();
  parse_see();
  parse_domain_bc();
  parse_chemistry();

  init_rng();                 // Initialize random number generators
  
  parse_initial_particles();  // Parse initial particles

  //  MayDay::Abort("stop");
}

air_eed::~air_eed() {

}

void air_eed::read_file_entries(lookup_table& a_table, const std::string a_string){
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

void air_eed::parse_transport_file(){
  ParmParse pp("air_eed");
  pp.get("transport_file",  m_transport_file);
  pp.get("uniform_tables",  m_uniform_entries);
  pp.get("min_energy",      m_min_energy);
  pp.get("max_energy",      m_max_energy);
  std::ifstream infile(m_transport_file);
  if(!infile.good()){
    MayDay::Abort("air_eed::parse_transport_file - could not find transport data");
  }
  else{
    infile.close();
  }
}

void air_eed::parse_transport(){
  ParmParse pp("air_eed");

  std::string str;

  pp.get("mobile_electrons", str);    m_mobile_electrons    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_diffusive_electrons = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);      m_diffusive_ions      = (str == "true") ? true : false;
  pp.get("mobile_ions", str);         m_mobile_ions         = (str == "true") ? true : false;
  pp.get("ion_mobility", m_ion_mobility);
  pp.get("photoi_gain", m_photoi_gain);

  m_ion_diffusion = m_ion_mobility*(units::s_kb*m_T)/units::s_Qe;
}

void air_eed::parse_gas_params(){
  ParmParse pp("air_eed");

  // Pressure form input script
  pp.get("pressure",    m_p);
  pp.get("temperature", m_T);
  pp.get("frac_N2",     m_N2frac);
  pp.get("frac_O2",     m_O2frac);

  m_p  *= units::s_atm2pascal;
  m_N = m_p*units::s_Na/(m_T*units::s_R);
}

void air_eed::parse_electron_mobility(){
  ParmParse pp("air_eed");

  read_file_entries(m_e_mobility, air_eed::s_bolsig_mobility);
  m_e_mobility.scale_y(1./m_N); 
  m_e_mobility.make_uniform(m_uniform_entries);
}

void air_eed::parse_electron_diffco(){
  ParmParse pp("air_eed");
  
  read_file_entries(m_e_diffco, air_eed::s_bolsig_diffco);
  m_e_diffco.scale_y(1./m_N);
  m_e_diffco.make_uniform(m_uniform_entries);
}

void air_eed::parse_alpha(){
  ParmParse pp("air_eed");
  read_file_entries(m_e_alpha, air_eed::s_bolsig_alpha);
  m_e_alpha.scale_y(m_N); 
  m_e_alpha.make_uniform(m_uniform_entries);
}

void air_eed::parse_eta(){
  ParmParse pp("air_eed");
  read_file_entries(m_e_eta, air_eed::s_bolsig_eta);

  m_e_eta.scale_y(m_N); 
  m_e_eta.make_uniform(m_uniform_entries);
}

void air_eed::parse_losses(){
  ParmParse pp("air_eed");

  // First read elastic losses, then inelastic losses
  lookup_table elastic_losses;
  lookup_table inelastic_losses;
  read_file_entries(elastic_losses,   air_eed::s_bolsig_elastic);
  read_file_entries(inelastic_losses, air_eed::s_bolsig_inelastic);

  m_e_losses = elastic_losses;
  m_e_losses.add_table(inelastic_losses, 1.0);
  m_e_losses.scale_y(m_N);
  m_e_losses.make_uniform(m_uniform_entries);

  m_ionization_loss = m_N2frac*15.6 + m_O2frac*12.06;
}

void air_eed::parse_lfa(){
  ParmParse pp("air_eed");
  read_file_entries(m_lfa, air_eed::s_bolsig_lfa);

  m_lfa.scale_y(m_N*units::s_Td);
  m_lfa.swap_xy();
  m_lfa.make_uniform(m_uniform_entries);

  m_lfa.dump_table();
}

void air_eed::parse_photoi(){

  ParmParse pp("air_eed");

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

  pp.get("quenching_pressure", m_pq);

  m_pq *= units::s_atm2pascal;
}

void air_eed::parse_see(){
  ParmParse pp("air_eed");
  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void air_eed::init_rng(){
  ParmParse pp("air_eed");
  pp.get("rng_seed", m_rng_seed);

  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  if(m_rng_seed < 0) {
    m_rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  m_udist11 = new std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_rng     = new std::mt19937_64(m_rng_seed);
}

void air_eed::instantiate_species(){
  m_num_species = 4;
  m_num_photons = 1;

  m_eed_idx  = 0;
  m_elec_idx = 1;
  m_plus_idx = 2;
  m_minu_idx = 3;

  m_photon_idx = 0;

  m_species.resize(m_num_species);
  m_species[m_eed_idx]   = RefCountedPtr<species>      (new air_eed::eed());
  m_species[m_elec_idx]  = RefCountedPtr<species>      (new air_eed::electron());
  m_species[m_plus_idx]  = RefCountedPtr<species>      (new air_eed::M_plus());
  m_species[m_minu_idx]  = RefCountedPtr<species>      (new air_eed::M_minus());

  m_photons.resize(m_num_photons);
  m_photons[m_photon_idx] = RefCountedPtr<photon_group> (new air_eed::agg_photon());
}

void air_eed::parse_initial_particles(){

  // Get some parameters from the input script
  Vector<Real> vec;
  Real weight, uniform_pairs, rad_pairs, gaussian_pairs, init_energy;
  RealVect center_pairs;

  ParmParse pp("air_eed");
  pp.get("particle_weight",            weight);
  pp.get("uniform_pairs",              uniform_pairs);
  pp.get("gaussian_pairs",             gaussian_pairs);
  pp.get("gaussian_pairs_radius",      rad_pairs);
  pp.get("init_energy",                init_energy);
  pp.getarr("gaussian_pairs_center",   vec, 0, SpaceDim); center_pairs   = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  
  List<Particle> electron_ion_pairs;

  // Add various types of particles
  add_uniform_particles(electron_ion_pairs,  round(uniform_pairs),    weight);
  add_gaussian_particles(electron_ion_pairs, round(gaussian_pairs),   weight, rad_pairs,   center_pairs);

  // Give electrons some energies
  List<Particle> electron_energies = electron_ion_pairs;
  for (ListIterator<Particle> li(electron_energies); li.ok(); ){
    Particle& p = li();
    Real& mass = p.mass();
    mass *= init_energy;
  }
  
  // Set initial particles
  m_species[m_eed_idx]->get_initial_particles()  = electron_energies;
  m_species[m_elec_idx]->get_initial_particles() = electron_ion_pairs;
  m_species[m_plus_idx]->get_initial_particles() = electron_ion_pairs;

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
    MayDay::Abort("air_eed::parse_initial_particles - unknown deposition type requested");
  }
  
  for (int i = 0; i < m_num_species; i++){
    m_species[i]->get_deposition() = deposition;
  }
}

void air_eed::parse_chemistry(){
  ParmParse pp("air_eed");

  std::string str;
  pp.get("chemistry_dt", m_chemistry_dt);
  pp.get("chemistry_algorithm", str);

  if(str == "euler"){
    m_chemistryAlgorithm = chemistryAlgorithm::euler;
  }
  else if(str == "rk2"){
    m_chemistryAlgorithm = chemistryAlgorithm::rk2;
  }
  else if(str == "rk4"){
    m_chemistryAlgorithm = chemistryAlgorithm::rk4;
  }
  else{
    MayDay::Abort("air_eed::parse_chemistry - unknown chemistry algorithm requested");
  }
}

void air_eed::add_uniform_particles(List<Particle>& a_particles, const int a_num, const Real a_weight){

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

void air_eed::add_gaussian_particles(List<Particle>& a_particles,
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

RealVect air_eed::random_gaussian(){

  const Real rad = m_gauss(*m_rng);
  return rad*random_direction();
}

RealVect air_eed::random_direction(){
#if CH_SPACEDIM == 2
  return random_direction2D();
#else
  return random_direction3D();
#endif
}

#if CH_SPACEDIM == 2
RealVect air_eed::random_direction2D(){
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
RealVect air_eed::random_direction3D(){
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

void air_eed::parse_domain_bc(){

  ParmParse pp("air_eed");
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

void air_eed::advance_reaction_network(Vector<Real>&          a_particle_sources,
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

  Vector<Real>     cdr_src(m_num_species, 0.0);
  Vector<Real>     rte_src(m_num_photons, 0.0);
  Vector<Real>     cdr_phi(m_num_species, 0.0);
  Vector<Real>     rte_phi(m_num_photons, 0.0);
  Vector<RealVect> cdr_grad(m_num_species, RealVect::Zero);


  // Copy starting data
  for (int i = 0; i < m_num_species; i++){
    cdr_phi[i]  = a_particle_densities[i];
    cdr_grad[i] = a_particle_gradients[i];
  }
  
  rte_phi[m_photon_idx]          = a_photon_densities[m_photon_idx];
  a_photon_sources[m_photon_idx] = 0.0;



  int nsteps = ceil(a_dt/m_chemistry_dt);
  Real dt    = a_dt/nsteps;
  Real time  = a_time;
  
  for (int istep = 0; istep < nsteps; istep++){
    if(m_chemistryAlgorithm == chemistryAlgorithm::euler){
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Increment
      for (int i = 0; i < m_num_species; i++){
	cdr_phi[i] += cdr_src[i]*dt;
	cdr_phi[i] = Max(0.0, cdr_phi[i]);
      }

      // Add photons produced in the substep
      a_photon_sources[m_photon_idx] += rte_src[m_photon_idx];
    }
    else if(m_chemistryAlgorithm == chemistryAlgorithm::rk2){

      // Compute slope at k
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      Vector<Real> k1 = cdr_src;

      // Euler update to k+1
      for (int i = 0; i < m_num_species; i++){
	cdr_phi[i] += cdr_src[i]*dt;
	cdr_phi[i] = Max(0.0, cdr_phi[i]);
      }

      // Photons only use the Euler update
      a_photon_sources[m_photon_idx] += rte_src[m_photon_idx];

      // Re-compute slope at k
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Funky notation, but checks out
      for (int i = 0; i < m_num_species; i++){
	cdr_phi[i] = cdr_phi[i] + dt*(0.5*cdr_src[i] - 0.5*k1[i]);
      }
    }
    else if(m_chemistryAlgorithm == chemistryAlgorithm::rk4){
      const Vector<Real> cdr_phi0 = cdr_phi;
      
      // Compute k1 slope
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k1 = cdr_src;

      // Only Euler update for photons. 
      a_photon_sources[m_photon_idx] += rte_src[m_photon_idx];

      // Compute k2 slope
      for (int i = 0; i < m_num_species; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k1[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k2 = cdr_src;

      // Compute k3 slope
      for (int i = 0; i < m_num_species; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k2[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k3 = cdr_src;

      // Compute k4 slope
      for (int i = 0; i < m_num_species; i++){
	cdr_phi[i] = cdr_phi0[i] + dt*k3[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k4 = cdr_src;

      for (int i = 0; i < m_num_species; i++){
	cdr_phi[i] = cdr_phi0[i] + dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.;
      }
    }
    else{
      MayDay::Abort("advance_reaction_network - not supporting other chemistry algorithm yet");
    }

    // Time at next substep
    time += dt;
  }

  // Linearize source term
  for (int i = 0; i < m_num_species; i++){
    a_particle_sources[i] = (cdr_phi[i] - a_particle_densities[i])/a_dt;
  }
}

void air_eed::advance_chemistry_euler(Vector<Real>&          a_particle_sources,
				      Vector<Real>&          a_photon_sources,
				      Vector<Real>&          a_particle_densities,
				      const Vector<RealVect> a_particle_gradients,
				      const Vector<Real>     a_photon_densities,
				      const RealVect         a_E,
				      const RealVect         a_pos,
				      const Real             a_dx,
				      const Real             a_dt,
				      const Real             a_time,
				      const Real             a_kappa) const{
  // R0: E.J -> heating
  // R1: e + M -> e + e + M+  alpha*Xe. Loss is m_ionization_loss, which is a bulk avlue
  // R2: e + M -> M-+         eta*Xe. Loss is equal to the mean energy.
  // R3: e + M -> e + M       "Friction" loss. Tabulated with bulk value
  // R4: Y + M -> e + M+      Energy gain equal to input value (1eV is fine)
  // R5: e + M -> e + M + Y   Excitation losses are a part of R3
  const Real volume = pow(a_dx, SpaceDim);
  const Real energy = compute_electron_energy(a_particle_densities[m_eed_idx], a_particle_densities[m_elec_idx]);

  const Real De     = m_e_diffco.get_entry(energy);
  const RealVect Ve = -m_e_mobility.get_entry(energy)*a_E;
  const RealVect Je = Ve*a_particle_densities[m_elec_idx] - De*a_particle_gradients[m_elec_idx];
  const Real ve     = Ve.vectorLength();
  
  // Ionization and attachment coefficients
  const Real alpha  = m_e_alpha.get_entry(energy);
  const Real eta    = m_e_eta.get_entry(energy);

  const Real R0 = -PolyGeom::dot(Je, a_E);                    // This is a energy/volume rate. Units of ev/(m^3s)
  const Real R1 = alpha*ve*a_particle_densities[m_elec_idx];  // This is a volume rate, units of 1/(m^3s)
  const Real R2 = eta*ve*a_particle_densities[m_elec_idx];    // This is a volume rate, units of 1/(m^3s)
  const Real R3 = m_e_losses.get_entry(energy)*a_particle_densities[m_elec_idx]; //Loss rate. Units of eV/(m^3s)
  const Real R4 = a_photon_densities[m_photon_idx]/a_dt;

  Real& SE = a_particle_sources[m_eed_idx];
  Real& Se = a_particle_sources[m_elec_idx];
  Real& Sp = a_particle_sources[m_plus_idx];
  Real& Sm = a_particle_sources[m_minu_idx];

  // Heating. 
  SE = 0.0;
  Se = 0.0;
  Sp = 0.0;
  Sm = 0.0;

  // R0: E.J -> Heating
  SE += R0;

  // R1: e + M => e + e + M+
  //SE -= R1*m_ionization_loss;
  Se += R1;
  Sp += R1;

  // R2: e + M => M-
  //  SE -= R2*energy;
  Se -= R2;
  Sm += R2;

  // R3: e + M => e + M + losses
  SE -= R3;

  // R4: Photoionization, M + y => e + M+
  SE += m_photoi_gain*R4;
  Se += R4;
  Sp += R4;

  
  // R5: Propensity functions for photon emission
  const Real quench         = m_pq/(m_pq+m_p);
  const Real prop_c4v0_X1v0 = quench*m_c4v0_X1v0_photoi_eff*m_c4v0_exc_eff*R1*volume;
  const Real prop_c4v0_X1v1 = quench*m_c4v0_X1v1_photoi_eff*m_c4v0_exc_eff*R1*volume;
  const Real prop_c4v1_X1v0 = quench*m_c4v1_X1v0_photoi_eff*m_c4v1_exc_eff*R1*volume;
  const Real prop_c4v1_X1v1 = quench*m_c4v1_X1v1_photoi_eff*m_c4v1_exc_eff*R1*volume;
  const Real prop_c4v1_X1v2 = quench*m_c4v1_X1v2_photoi_eff*m_c4v1_exc_eff*R1*volume;
  const Real prop_c4v1_X1v3 = quench*m_c4v1_X1v3_photoi_eff*m_c4v1_exc_eff*R1*volume;
  const Real prop_b1v1_X1v0 = quench*m_b1v1_X1v0_photoi_eff*m_b1v1_exc_eff*R1*volume;
  const Real prop_b1v1_X1v1 = quench*m_b1v1_X1v1_photoi_eff*m_b1v1_exc_eff*R1*volume;

  // R5: Draw total number of photons. The photon type is late-resolved. 
  int num_photons = 0;
  num_photons += poisson_reaction(prop_c4v0_X1v0, a_dt);
  num_photons += poisson_reaction(prop_c4v0_X1v1, a_dt);
  num_photons += poisson_reaction(prop_c4v1_X1v0, a_dt);
  num_photons += poisson_reaction(prop_c4v1_X1v1, a_dt);
  num_photons += poisson_reaction(prop_c4v1_X1v2, a_dt);
  num_photons += poisson_reaction(prop_c4v1_X1v3, a_dt);
  num_photons += poisson_reaction(prop_b1v1_X1v0, a_dt);
  num_photons += poisson_reaction(prop_b1v1_X1v1, a_dt);

  a_photon_sources[m_photon_idx] = 1.0*num_photons;

  return;
}

int air_eed::poisson_reaction(const Real a_propensity, const Real a_dt) const{
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

Real air_eed::compute_electron_energy(const Real a_electron_energy, const Real a_electron_density) const {
  Real energy = a_electron_energy/(Max(1E3, a_electron_density));
  energy = Min(Max(energy, m_min_energy), m_max_energy);

  return energy;
}

Real air_eed::compute_alpha_eff(const Real a_E) const{
  const Real energy = m_lfa.get_entry(a_E); // Get energy in the LFA approximation
  
  const Real alpha = m_e_alpha.get_entry(energy);
  const Real eta   = m_e_eta.get_entry(energy);

  return (alpha-eta);
}


Vector<Real> air_eed::compute_cdr_diffusion_coefficients(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_E,
							 const Vector<Real> a_cdr_densities) const {

  const Real electron_energy_density = a_cdr_densities[m_eed_idx];
  const Real electron_density        = a_cdr_densities[m_elec_idx];
  const Real electron_energy         = compute_electron_energy(electron_energy_density, electron_density);
  
  Vector<Real> dco(m_num_species, 0.0);
  dco[m_eed_idx]  = m_e_diffco.get_entry(electron_energy)*5.0/3.0;
  dco[m_elec_idx] = m_e_diffco.get_entry(electron_energy);
  dco[m_plus_idx] = m_ion_diffusion;
  dco[m_minu_idx] = m_ion_diffusion;
  
  return dco;
}
  
Vector<RealVect> air_eed::compute_cdr_velocities(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_num_species, RealVect::Zero);
  
  const Real electron_energy_density = a_cdr_densities[m_eed_idx];
  const Real electron_density        = a_cdr_densities[m_elec_idx];
  const Real electron_energy         = compute_electron_energy(electron_energy_density, electron_density);
  Real electron_mobility       = m_e_mobility.get_entry(electron_energy);


  vel[m_eed_idx]  = -(5.0/3.0)*electron_mobility*a_E;
  vel[m_elec_idx] = -electron_mobility*a_E;
  vel[m_plus_idx] =  a_E*m_ion_mobility;
  vel[m_minu_idx] = -a_E*m_ion_mobility;
  
  return vel;
}
  
Vector<Real> air_eed::compute_cdr_domain_fluxes(const Real           a_time,
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
  
Vector<Real> air_eed::compute_cdr_electrode_fluxes(const Real         a_time,
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

Vector<Real> air_eed::compute_cdr_dielectric_fluxes(const Real         a_time,
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

Vector<Real> air_eed::compute_cdr_fluxes(const Real         a_time,
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

  // Drift outflow for now
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = Max(0.0, a_extrap_cdr_fluxes[i]);
  }

  return fluxes;
}

Real air_eed::initial_sigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}


