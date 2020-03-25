/*!
  @file   air7_stephens.H
  @brief  7-species and 8-photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air7_stephens.H"
#include "air7_stephens_species.H"
#include "data_ops.H"
#include "units.H" 

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include <ParmParse.H>
#include <PolyGeom.H>

using namespace physics::cdr_plasma;

std::string air7_stephens::s_bolsig_energy   = "# Mean energy (E/N, eV)";
std::string air7_stephens::s_bolsig_mobility = "# Electron mobility (E/N, mu*N)";
std::string air7_stephens::s_bolsig_diffco   = "# Electron diffusion coefficient (E/N, D*N)";
std::string air7_stephens::s_bolsig_alpha    = "# Townsend alpha (E/N, alpha/N)";
std::string air7_stephens::s_bolsig_eta      = "# Townsend eta (E/N, eta/N)";
std::string air7_stephens::s_bolsig_b1_exc   = "# b1 excitation (E/N, rate/N)";
std::string air7_stephens::s_bolsig_c4_exc   = "# c4 excitation (E/N, rate/N)";
std::string air7_stephens::s_bolsig_alphaN2  = "# N2 Ionization (E/N, rate/N)";
std::string air7_stephens::s_bolsig_alphaO2  = "# O2 Ionization (E/N, rate/N)";

air7_stephens::air7_stephens() {

  instantiate_species();
  parse_transport_file();
  parse_transport();
  parse_chemistry();
  parse_gas_params();
  parse_electron_mobility();
  parse_electron_diffco();
  parse_alpha();
  parse_eta();
  parse_excitations();
  parse_photoi();
  parse_temperature();
  parse_see();
  parse_domain_bc();
  init_rng();                 // Initialize random number generators
}

air7_stephens::~air7_stephens() {

}

void air7_stephens::read_file_entries(lookup_table& a_table, const std::string a_string){
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

void air7_stephens::parse_chemistry(){
  ParmParse pp("air7_stephens");

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

void air7_stephens::parse_transport_file(){
  ParmParse pp("air7_stephens");
  pp.get("transport_file",  m_transport_file);
  pp.get("uniform_tables",  m_uniform_entries);
  std::ifstream infile(m_transport_file);
  if(!infile.good()){
    MayDay::Abort("air7_stephens::parse_transport_file - could not find transport data");
  }
  else{
    infile.close();
  }
}

void air7_stephens::parse_transport(){
  ParmParse pp("air7_stephens");

  std::string str;

  pp.get("use_alpha_corr", str);      m_alpha_corr          = (str == "true") ? true : false;
  pp.get("mobile_electrons", str);    m_mobile_electrons    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_diffusive_electrons = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);      m_diffusive_ions      = (str == "true") ? true : false;
  pp.get("mobile_ions", str);         m_mobile_ions         = (str == "true") ? true : false;
  
  pp.get("ion_mobility", m_ion_mobility);

  m_ion_diffusion = m_ion_mobility*(units::s_kb*m_T)/units::s_Qe;
}

void air7_stephens::parse_gas_params(){
  ParmParse pp("air7_stephens");

  // Pressure form input script
  pp.get("pressure",           m_p);
  pp.get("temperature",        m_T);
  pp.get("frac_N2",            m_N2frac);
  pp.get("frac_O2",            m_O2frac);
  pp.get("photoi_factor",      m_photoi_factor);

  m_p  *= units::s_atm2pascal;
  m_N = m_p*units::s_Na/(m_T*units::s_R);
}

void air7_stephens::parse_electron_mobility(){
  ParmParse pp("air7_stephens");

  read_file_entries(m_e_mobility, air7_stephens::s_bolsig_mobility);
  m_e_mobility.scale_x(m_N*units::s_Td);
  m_e_mobility.scale_y(1./m_N); 
  m_e_mobility.make_uniform(m_uniform_entries);
}

void air7_stephens::parse_electron_diffco(){
  ParmParse pp("air7_stephens");
  
  read_file_entries(m_e_diffco, air7_stephens::s_bolsig_diffco);
  m_e_diffco.scale_x(m_N*units::s_Td);
  m_e_diffco.scale_y(1./m_N); 
  m_e_diffco.make_uniform(m_uniform_entries);
}

void air7_stephens::parse_alpha(){
  ParmParse pp("air7_stephens");
  read_file_entries(m_e_alpha,   air7_stephens::s_bolsig_alpha);
  read_file_entries(m_e_alphaN2, air7_stephens::s_bolsig_alphaN2);
  read_file_entries(m_e_alphaO2, air7_stephens::s_bolsig_alphaO2);
  
  m_e_alpha.scale_x(m_N*units::s_Td);
  m_e_alphaN2.scale_x(m_N*units::s_Td);
  m_e_alphaO2.scale_x(m_N*units::s_Td);
  
  m_e_alpha.scale_y(m_N);
  m_e_alphaN2.scale_y(m_N*m_N2frac);
  m_e_alphaO2.scale_y(m_N*m_O2frac);
  
  m_e_alpha.make_uniform(m_uniform_entries);
  m_e_alphaN2.make_uniform(m_uniform_entries);
  m_e_alphaO2.make_uniform(m_uniform_entries);
}

void air7_stephens::parse_eta(){
  ParmParse pp("air7_stephens");
  read_file_entries(m_e_eta, air7_stephens::s_bolsig_eta);
  m_e_eta.scale_x(m_N*units::s_Td);
  m_e_eta.scale_y(m_N);
  m_e_eta.make_uniform(m_uniform_entries);
}

void air7_stephens::parse_temperature(){
  ParmParse pp("air7_stephens");
  read_file_entries(m_e_temperature, air7_stephens::s_bolsig_energy);

  m_e_temperature.dump_table();
  m_e_temperature.scale_x(m_N*units::s_Td);
  m_e_temperature.scale_y(2.0*units::s_Qe/(3.0*units::s_kb));
  m_e_temperature.make_uniform(m_uniform_entries);
}

void air7_stephens::parse_excitations(){
  ParmParse pp("air7_mc8");
  
  read_file_entries(m_b1_exc, air7_stephens::s_bolsig_b1_exc);
  m_b1_exc.scale_x(m_N*units::s_Td);
  m_b1_exc.scale_y(m_N); 
  m_b1_exc.make_uniform(m_uniform_entries);

  read_file_entries(m_c4_exc, air7_stephens::s_bolsig_c4_exc);

  m_c4_exc.scale_x(m_N*units::s_Td);
  m_c4_exc.scale_y(m_N);

  m_c4_exc.make_uniform(m_uniform_entries);
}

void air7_stephens::parse_see(){
  ParmParse pp("air7_stephens");
  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void air7_stephens::parse_photoi(){
  ParmParse pp("air7_stephens");

  pp.get("c4v0_exc_rep", m_c4v0_exc_rep);
  pp.get("c4v1_exc_rep", m_c4v1_exc_rep);
  pp.get("b1v1_exc_rep", m_b1v1_exc_rep);

  pp.get("c4v0_kp", m_c4v0_kp);
  pp.get("c4v1_kp", m_c4v1_kp);
  pp.get("b1v1_kp", m_b1v1_kp);

  pp.get("c4v0_X1v0_kr", m_c4v0_X1v0_kr); 
  pp.get("c4v0_X1v1_kr", m_c4v0_X1v1_kr);
  pp.get("c4v1_X1v0_kr", m_c4v1_X1v0_kr);
  pp.get("c4v1_X1v1_kr", m_c4v1_X1v1_kr);
  pp.get("c4v1_X1v2_kr", m_c4v1_X1v2_kr);
  pp.get("c4v1_X1v3_kr", m_c4v1_X1v3_kr);
  pp.get("b1v1_X1v0_kr", m_b1v1_X1v0_kr); 
  pp.get("b1v1_X1v1_kr", m_b1v1_X1v1_kr);

  pp.get("c4v0_X1v0_photoi_eff", m_c4v0_X1v0_photoi_eff);
  pp.get("c4v0_X1v1_photoi_eff", m_c4v0_X1v1_photoi_eff);
  pp.get("c4v1_X1v0_photoi_eff", m_c4v1_X1v0_photoi_eff);
  pp.get("c4v1_X1v1_photoi_eff", m_c4v1_X1v1_photoi_eff);
  pp.get("c4v1_X1v2_photoi_eff", m_c4v1_X1v2_photoi_eff);
  pp.get("c4v1_X1v3_photoi_eff", m_c4v1_X1v3_photoi_eff);
  pp.get("b1v1_X1v0_photoi_eff", m_b1v1_X1v0_photoi_eff);
  pp.get("b1v1_X1v1_photoi_eff", m_b1v1_X1v1_photoi_eff);

  pp.get("k_quench",  m_kq);
  pp.get("photoi_factor", m_photoi_factor);

  m_kq *= m_N;

  m_c4v0_kr = m_c4v0_X1v0_kr + m_c4v0_X1v1_kr;
  m_c4v0_k  = m_c4v0_kr + m_c4v0_kp + m_kq;

  m_c4v1_kr = m_c4v1_X1v0_kr + m_c4v1_X1v1_kr + m_c4v1_X1v2_kr + m_c4v1_X1v3_kr;
  m_c4v1_k  = m_c4v1_kr + m_c4v1_kp + m_kq;

  m_b1v1_kr = m_b1v1_X1v0_kr + m_b1v1_X1v1_kr;
  m_b1v1_k  = m_b1v1_kr + m_b1v1_kp + m_kq;
}



void air7_stephens::init_rng(){
  ParmParse pp("air7_stephens");
  pp.get("rng_seed", m_rng_seed);

  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  if(m_rng_seed < 0) {
    m_rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  m_udist11 = new std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_rng     = new std::mt19937_64(m_rng_seed);
}

void air7_stephens::instantiate_species(){
  m_num_cdr_species = 7;
  m_num_rte_species = 8;

  m_elec_idx     = 0;
  m_N2plus_idx   = 1;
  m_O2plus_idx   = 2;
  m_N4plus_idx   = 3;
  m_O4plus_idx   = 4;
  m_O2plusN2_idx = 5;
  m_O2minus_idx  = 6;

  m_c4v0_X1v0_idx = 0;
  m_c4v0_X1v1_idx = 1;
  m_c4v1_X1v0_idx = 2;
  m_c4v1_X1v1_idx = 3;
  m_c4v1_X1v2_idx = 4;
  m_c4v1_X1v3_idx = 5;
  m_b1v1_X1v0_idx = 6;
  m_b1v1_X1v1_idx = 7;

  m_cdr_species.resize(m_num_cdr_species);
  m_cdr_species[m_elec_idx]      = RefCountedPtr<cdr_species>      (new air7_stephens::electron());
  m_cdr_species[m_N2plus_idx]    = RefCountedPtr<cdr_species>      (new air7_stephens::N2plus());
  m_cdr_species[m_O2plus_idx]    = RefCountedPtr<cdr_species>      (new air7_stephens::O2plus());
  m_cdr_species[m_N4plus_idx]    = RefCountedPtr<cdr_species>      (new air7_stephens::N4plus());
  m_cdr_species[m_O4plus_idx]    = RefCountedPtr<cdr_species>      (new air7_stephens::O4plus());
  m_cdr_species[m_O2plusN2_idx]  = RefCountedPtr<cdr_species>      (new air7_stephens::O2plusN2());
  m_cdr_species[m_O2minus_idx]   = RefCountedPtr<cdr_species>      (new air7_stephens::O2minus());

  m_rte_species.resize(m_num_rte_species);
  m_rte_species[m_c4v0_X1v0_idx] = RefCountedPtr<rte_species> (new air7_stephens::phot_c4v0_X1v0());
  m_rte_species[m_c4v0_X1v1_idx] = RefCountedPtr<rte_species> (new air7_stephens::phot_c4v0_X1v1());
  m_rte_species[m_c4v1_X1v0_idx] = RefCountedPtr<rte_species> (new air7_stephens::phot_c4v1_X1v0());
  m_rte_species[m_c4v1_X1v1_idx] = RefCountedPtr<rte_species> (new air7_stephens::phot_c4v1_X1v1());
  m_rte_species[m_c4v1_X1v2_idx] = RefCountedPtr<rte_species> (new air7_stephens::phot_c4v1_X1v2());
  m_rte_species[m_c4v1_X1v3_idx] = RefCountedPtr<rte_species> (new air7_stephens::phot_c4v1_X1v3());
  m_rte_species[m_b1v1_X1v0_idx] = RefCountedPtr<rte_species> (new air7_stephens::phot_b1v1_X1v0());
  m_rte_species[m_b1v1_X1v1_idx] = RefCountedPtr<rte_species> (new air7_stephens::phot_b1v1_X1v1());

}

void air7_stephens::parse_domain_bc(){

  ParmParse pp("air7_stephens");
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

void air7_stephens::advance_reaction_network(Vector<Real>&          a_particle_sources,
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
  Vector<Real>     cdr_src(m_num_cdr_species, 0.0);
  Vector<Real>     rte_src(m_num_rte_species, 0.0);
  Vector<Real>     cdr_phi(m_num_cdr_species, 0.0);
  Vector<Real>     rte_phi(m_num_rte_species, 0.0);
  Vector<RealVect> cdr_grad(m_num_cdr_species, RealVect::Zero);


  // Copy starting data
  for (int i = 0; i < m_num_cdr_species; i++){
    cdr_phi[i]  = a_particle_densities[i];
    cdr_grad[i] = a_particle_gradients[i];
  }

  for (int i = 0; i < m_num_rte_species; i++){
    rte_phi[i]          = a_photon_densities[i];
    a_photon_sources[i] = 0.0;
  }

  int nsteps = ceil(a_dt/m_chemistry_dt);
  Real dt    = a_dt/nsteps;
  Real time  = a_time;
  
  for (int istep = 0; istep < nsteps; istep++){
    if(m_chemistryAlgorithm == chemistryAlgorithm::euler){
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Increment
      for (int i = 0; i < m_num_cdr_species; i++){
	cdr_phi[i] = cdr_phi[i] + cdr_src[i]*dt;
      }

      // Add photons produced in the substep
      for (int i = 0; i < m_num_rte_species; i++){
	a_photon_sources[i] += rte_src[i];
      }
    }
    else if(m_chemistryAlgorithm == chemistryAlgorithm::rk2){

      // Compute slope at k
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      Vector<Real> k1 = cdr_src;

      // Euler update to k+1
      for (int i = 0; i < m_num_cdr_species; i++){
	cdr_phi[i] += cdr_src[i]*dt;
	cdr_phi[i] = Max(0.0, cdr_phi[i]);
      }

      // Photons only use the Euler update
      for (int i = 0; i < m_num_rte_species; i++){
	a_photon_sources[i] += rte_src[i];
      }

      // Re-compute slope at k
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Funky notation, but checks out
      for (int i = 0; i < m_num_cdr_species; i++){
	cdr_phi[i] = cdr_phi[i] + dt*(0.5*cdr_src[i] - 0.5*k1[i]);
      }
    }
    else if(m_chemistryAlgorithm == chemistryAlgorithm::rk4){
      const Vector<Real> cdr_phi0 = cdr_phi;
      
      // Compute k1 slope
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k1 = cdr_src;

      // Only Euler update for photons.
      for (int i = 0; i < m_num_rte_species; i++){
	a_photon_sources[i] += rte_src[i];
      }

      // Compute k2 slope
      for (int i = 0; i < m_num_cdr_species; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k1[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k2 = cdr_src;

      // Compute k3 slope
      for (int i = 0; i < m_num_cdr_species; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k2[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k3 = cdr_src;

      // Compute k4 slope
      for (int i = 0; i < m_num_cdr_species; i++){
	cdr_phi[i] = cdr_phi0[i] + dt*k3[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k4 = cdr_src;

      for (int i = 0; i < m_num_cdr_species; i++){
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
  for (int i = 0; i < m_num_cdr_species; i++){
    a_particle_sources[i] = (cdr_phi[i] - a_particle_densities[i])/a_dt;
  }

}

void air7_stephens::advance_chemistry_euler(Vector<Real>&          a_particle_sources,
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
  const Real volume = pow(a_dx, SpaceDim);
  const Real E      = a_E.vectorLength();
  const Real ve     = E*m_e_mobility.get_entry(E);
  const Real Te     = Max(m_e_temperature.get_entry(E), 300.);

  const Real Ne    = a_particle_densities[m_elec_idx];
  const Real N2p   = a_particle_densities[m_N2plus_idx];
  const Real O2p   = a_particle_densities[m_O2plus_idx];
  const Real N4p   = a_particle_densities[m_N4plus_idx];
  const Real O4p   = a_particle_densities[m_O4plus_idx];
  const Real O2pN2 = a_particle_densities[m_O2plusN2_idx];
  const Real O2m   = a_particle_densities[m_O2minus_idx];

  const Real M  = m_N;
  const Real N2 = m_N*m_N2frac;
  const Real O2 = m_N*m_O2frac;

  Real& S_e     = a_particle_sources[m_elec_idx];
  Real& S_N2p   = a_particle_sources[m_N2plus_idx];
  Real& S_O2p   = a_particle_sources[m_O2plus_idx];
  Real& S_N4p   = a_particle_sources[m_N4plus_idx];
  Real& S_O4p   = a_particle_sources[m_O4plus_idx];
  Real& S_O2pN2 = a_particle_sources[m_O2plusN2_idx];
  Real& S_O2m   = a_particle_sources[m_O2minus_idx];

  S_e     = 0.0;
  S_N2p   = 0.0;
  S_O2p   = 0.0;
  S_N4p   = 0.0;
  S_O4p   = 0.0;
  S_O2pN2 = 0.0;
  S_O2m   = 0.0;

  Real fcorr = 1.0;
  if(m_alpha_corr){
    const RealVect Eunit = a_E/a_E.vectorLength();
    const Real De        = m_e_diffco.get_entry(E);
    const RealVect gNe   = a_particle_gradients[m_elec_idx];

    fcorr = 1.0 + PolyGeom::dot(Eunit, De*gNe)/(1.0+a_particle_densities[m_elec_idx]*ve);
    fcorr = Min(fcorr, 1.0);
    fcorr = Max(0.0, fcorr);
  }

  const Real R1  = fcorr*m_e_alphaN2.get_entry(E)*Ne;
  const Real R2  = fcorr*m_e_alphaO2.get_entry(E)*Ne;
  const Real R3  = (5.E-41)*N2p*N2*M;
  const Real R4  = 2.5E-16*N4p*O2;
  const Real R5  = 6E-17*N2p*O2;
  const Real R6  = 9E-43*O2p*N2*N2;
  const Real R7  = 4.3E-16*O2pN2*N2;
  const Real R8  = 1E-15*O2pN2*O2;
  const Real R9  = 2.4E-42*O2p*O2*M;
  const Real R10 = 1.4E-12*sqrt(300./Te)*Ne*O4p;
  const Real R11 = 2E-13*(300./Te)*Ne*O2p;
  const Real R12 = 2E-41*(300./Te)*Ne*O2*O2;
  const Real R13 = 1E-13*O2m*O4p;
  const Real R14 = 2E-37*O2m*O4p*M;
  const Real R15 = 2E-37*O2m*O2p*M;

  // R1: e + N2 -> e + e + N2plus
  S_e   += R1;
  S_N2p += R1;

  // R2: e + O2 -> e + e + O2plus
  S_e   += R2;
  S_O2p += R2;

  // R3:  N2plus + N2 + M -> N4plus + M
  S_N2p -= R3;
  S_N4p += R3;

  // R4:  N4plus + O2 -> O2plus + N2 + N2
  S_N4p -= R4;
  S_O2p += R4;

  // R5:  N2plus + O2 -> O2plus + N2
  S_N2p -= R5;
  S_O2p += R5;

  // R6:  O2plus + N2 + N -> O2plusN2 + N2
  S_O2p   -= R6;
  S_O2pN2 += R6;

  // R7:  O2plusN2 + N2 -> O2plus + N2 + N2
  S_O2pN2 -= R7;
  S_O2p   += R7;

  // R8:  O2plusN2 + O2 -> O4plus + N2
  S_O2pN2 -= R8;
  S_O4p   += R8;

  // R9:  O2plus + O2 + M -> O4plus + M
  S_O2p -= R9;
  S_O4p += R9;

  // R10: e + O4plus -> O2 + O2
  S_e   -= R10;
  S_O4p -= R10;

  // R11: e + O2plus -> O + O
  S_e   -= R11;
  S_O2p -= R11;

  // R12: e + O2 + O2 -> O2minus + O2
  S_e   -= R12;
  S_O2m += R12;

  // R13: O2minus + O4plus -> O2 + O2 + O2
  S_O2m -= R13;
  S_O4p -= R13;

  // R14: O2minus + O4plus + M -> O2 + O2 + O2 + M
  S_O2m -= R14;
  S_O4p -= R14;

  // R15: O2minus + O2plus + M -> O2 + O2 + M
  S_O2m -= R15;
  S_O2p -= R15;

  // Photoionization
  for (int i = 0; i < m_num_rte_species; i++){
    S_e   += a_photon_densities[i]/a_dt;
    S_O2p += a_photon_densities[i]/a_dt;
  }

  // Photoionization code below
  const Real Rb1    = m_b1_exc.get_entry(E)*a_particle_densities[m_elec_idx];  // Excitation rate for b1
  const Real Rc4    = m_c4_exc.get_entry(E)*a_particle_densities[m_elec_idx];  // Excitation rate for c4
  
  const int num_exc_c4v0 = poisson_reaction(m_photoi_factor*Rc4*volume*m_c4v0_exc_rep, a_dt);
  const int num_exc_c4v1 = poisson_reaction(m_photoi_factor*Rc4*volume*m_c4v1_exc_rep, a_dt);
  const int num_exc_b1v1 = poisson_reaction(m_photoi_factor*Rb1*volume*m_b1v1_exc_rep, a_dt);
  
  // Determine number of radiative de-excitations
  const int num_c4v0_rad = binomial_trials(num_exc_c4v0, m_c4v0_kr/m_c4v0_k);
  const int num_c4v1_rad = binomial_trials(num_exc_c4v1, m_c4v1_kr/m_c4v1_k);
  const int num_b1v1_rad = binomial_trials(num_exc_b1v1, m_b1v1_kr/m_b1v1_k);

  // For c4v0, determine distribution of radiative de-excitations
  int num_c4v0_X1v0 = binomial_trials(num_c4v0_rad, m_c4v0_X1v0_kr/m_c4v0_kr);
  int num_c4v0_X1v1 = num_c4v0_rad - num_c4v0_X1v0; // Rest must be other transition

  // For b1v1, determine distribution of radiative de-excitations
  int num_b1v1_X1v0 = binomial_trials(num_b1v1_rad, m_b1v1_X1v0_kr/m_b1v1_kr);
  int num_b1v1_X1v1 = num_b1v1_rad - num_b1v1_X1v0; // Rest must be other transition

  // 4 transitions for c4v1. C++ doesn't have a multinomial implementation, so do some nested binomial magic instead
  int num_c4v1_X1v0 = binomial_trials(num_c4v1_rad, m_c4v1_X1v0_kr/m_c4v1_kr);
  int num_c4v1_X1v1 = binomial_trials(num_c4v1_rad - num_c4v1_X1v0,
				      m_c4v1_X1v1_kr/(m_c4v1_X1v1_kr + m_c4v1_X1v2_kr, m_c4v1_X1v3_kr));
  int num_c4v1_X1v2 = binomial_trials(num_c4v1_rad - (num_c4v1_X1v0 + num_c4v1_X1v1),
				      m_c4v1_X1v2_kr/(m_c4v1_X1v2_kr, m_c4v1_X1v3_kr));
  int num_c4v1_X1v3 = num_c4v1_rad - (num_c4v1_X1v0 + num_c4v1_X1v1 + num_c4v1_X1v2);


  // Check if the transitions lead to photoionization
  num_c4v0_X1v0 = binomial_trials(num_c4v0_X1v0, m_c4v0_X1v0_photoi_eff);
  num_c4v0_X1v1 = binomial_trials(num_c4v0_X1v1, m_c4v0_X1v1_photoi_eff);
  num_c4v1_X1v0 = binomial_trials(num_c4v1_X1v0, m_c4v1_X1v0_photoi_eff);
  num_c4v1_X1v1 = binomial_trials(num_c4v1_X1v1, m_c4v1_X1v1_photoi_eff);
  num_c4v1_X1v2 = binomial_trials(num_c4v1_X1v2, m_c4v1_X1v2_photoi_eff);
  num_c4v1_X1v3 = binomial_trials(num_c4v1_X1v3, m_c4v1_X1v3_photoi_eff);
  num_b1v1_X1v0 = binomial_trials(num_b1v1_X1v0, m_b1v1_X1v0_photoi_eff);
  num_b1v1_X1v1 = binomial_trials(num_b1v1_X1v1, m_b1v1_X1v1_photoi_eff);

  a_photon_sources[m_c4v0_X1v0_idx] = 1.0*num_c4v0_X1v0;
  a_photon_sources[m_c4v0_X1v1_idx] = 1.0*num_c4v0_X1v1;
  a_photon_sources[m_c4v1_X1v0_idx] = 1.0*num_c4v1_X1v0;
  a_photon_sources[m_c4v1_X1v1_idx] = 1.0*num_c4v1_X1v1;
  a_photon_sources[m_c4v1_X1v2_idx] = 1.0*num_c4v1_X1v2;
  a_photon_sources[m_c4v1_X1v3_idx] = 1.0*num_c4v1_X1v3;
  a_photon_sources[m_b1v1_X1v0_idx] = 1.0*num_b1v1_X1v0;
  a_photon_sources[m_b1v1_X1v1_idx] = 1.0*num_b1v1_X1v1;

  return;
}

int air7_stephens::poisson_reaction(const Real a_propensity, const Real a_dt) const{
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

int air7_stephens::binomial_trials(const int a_trials, const Real a_p) const{
  std::binomial_distribution<int> dist(a_trials, a_p);
  return dist(*m_rng);
}


Vector<Real> air7_stephens::compute_cdr_diffusion_coefficients(const Real         a_time,
							       const RealVect     a_pos,
							       const RealVect     a_E,
							       const Vector<Real> a_cdr_densities) const {

  Vector<Real> dco(m_num_cdr_species, 0.0);
  dco[m_elec_idx] = m_e_diffco.get_entry(a_E.vectorLength());
  if(m_diffusive_ions){
    dco[m_N2plus_idx]   = m_ion_diffusion;
    dco[m_O2plus_idx]   = m_ion_diffusion;
    dco[m_N4plus_idx]   = m_ion_diffusion;
    dco[m_O4plus_idx]   = m_ion_diffusion;
    dco[m_O2plusN2_idx] = m_ion_diffusion;
    dco[m_O2minus_idx]  = m_ion_diffusion;
  }

  return dco;

}
  
Vector<RealVect> air7_stephens::compute_cdr_velocities(const Real         a_time,
						       const RealVect     a_pos,
						       const RealVect     a_E,
						       const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_num_cdr_species, RealVect::Zero);

  vel[m_elec_idx] = -a_E*m_e_mobility.get_entry(a_E.vectorLength());
  if(m_mobile_ions){
    vel[m_N2plus_idx]   =  a_E*m_ion_mobility;
    vel[m_O2plus_idx]   =  a_E*m_ion_mobility;
    vel[m_N4plus_idx]   =  a_E*m_ion_mobility;
    vel[m_O4plus_idx]   =  a_E*m_ion_mobility;
    vel[m_O2plusN2_idx] =  a_E*m_ion_mobility;
    vel[m_O2minus_idx]  = -a_E*m_ion_mobility;
  }
  
  return vel;
}
  
Vector<Real> air7_stephens::compute_cdr_domain_fluxes(const Real           a_time,
						      const RealVect       a_pos,
						      const int            a_dir,
						      const Side::LoHiSide a_side,
						      const RealVect       a_E,
						      const Vector<Real>   a_cdr_densities,
						      const Vector<Real>   a_cdr_velocities,
						      const Vector<Real>   a_cdr_gradients,
						      const Vector<Real>   a_rte_fluxes,
						      const Vector<Real>   a_extrap_cdr_fluxes) const{
  Vector<Real> fluxes(m_num_cdr_species, 0.0);

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
  
Vector<Real> air7_stephens::compute_cdr_electrode_fluxes(const Real         a_time,
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

Vector<Real> air7_stephens::compute_cdr_dielectric_fluxes(const Real         a_time,
							  const RealVect     a_pos,
							  const RealVect     a_normal,
							  const RealVect     a_E,
							  const Vector<Real> a_cdr_densities,
							  const Vector<Real> a_cdr_velocities,
							  const Vector<Real> a_cdr_gradients,
							  const Vector<Real> a_rte_fluxes,
							  const Vector<Real> a_extrap_cdr_fluxes) const{
  return Vector<Real>(m_num_cdr_species, 0.0);
  return compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			    a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

Vector<Real> air7_stephens::compute_cdr_fluxes(const Real         a_time,
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
  Vector<Real> fluxes(m_num_cdr_species, 0.0);

  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;

  // Switch for setting drift flux to zero for charge species
  Vector<Real> aj(m_num_cdr_species, 0.0);
  for (int i = 0; i < m_num_cdr_species; i++){
    if(data_ops::sgn(m_cdr_species[i]->get_charge())*PolyGeom::dot(a_E, a_normal) < 0){
      aj[i] = 1.0;
    }
    else {
      aj[i] = 0.0;
    }
  }

  // Drift outflow for now
  for (int i = 0; i < m_num_cdr_species; i++){
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  return fluxes;
}

Real air7_stephens::initial_sigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

Real air7_stephens::compute_alpha(const RealVect a_E) const{
  const Real E     = a_E.vectorLength();
  const Real alpha = m_e_alpha.get_entry(E);
  const Real eta   = m_e_eta.get_entry(E);

  return alpha;
}
