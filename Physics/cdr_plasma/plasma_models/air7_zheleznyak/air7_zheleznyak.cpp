/*!
  @file   air7_zheleznyak.H
  @brief  3-species and 8-Photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air7_zheleznyak.H"
#include "air7_zheleznyak_species.H"
#include "data_ops.H"
#include "units.H" 

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include <ParmParse.H>
#include <PolyGeom.H>

#include "CD_NamespaceHeader.H"
using namespace physics::cdr_plasma;

std::string air7_zheleznyak::s_bolsig_energy   = "# Electron mean energy (E/N, eV)";
std::string air7_zheleznyak::s_bolsig_mobility = "# Electron mobility (E/N, mu*N)";
std::string air7_zheleznyak::s_bolsig_diffco   = "# Electron diffusion coefficient (E/N, D*N)";
std::string air7_zheleznyak::s_bolsig_alpha    = "# Townsend ionization coeff (E/N, alpha/N)";
std::string air7_zheleznyak::s_bolsig_eta      = "# Townsend attachment coeff (E/N, eta/N)";
std::string air7_zheleznyak::s_bolsig_alphaN2  = "# N2 ionization (E/N, rate/N)";
std::string air7_zheleznyak::s_bolsig_alphaO2  = "# O2 ionization (E/N, rate/N)";

air7_zheleznyak::air7_zheleznyak() {

  instantiate_species();      
  parse_transport_file();
  parse_transport();
  parse_chemistry();
  parse_gas_params();
  parse_electron_mobility();
  parse_electron_diffco();
  parse_alpha();
  parse_eta();
  parse_temperature();
  parse_see();
  parseDomainBc();
  init_rng();                 // Initialize random number generators
}

air7_zheleznyak::~air7_zheleznyak() {

}

void air7_zheleznyak::read_file_entries(lookup_table& a_table, const std::string a_string){
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

void air7_zheleznyak::parse_chemistry(){
  ParmParse pp("air7_zheleznyak");

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

void air7_zheleznyak::parse_transport_file(){
  ParmParse pp("air7_zheleznyak");
  pp.get("transport_file",  m_transport_file);
  pp.get("uniform_tables",  m_uniform_entries);
  std::ifstream infile(m_transport_file);
  if(!infile.good()){
    MayDay::Abort("air7_zheleznyak::parse_transport_file - could not find transport data");
  }
  else{
    infile.close();
  }
}

void air7_zheleznyak::parse_transport(){
  ParmParse pp("air7_zheleznyak");

  std::string str;

  pp.get("use_alpha_corr", str);      m_alpha_corr          = (str == "true") ? true : false;
  pp.get("mobile_electrons", str);    m_isMobile_electrons    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_isDiffusive_electrons = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);      m_isDiffusive_ions      = (str == "true") ? true : false;
  pp.get("mobile_ions", str);         m_isMobile_ions         = (str == "true") ? true : false;
  
  pp.get("ion_mobility", m_ion_mobility);

  m_ion_diffusion = m_ion_mobility*(units::s_kb*m_T)/units::s_Qe;
}

void air7_zheleznyak::parse_gas_params(){
  ParmParse pp("air7_zheleznyak");

  // Pressure form input script
  pp.get("pressure",           m_p);
  pp.get("quenching_pressure", m_pq);
  pp.get("temperature",        m_T);
  pp.get("frac_N2",            m_N2frac);
  pp.get("frac_O2",            m_O2frac);
  pp.get("photoi_factor",      m_photoi_factor);

  m_p  *= units::s_atm2pascal;
  m_pq *= units::s_atm2pascal;
  m_N = m_p*units::s_Na/(m_T*units::s_R);
}

void air7_zheleznyak::parse_electron_mobility(){
  ParmParse pp("air7_zheleznyak");

  read_file_entries(m_e_mobility, air7_zheleznyak::s_bolsig_mobility);
  m_e_mobility.scale_x(m_N*units::s_Td);
  m_e_mobility.scale_y(1./m_N); 
  m_e_mobility.make_uniform(m_uniform_entries);
}

void air7_zheleznyak::parse_electron_diffco(){
  ParmParse pp("air7_zheleznyak");
  
  read_file_entries(m_e_diffco, air7_zheleznyak::s_bolsig_diffco);
  m_e_diffco.scale_x(m_N*units::s_Td);
  m_e_diffco.scale_y(1./m_N); 
  m_e_diffco.make_uniform(m_uniform_entries);
}

void air7_zheleznyak::parse_alpha(){
  ParmParse pp("air7_zheleznyak");
  read_file_entries(m_e_alpha,   air7_zheleznyak::s_bolsig_alpha);
  read_file_entries(m_e_alphaN2, air7_zheleznyak::s_bolsig_alphaN2);
  read_file_entries(m_e_alphaO2, air7_zheleznyak::s_bolsig_alphaO2);
  
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

void air7_zheleznyak::parse_eta(){
  ParmParse pp("air7_zheleznyak");
  read_file_entries(m_e_eta, air7_zheleznyak::s_bolsig_eta);
  m_e_eta.scale_x(m_N*units::s_Td);
  m_e_eta.scale_y(m_N);
  m_e_eta.make_uniform(m_uniform_entries);
}

void air7_zheleznyak::parse_temperature(){
  ParmParse pp("air7_zheleznyak");
  read_file_entries(m_e_temperature, air7_zheleznyak::s_bolsig_energy);
  m_e_temperature.scale_x(m_N*units::s_Td);
  m_e_temperature.scale_y(2.0*units::s_Qe/(3.0*units::s_kb));
  m_e_temperature.make_uniform(m_uniform_entries);
}

void air7_zheleznyak::parse_see(){
  ParmParse pp("air7_zheleznyak");
  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void air7_zheleznyak::init_rng(){
  ParmParse pp("air7_zheleznyak");
  pp.get("rng_seed", m_rng_seed);

  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  if(m_rng_seed < 0) {
    m_rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  m_udist11 = new std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_rng     = new std::mt19937_64(m_rng_seed);
}

void air7_zheleznyak::instantiate_species(){
  m_num_CdrSpecies = 7;
  m_num_RtSpecies = 1;

  m_elec_idx     = 0;
  m_N2plus_idx   = 1;
  m_O2plus_idx   = 2;
  m_N4plus_idx   = 3;
  m_O4plus_idx   = 4;
  m_O2plusN2_idx = 5;
  m_O2minus_idx  = 6;
  
  m_phot_idx         = 0;

  m_CdrSpecies.resize(m_num_CdrSpecies);
  m_CdrSpecies[m_elec_idx]      = RefCountedPtr<CdrSpecies>  (new air7_zheleznyak::electron());
  m_CdrSpecies[m_N2plus_idx]    = RefCountedPtr<CdrSpecies> (new air7_zheleznyak::N2plus());
  m_CdrSpecies[m_O2plus_idx]    = RefCountedPtr<CdrSpecies> (new air7_zheleznyak::O2plus());
  m_CdrSpecies[m_N4plus_idx]    = RefCountedPtr<CdrSpecies> (new air7_zheleznyak::N4plus());
  m_CdrSpecies[m_O4plus_idx]    = RefCountedPtr<CdrSpecies> (new air7_zheleznyak::O4plus());
  m_CdrSpecies[m_O2plusN2_idx]  = RefCountedPtr<CdrSpecies> (new air7_zheleznyak::O2plusN2());
  m_CdrSpecies[m_O2minus_idx]   = RefCountedPtr<CdrSpecies> (new air7_zheleznyak::O2minus());

  m_RtSpecies.resize(m_num_RtSpecies);
  m_RtSpecies[m_phot_idx] = RefCountedPtr<RtSpecies> (new air7_zheleznyak::uv_Photon());

}

void air7_zheleznyak::parseDomainBc(){

  ParmParse pp("air7_zheleznyak");
  std::string str;

  m_wallBc.resize(2*SpaceDim, 0); 
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
	    m_wallBc[idx] = 1;
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
	    m_wallBc[idx] = 1;
	  }
	}
      }
    }
  }
}

void air7_zheleznyak::advance_reaction_network(Vector<Real>&          a_particle_sources,
					       Vector<Real>&          a_Photon_sources,
					       const Vector<Real>     a_particle_densities,
					       const Vector<RealVect> a_particle_gradients,
					       const Vector<Real>     a_Photon_densities,
					       const RealVect         a_E,
					       const RealVect         a_pos,
					       const Real             a_dx,
					       const Real             a_dt,
					       const Real             a_time,
					       const Real             a_kappa) const{
  Vector<Real>     cdr_src(m_num_CdrSpecies, 0.0);
  Vector<Real>     rte_src(m_num_RtSpecies, 0.0);
  Vector<Real>     cdr_phi(m_num_CdrSpecies, 0.0);
  Vector<Real>     rte_phi(m_num_RtSpecies, 0.0);
  Vector<RealVect> cdr_grad(m_num_CdrSpecies, RealVect::Zero);


  // Copy starting data
  for (int i = 0; i < m_num_CdrSpecies; i++){
    cdr_phi[i]  = a_particle_densities[i];
    cdr_grad[i] = a_particle_gradients[i];
  }

  for (int i = 0; i < m_num_RtSpecies; i++){
    rte_phi[i]          = a_Photon_densities[i];
    a_Photon_sources[i] = 0.0;
  }

  int nsteps = ceil(a_dt/m_chemistry_dt);
  Real dt    = a_dt/nsteps;
  Real time  = a_time;
  
  for (int istep = 0; istep < nsteps; istep++){
    if(m_chemistryAlgorithm == chemistryAlgorithm::euler){
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Increment
      for (int i = 0; i < m_num_CdrSpecies; i++){
	cdr_phi[i] = cdr_phi[i] + cdr_src[i]*dt;
      }

      // Add Photons produced in the substep
      for (int i = 0; i < m_num_RtSpecies; i++){
	a_Photon_sources[i] += rte_src[i];
      }
    }
    else if(m_chemistryAlgorithm == chemistryAlgorithm::rk2){

      // Compute slope at k
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      Vector<Real> k1 = cdr_src;

      // Euler update to k+1
      for (int i = 0; i < m_num_CdrSpecies; i++){
	cdr_phi[i] += cdr_src[i]*dt;
	cdr_phi[i] = Max(0.0, cdr_phi[i]);
      }

      // Photons only use the Euler update
      for (int i = 0; i < m_num_RtSpecies; i++){
	a_Photon_sources[i] += rte_src[i];
      }

      // Re-compute slope at k
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Funky notation, but checks out
      for (int i = 0; i < m_num_CdrSpecies; i++){
	cdr_phi[i] = cdr_phi[i] + dt*(0.5*cdr_src[i] - 0.5*k1[i]);
      }
    }
    else if(m_chemistryAlgorithm == chemistryAlgorithm::rk4){
      const Vector<Real> cdr_phi0 = cdr_phi;
      
      // Compute k1 slope
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k1 = cdr_src;

      // Only Euler update for Photons.
      for (int i = 0; i < m_num_RtSpecies; i++){
	a_Photon_sources[i] += rte_src[i];
      }

      // Compute k2 slope
      for (int i = 0; i < m_num_CdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k1[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k2 = cdr_src;

      // Compute k3 slope
      for (int i = 0; i < m_num_CdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k2[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k3 = cdr_src;

      // Compute k4 slope
      for (int i = 0; i < m_num_CdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + dt*k3[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k4 = cdr_src;

      for (int i = 0; i < m_num_CdrSpecies; i++){
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
  for (int i = 0; i < m_num_CdrSpecies; i++){
    a_particle_sources[i] = (cdr_phi[i] - a_particle_densities[i])/a_dt;
  }

}

void air7_zheleznyak::advance_chemistry_euler(Vector<Real>&          a_particle_sources,
					      Vector<Real>&          a_Photon_sources,
					      Vector<Real>&          a_particle_densities,
					      const Vector<RealVect> a_particle_gradients,
					      const Vector<Real>     a_Photon_densities,
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
  const Real R11 = 2.E-13*(300./Te)*Ne*O2p;
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

  // Photoionization, M + y => e + M+
  for (int i = 0; i < a_Photon_densities.size(); i++){
    S_e   += a_Photon_densities[i]/a_dt;
    S_O2p += a_Photon_densities[i]/a_dt;
  }

  // Photon generation
  const Real xfactor = excitation_rates(E)*sergey_factor(m_O2frac)*m_photoi_factor;
  const Real Rgamma  = R1*xfactor;
  const int num_phot = poisson_reaction(Rgamma*volume, a_dt);
  a_Photon_sources[m_phot_idx] = 1.0*num_phot;

  return;
}

int air7_zheleznyak::poisson_reaction(const Real a_propensity, const Real a_dt) const{
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


Vector<Real> air7_zheleznyak::compute_cdr_diffusion_coefficients(const Real         a_time,
								 const RealVect     a_pos,
								 const RealVect     a_E,
								 const Vector<Real> a_cdr_densities) const {

  Vector<Real> dco(m_num_CdrSpecies, 0.0);
  dco[m_elec_idx] = m_e_diffco.get_entry(a_E.vectorLength());
  if(m_isDiffusive_ions){
    dco[m_N2plus_idx]   = m_ion_diffusion;
    dco[m_O2plus_idx]   = m_ion_diffusion;
    dco[m_N4plus_idx]   = m_ion_diffusion;
    dco[m_O4plus_idx]   = m_ion_diffusion;
    dco[m_O2plusN2_idx] = m_ion_diffusion;
    dco[m_O2minus_idx]  = m_ion_diffusion;
  }

  return dco;

}
  
Vector<RealVect> air7_zheleznyak::compute_cdr_velocities(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_E,
							 const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_num_CdrSpecies, RealVect::Zero);

  vel[m_elec_idx] = -a_E*m_e_mobility.get_entry(a_E.vectorLength());
  if(m_isMobile_ions){
    vel[m_N2plus_idx]   =  a_E*m_ion_mobility;
    vel[m_O2plus_idx]   =  a_E*m_ion_mobility;
    vel[m_N4plus_idx]   =  a_E*m_ion_mobility;
    vel[m_O4plus_idx]   =  a_E*m_ion_mobility;
    vel[m_O2plusN2_idx] =  a_E*m_ion_mobility;
    vel[m_O2minus_idx]  = -a_E*m_ion_mobility;
  }
  
  return vel;
}
  
Vector<Real> air7_zheleznyak::compute_cdr_domain_fluxes(const Real           a_time,
							const RealVect       a_pos,
							const int            a_dir,
							const Side::LoHiSide a_side,
							const RealVect       a_E,
							const Vector<Real>   a_cdr_densities,
							const Vector<Real>   a_cdr_velocities,
							const Vector<Real>   a_cdr_gradients,
							const Vector<Real>   a_rte_fluxes,
							const Vector<Real>   a_extrap_cdr_fluxes) const{
  Vector<Real> fluxes(m_num_CdrSpecies, 0.0);

  int idx, sgn;
  if(a_side == Side::Lo){
    sgn = -1;
    idx = 2*a_dir;
  }
  else{
    sgn = 1;
    idx = 2*a_dir + 1;
  }

  if(m_wallBc[idx] == 0){ // Inflow/outflow
    for (int i = 0; i < fluxes.size(); i++){
      fluxes[i] = a_extrap_cdr_fluxes[i];
    }
  }
  else if(m_wallBc[idx] == 1){ // wall
    for (int i = 0; i < fluxes.size(); i++){
      fluxes[i] = 0.0;
    }
  }
  else{
    MayDay::Abort("morrow_jiang::compute_cdr_domain_fluxes - uknown domain bc requested");
  }

  
  return fluxes;
}
  
Vector<Real> air7_zheleznyak::compute_cdr_electrode_fluxes(const Real         a_time,
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

Vector<Real> air7_zheleznyak::compute_cdr_dielectric_fluxes(const Real         a_time,
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

Vector<Real> air7_zheleznyak::compute_cdr_fluxes(const Real         a_time,
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
  Vector<Real> fluxes(m_num_CdrSpecies, 0.0);


  // Drift outflow for now
  for (int i = 0; i < m_num_CdrSpecies; i++){
    fluxes[i] = Max(0.0, a_extrap_cdr_fluxes[i]);
  }

  return fluxes;
}

Real air7_zheleznyak::initial_sigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

Real air7_zheleznyak::compute_alpha(const RealVect a_E) const{
  const Real E     = a_E.vectorLength();
  const Real alpha = m_e_alpha.get_entry(E);
  const Real eta   = m_e_eta.get_entry(E);

  return alpha;
}

Real air7_zheleznyak::excitation_rates(const Real a_E) const{
  const Real Etd = a_E/(m_N*units::s_Td);

  Real y = 1.0;
  if(Etd > 100){
    y = 0.1*exp(233/Etd);
  }

  return y;
}

Real air7_zheleznyak::sergey_factor(const Real a_O2frac) const{
  return 3E-2 + 0.4*pow(a_O2frac, 0.6);
}
#include "CD_NamespaceFooter.H"
