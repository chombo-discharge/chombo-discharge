/*!
  @file   air6_mc8.H
  @brief  6-species (3/3 charged/excited) and 8-photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air6_mc8.H"
#include "air6_mc8_species.H"
#include "units.H" 

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include <ParmParse.H>

std::string air6_mc8::s_bolsig_mobility = "Energy (eV)	Mobility *N (1/m/V/s)";
std::string air6_mc8::s_bolsig_diffco   = "E/N (Td)	Diffusion coefficient *N (1/m/s)";
std::string air6_mc8::s_bolsig_alpha    = "E/N (Td)	Total ionization freq. /N (m3/s)";
std::string air6_mc8::s_bolsig_eta      = "E/N (Td)	Total attachment freq. /N (m3/s)";


air6_mc8::air6_mc8() {

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
  //  parse_domain_bc();

  init_rng();                 // Initialize random number generators
  //  instantiate_species();      // Instantiate species
  //   parse_initial_particles();  // Parse initial particles

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

  m_N = m_p*units::s_Na/(m_T*units::s_R);

}

void air6_mc8::parse_electron_mobility(){
  ParmParse pp("air6_mc8");
  read_file_entries(m_e_mobility, air6_mc8::s_bolsig_mobility);
  m_e_mobility.scale_y(1./m_N); 
  m_e_mobility.make_uniform(m_uniform_entries);
}

void air6_mc8::parse_electron_diffco(){
  ParmParse pp("air6_mc8");
  
  read_file_entries(m_e_diffco, air6_mc8::s_bolsig_mobility);
  m_e_diffco.scale_y(1./m_N); 
  m_e_diffco.make_uniform(m_uniform_entries);
}

void air6_mc8::parse_alpha(){
  ParmParse pp("air6_mc8");
  read_file_entries(m_e_alpha, air6_mc8::s_bolsig_alpha);
  m_e_alpha.scale_y(m_N); 
  m_e_alpha.make_uniform(m_uniform_entries);
}

void air6_mc8::parse_eta(){
  ParmParse pp("air6_mc8");
  read_file_entries(m_e_eta, air6_mc8::s_bolsig_eta);
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

  pp.get("c4v0_X1v0_tau", m_c4v0_X1v0_tau_r);
  pp.get("c4v0_X1v1_tau", m_c4v0_X1v1_tau_r);
  pp.get("c4v1_X1v0_tau", m_c4v1_X1v0_tau_r);
  pp.get("c4v1_X1v1_tau", m_c4v1_X1v1_tau_r);
  pp.get("c4v1_X1v2_tau", m_c4v1_X1v2_tau_r);
  pp.get("c4v1_X1v3_tau", m_c4v1_X1v3_tau_r);
  pp.get("b1v1_X1v0_tau", m_b1v1_X1v0_tau_r);
  pp.get("b1v1_X1v1_tau", m_b1v1_X1v1_tau_r);

  pp.get("c4v0_X1v0_pre", m_c4v0_X1v0_tau_p);
  pp.get("c4v0_X1v1_pre", m_c4v0_X1v1_tau_p);
  pp.get("c4v1_X1v0_pre", m_c4v1_X1v0_tau_p);
  pp.get("c4v1_X1v1_pre", m_c4v1_X1v1_tau_p);
  pp.get("c4v1_X1v2_pre", m_c4v1_X1v2_tau_p);
  pp.get("c4v1_X1v3_pre", m_c4v1_X1v3_tau_p);
  pp.get("b1v1_X1v0_pre", m_b1v1_X1v0_tau_p);
  pp.get("b1v1_X1v1_pre", m_b1v1_X1v1_tau_p);

  pp.get("quenching_photoi_effssure", m_pq);

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

  pp.get("poiss_exp_swap", m_poiss_exp_swap);
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

  m_species[m_elec_idx]  = RefCountedPtr<species>      (new air6_mc8::electron());
  m_species[m_plus_idx]  = RefCountedPtr<species>      (new air6_mc8::M_plus());
  m_species[m_minu_idx]  = RefCountedPtr<species>      (new air6_mc8::M_minus());
  m_species[m_c4v0_idx]  = RefCountedPtr<species>      (new air6_mc8::N2_c4v0());
  m_species[m_c4v1_idx]  = RefCountedPtr<species>      (new air6_mc8::N2_c4v1());
  m_species[m_b1v1_idx]  = RefCountedPtr<species>      (new air6_mc8::N2_b1v1());

  m_photons[m_c4v0_X1v0_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v0_X1v0());
  m_photons[m_c4v0_X1v1_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v0_X1v0());
  m_photons[m_c4v1_X1v0_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v1_X1v0());
  m_photons[m_c4v1_X1v1_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v1_X1v1());
  m_photons[m_c4v1_X1v2_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v1_X1v2());
  m_photons[m_c4v1_X1v3_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_c4v1_X1v3());
  m_photons[m_b1v1_X1v0_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_b1v1_X1v0());
  m_photons[m_b1v1_X1v1_idx] = RefCountedPtr<photon_group> (new air6_mc8::phot_b1v1_X1v1());
}
