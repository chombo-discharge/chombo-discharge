/*!
  @file   air6_mc8.H
  @brief  6-species (3/3 charged/excited) and 8-photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air6_mc8.H"
#include "units.H" 

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include <ParmParse.H>

air6_mc8::air6_mc8() {

  parse_transport_file();
  parse_transport();
  //  parse_gas_params();
  //  parse_electron_mobility();
  //  parse_alpha();
  //  parse_eta();
  
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
