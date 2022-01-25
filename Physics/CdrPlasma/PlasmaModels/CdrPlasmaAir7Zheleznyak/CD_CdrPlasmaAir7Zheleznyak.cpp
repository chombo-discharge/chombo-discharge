/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir7Zheleznyak.cpp
  @brief  Implementation of CD_CdrPlasmaAir7Zheleznyak.H
  @author Robert Marskar
*/

// Std includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

// Chombo includes
#include <ParmParse.H>
#include <PolyGeom.H>

// Our includes
#include <CD_CdrPlasmaAir7Zheleznyak.H>
#include <CD_CdrPlasmaAir7ZheleznyakSpecies.H>
#include <CD_DataOps.H>
#include <CD_Units.H> 
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

std::string CdrPlasmaAir7Zheleznyak::s_bolsig_energy   = "# Electron mean energy (E/N, eV)";
std::string CdrPlasmaAir7Zheleznyak::s_bolsig_mobility = "# Electron mobility (E/N, mu*N)";
std::string CdrPlasmaAir7Zheleznyak::s_bolsig_diffco   = "# Electron diffusion coefficient (E/N, D*N)";
std::string CdrPlasmaAir7Zheleznyak::s_bolsig_alpha    = "# Townsend ionization coeff (E/N, alpha/N)";
std::string CdrPlasmaAir7Zheleznyak::s_bolsig_eta      = "# Townsend attachment coeff (E/N, eta/N)";
std::string CdrPlasmaAir7Zheleznyak::s_bolsig_alphaN2  = "# N2 ionization (E/N, rate/N)";
std::string CdrPlasmaAir7Zheleznyak::s_bolsig_alphaO2  = "# O2 ionization (E/N, rate/N)";

CdrPlasmaAir7Zheleznyak::CdrPlasmaAir7Zheleznyak() {

  initSpecies();      
  parseTransportFile();
  parseTransport();
  parseChemistry();
  parseGasParameters();
  parseElectronMobility();
  parseElectronDiffusionCoefficient();
  parseAlpha();
  parseEta();
  parseTemperature();
  parseSEE();
  parseDomainBc();
  initRNG();                 // Initialize random number generators
}

CdrPlasmaAir7Zheleznyak::~CdrPlasmaAir7Zheleznyak() {

}

void CdrPlasmaAir7Zheleznyak::readFileEntries(LookupTable<2>& a_table, const std::string a_string){
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
      a_table.addEntry(x, y);
    }
  }
  infile.close();
}

void CdrPlasmaAir7Zheleznyak::parseChemistry(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");

  std::string str;
  pp.get("chemistry_dt", m_chemistry_dt);
  pp.get("chemistry_algorithm", str);

  if(str == "euler"){
    m_ChemistryAlgorithm = ChemistryAlgorithm::euler;
  }
  else if(str == "rk2"){
    m_ChemistryAlgorithm = ChemistryAlgorithm::rk2;
  }
  else if(str == "rk4"){
    m_ChemistryAlgorithm = ChemistryAlgorithm::rk4;
  }
  else{
    MayDay::Abort("air_eed::parseChemistry - unknown chemistry algorithm requested");
  }
}

void CdrPlasmaAir7Zheleznyak::parseTransportFile(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  pp.get("transport_file",  m_transport_file);
  pp.get("uniform_tables",  m_uniform_entries);
  std::ifstream infile(m_transport_file);
  if(!infile.good()){
    MayDay::Abort("CdrPlasmaAir7Zheleznyak::parseTransportFile - could not find transport data");
  }
  else{
    infile.close();
  }
}

void CdrPlasmaAir7Zheleznyak::parseTransport(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");

  std::string str;

  pp.get("use_alpha_corr", str);      m_alpha_corr          = (str == "true") ? true : false;
  pp.get("mobile_electrons", str);    m_isMobile_Electrons    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_isDiffusive_Electrons = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);      m_isDiffusive_ions      = (str == "true") ? true : false;
  pp.get("mobile_ions", str);         m_isMobile_ions         = (str == "true") ? true : false;
  
  pp.get("ion_mobility", m_ion_mobility);

  m_ion_diffusion = m_ion_mobility*(Units::kb*m_T)/Units::Qe;
}

void CdrPlasmaAir7Zheleznyak::parseGasParameters(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");

  // Pressure form input script
  pp.get("pressure",           m_p);
  pp.get("quenching_pressure", m_pq);
  pp.get("temperature",        m_T);
  pp.get("frac_N2",            m_N2frac);
  pp.get("frac_O2",            m_O2frac);
  pp.get("photoi_factor",      m_photoi_factor);

  m_p  *= Units::atm2pascal;
  m_pq *= Units::atm2pascal;
  m_N = m_p*Units::Na/(m_T*Units::R);
}

void CdrPlasmaAir7Zheleznyak::parseElectronMobility(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");

  readFileEntries(m_e_mobility, CdrPlasmaAir7Zheleznyak::s_bolsig_mobility);
  m_e_mobility.scale<0>(m_N*Units::Td);
  m_e_mobility.scale<1>(1./m_N); 
  m_e_mobility.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir7Zheleznyak::parseElectronDiffusionCoefficient(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  
  readFileEntries(m_e_diffco, CdrPlasmaAir7Zheleznyak::s_bolsig_diffco);
  m_e_diffco.scale<0>(m_N*Units::Td);
  m_e_diffco.scale<1>(1./m_N); 
  m_e_diffco.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir7Zheleznyak::parseAlpha(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  readFileEntries(m_e_alpha,   CdrPlasmaAir7Zheleznyak::s_bolsig_alpha);
  readFileEntries(m_e_alphaN2, CdrPlasmaAir7Zheleznyak::s_bolsig_alphaN2);
  readFileEntries(m_e_alphaO2, CdrPlasmaAir7Zheleznyak::s_bolsig_alphaO2);
  
  m_e_alpha.scale<0>(m_N*Units::Td);
  m_e_alphaN2.scale<0>(m_N*Units::Td);
  m_e_alphaO2.scale<0>(m_N*Units::Td);
  
  m_e_alpha.scale<1>(m_N);
  m_e_alphaN2.scale<1>(m_N*m_N2frac);
  m_e_alphaO2.scale<1>(m_N*m_O2frac);
  
  m_e_alpha.makeUniform(m_uniform_entries);
  m_e_alphaN2.makeUniform(m_uniform_entries);
  m_e_alphaO2.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir7Zheleznyak::parseEta(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  readFileEntries(m_e_eta, CdrPlasmaAir7Zheleznyak::s_bolsig_eta);
  m_e_eta.scale<0>(m_N*Units::Td);
  m_e_eta.scale<1>(m_N);
  m_e_eta.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir7Zheleznyak::parseTemperature(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  readFileEntries(m_e_temperature, CdrPlasmaAir7Zheleznyak::s_bolsig_energy);
  m_e_temperature.scale<0>(m_N*Units::Td);
  m_e_temperature.scale<1>(2.0*Units::Qe/(3.0*Units::kb));
  m_e_temperature.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir7Zheleznyak::parseSEE(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void CdrPlasmaAir7Zheleznyak::initRNG(){
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  pp.get("rng_seed", m_rng_seed);

  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  if(m_rng_seed < 0) {
    m_rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  m_udist11 = new std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_rng     = new std::mt19937_64(m_rng_seed);
}

void CdrPlasmaAir7Zheleznyak::initSpecies(){
  m_numCdrSpecies = 7;
  m_numRtSpecies = 1;

  m_elec_idx     = 0;
  m_N2plus_idx   = 1;
  m_O2plus_idx   = 2;
  m_N4plus_idx   = 3;
  m_O4plus_idx   = 4;
  m_O2plusN2_idx = 5;
  m_O2minus_idx  = 6;
  
  m_phot_idx         = 0;

  m_cdrSpecies.resize(m_numCdrSpecies);
  m_cdrSpecies[m_elec_idx]      = RefCountedPtr<CdrSpecies>  (new CdrPlasmaAir7Zheleznyak::Electron());
  m_cdrSpecies[m_N2plus_idx]    = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir7Zheleznyak::N2plus());
  m_cdrSpecies[m_O2plus_idx]    = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir7Zheleznyak::O2plus());
  m_cdrSpecies[m_N4plus_idx]    = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir7Zheleznyak::N4plus());
  m_cdrSpecies[m_O4plus_idx]    = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir7Zheleznyak::O4plus());
  m_cdrSpecies[m_O2plusN2_idx]  = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir7Zheleznyak::O2plusN2());
  m_cdrSpecies[m_O2minus_idx]   = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir7Zheleznyak::O2minus());

  m_rtSpecies.resize(m_numRtSpecies);
  m_rtSpecies[m_phot_idx] = RefCountedPtr<RtSpecies> (new CdrPlasmaAir7Zheleznyak::uv_Photon());

}

void CdrPlasmaAir7Zheleznyak::parseDomainBc(){

  ParmParse pp("CdrPlasmaAir7Zheleznyak");
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

void CdrPlasmaAir7Zheleznyak::advanceReactionNetwork(Vector<Real>&          a_particle_sources,
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
  Vector<Real>     cdr_src(m_numCdrSpecies, 0.0);
  Vector<Real>     rte_src(m_numRtSpecies, 0.0);
  Vector<Real>     cdr_phi(m_numCdrSpecies, 0.0);
  Vector<Real>     rte_phi(m_numRtSpecies, 0.0);
  Vector<RealVect> cdr_grad(m_numCdrSpecies, RealVect::Zero);


  // Copy starting data
  for (int i = 0; i < m_numCdrSpecies; i++){
    cdr_phi[i]  = a_particle_densities[i];
    cdr_grad[i] = a_particle_gradients[i];
  }

  for (int i = 0; i < m_numRtSpecies; i++){
    rte_phi[i]          = a_Photon_densities[i];
    a_Photon_sources[i] = 0.0;
  }

  int nsteps = ceil(a_dt/m_chemistry_dt);
  Real dt    = a_dt/nsteps;
  Real time  = a_time;
  
  for (int istep = 0; istep < nsteps; istep++){
    if(m_ChemistryAlgorithm == ChemistryAlgorithm::euler){
      advanceChemistryEuler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Increment
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi[i] + cdr_src[i]*dt;
      }

      // Add Photons produced in the substep
      for (int i = 0; i < m_numRtSpecies; i++){
	a_Photon_sources[i] += rte_src[i];
      }
    }
    else if(m_ChemistryAlgorithm == ChemistryAlgorithm::rk2){

      // Compute slope at k
      advanceChemistryEuler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      Vector<Real> k1 = cdr_src;

      // Euler update to k+1
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] += cdr_src[i]*dt;
	cdr_phi[i] = Max(0.0, cdr_phi[i]);
      }

      // Photons only use the Euler update
      for (int i = 0; i < m_numRtSpecies; i++){
	a_Photon_sources[i] += rte_src[i];
      }

      // Re-compute slope at k
      advanceChemistryEuler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Funky notation, but checks out
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi[i] + dt*(0.5*cdr_src[i] - 0.5*k1[i]);
      }
    }
    else if(m_ChemistryAlgorithm == ChemistryAlgorithm::rk4){
      const Vector<Real> cdr_phi0 = cdr_phi;
      
      // Compute k1 slope
      advanceChemistryEuler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k1 = cdr_src;

      // Only Euler update for Photons.
      for (int i = 0; i < m_numRtSpecies; i++){
	a_Photon_sources[i] += rte_src[i];
      }

      // Compute k2 slope
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k1[i];
      }
      advanceChemistryEuler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k2 = cdr_src;

      // Compute k3 slope
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k2[i];
      }
      advanceChemistryEuler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k3 = cdr_src;

      // Compute k4 slope
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + dt*k3[i];
      }
      advanceChemistryEuler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k4 = cdr_src;

      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.;
      }
    }
    else{
      MayDay::Abort("advanceReactionNetwork - not supporting other chemistry algorithm yet");
    }

    // Time at next substep
    time += dt;
  }

  // Linearize source term
  for (int i = 0; i < m_numCdrSpecies; i++){
    a_particle_sources[i] = (cdr_phi[i] - a_particle_densities[i])/a_dt;
  }

}

void CdrPlasmaAir7Zheleznyak::advanceChemistryEuler(Vector<Real>&          a_particle_sources,
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
  const Real ve     = E*m_e_mobility.getEntry<1>(E);
  const Real Te     = Max(m_e_temperature.getEntry<1>(E), 300.);

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
    const Real De        = m_e_diffco.getEntry<1>(E);
    const RealVect gNe   = a_particle_gradients[m_elec_idx];

    fcorr = 1.0 + PolyGeom::dot(Eunit, De*gNe)/(1.0+a_particle_densities[m_elec_idx]*ve);
    fcorr = Min(fcorr, 1.0);
    fcorr = Max(0.0, fcorr);
  }

  const Real R1  = fcorr*m_e_alphaN2.getEntry<1>(E)*Ne;
  const Real R2  = fcorr*m_e_alphaO2.getEntry<1>(E)*Ne;
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
  const Real xfactor = excitationRates(E)*sergeyFactor(m_O2frac)*m_photoi_factor;
  const Real Rgamma  = R1*xfactor;
  const int num_phot = poissonReaction(Rgamma*volume, a_dt);
  a_Photon_sources[m_phot_idx] = 1.0*num_phot;

  return;
}

int CdrPlasmaAir7Zheleznyak::poissonReaction(const Real a_propensity, const Real a_dt) const{
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


Vector<Real> CdrPlasmaAir7Zheleznyak::computeCdrDiffusionCoefficients(const Real         a_time,
								      const RealVect     a_pos,
								      const RealVect     a_E,
								      const Vector<Real> a_cdr_densities) const {

  Vector<Real> dco(m_numCdrSpecies, 0.0);
  dco[m_elec_idx] = m_e_diffco.getEntry<1>(a_E.vectorLength());
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
  
Vector<RealVect> CdrPlasmaAir7Zheleznyak::computeCdrDriftVelocities(const Real         a_time,
								    const RealVect     a_pos,
								    const RealVect     a_E,
								    const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_numCdrSpecies, RealVect::Zero);

  vel[m_elec_idx] = -a_E*m_e_mobility.getEntry<1>(a_E.vectorLength());
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
  
Vector<Real> CdrPlasmaAir7Zheleznyak::computeCdrDomainFluxes(const Real           a_time,
							     const RealVect       a_pos,
							     const int            a_dir,
							     const Side::LoHiSide a_side,
							     const RealVect       a_E,
							     const Vector<Real>   a_cdr_densities,
							     const Vector<Real>   a_cdr_velocities,
							     const Vector<Real>   a_cdr_gradients,
							     const Vector<Real>   a_rte_fluxes,
							     const Vector<Real>   a_extrap_cdr_fluxes) const{
  Vector<Real> fluxes(m_numCdrSpecies, 0.0);

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
    MayDay::Abort("morrow_jiang::computeCdrDomainFluxes - uknown domain bc requested");
  }

  
  return fluxes;
}
  
Vector<Real> CdrPlasmaAir7Zheleznyak::computeCdrElectrodeFluxes(const Real         a_time,
								const RealVect     a_pos,
								const RealVect     a_normal,
								const RealVect     a_E,
								const Vector<Real> a_cdr_densities,
								const Vector<Real> a_cdr_velocities,
								const Vector<Real> a_cdr_gradients,
								const Vector<Real> a_rte_fluxes,
								const Vector<Real> a_extrap_cdr_fluxes) const{
  return computeCdrFluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			  a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
}

Vector<Real> CdrPlasmaAir7Zheleznyak::computeCdrDielectricFluxes(const Real         a_time,
								 const RealVect     a_pos,
								 const RealVect     a_normal,
								 const RealVect     a_E,
								 const Vector<Real> a_cdr_densities,
								 const Vector<Real> a_cdr_velocities,
								 const Vector<Real> a_cdr_gradients,
								 const Vector<Real> a_rte_fluxes,
								 const Vector<Real> a_extrap_cdr_fluxes) const{
  return computeCdrFluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			  a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

Vector<Real> CdrPlasmaAir7Zheleznyak::computeCdrFluxes(const Real         a_time,
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
  Vector<Real> fluxes(m_numCdrSpecies, 0.0);


  // Drift outflow for now
  for (int i = 0; i < m_numCdrSpecies; i++){
    fluxes[i] = Max(0.0, a_extrap_cdr_fluxes[i]);
  }

  return fluxes;
}

Real CdrPlasmaAir7Zheleznyak::initialSigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

Real CdrPlasmaAir7Zheleznyak::computeAlpha(const Real a_E, const RealVect a_pos) const{
  return m_e_alpha.getEntry<1>(a_E);
}

Real CdrPlasmaAir7Zheleznyak::excitationRates(const Real a_E) const{
  const Real Etd = a_E/(m_N*Units::Td);

  Real y = 1.0;
  if(Etd > 100){
    y = 0.1*exp(233/Etd);
  }

  return y;
}

Real CdrPlasmaAir7Zheleznyak::sergeyFactor(const Real a_O2frac) const{
  return 3E-2 + 0.4*pow(a_O2frac, 0.6);
}

#include <CD_NamespaceFooter.H>
