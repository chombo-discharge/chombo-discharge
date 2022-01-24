/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir3Zheleznyak.cpp
  @brief  Implementation of CD_CdrPlasmaAir3Zheleznyak.H
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
#include <CD_CdrPlasmaAir3Zheleznyak.H>
#include <CD_CdrPlasmaAir3ZheleznyakSpecies.H>
#include <CD_DataOps.H>
#include <CD_DataParser.H>
#include <CD_Units.H> 
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

std::string CdrPlasmaAir3Zheleznyak::s_bolsig_mobility = "# Electron mobility (E/N, mu*N)";
std::string CdrPlasmaAir3Zheleznyak::s_bolsig_diffco   = "# Electron diffusion coefficient (E/N, D*N)";
std::string CdrPlasmaAir3Zheleznyak::s_bolsig_alpha    = "# Townsend ionization coeff (E/N, alpha/N)";
std::string CdrPlasmaAir3Zheleznyak::s_bolsig_eta      = "# Townsend attachment coeff (E/N, eta/N)";

CdrPlasmaAir3Zheleznyak::CdrPlasmaAir3Zheleznyak() {

  initSpecies();      

  parseTransportFile(); 
  parseTransport();
  parseChemistry();
  parseGasParameters();
  parseElectronMobility();
  parseElectronDiffusionCoefficient();
  parseAlpha();
  parseEta();
  parseSEE();
  parseDomainBc();

  initRNG();                 // Initialize random number generators
  
  parse_initial_particles();  // Parse initial particles
}

CdrPlasmaAir3Zheleznyak::~CdrPlasmaAir3Zheleznyak() {

}

void CdrPlasmaAir3Zheleznyak::parseChemistry(){
  ParmParse pp("CdrPlasmaAir3Zheleznyak");

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

void CdrPlasmaAir3Zheleznyak::parseTransportFile(){
  ParmParse pp("CdrPlasmaAir3Zheleznyak");
  pp.get("transport_file",  m_transport_file);
  pp.get("uniform_tables",  m_uniform_entries);
  std::ifstream infile(m_transport_file);
  if(!infile.good()){
    MayDay::Abort("CdrPlasmaAir3Zheleznyak::parseTransportFile - could not find transport data");
  }
  else{
    infile.close();
  }
}

void CdrPlasmaAir3Zheleznyak::parseTransport(){
  ParmParse pp("CdrPlasmaAir3Zheleznyak");

  std::string str;

  pp.get("use_alpha_corr", str);      m_alpha_corr          = (str == "true") ? true : false;
  pp.get("mobile_electrons", str);    m_isMobile_Electrons    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_isDiffusive_Electrons = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);      m_isDiffusive_ions      = (str == "true") ? true : false;
  pp.get("mobile_ions", str);         m_isMobile_ions         = (str == "true") ? true : false;
  
  pp.get("ion_mobility", m_ion_mobility);

  m_ion_diffusion = m_ion_mobility*(Units::kb*m_T)/Units::Qe;
}

void CdrPlasmaAir3Zheleznyak::parseGasParameters(){
  ParmParse pp("CdrPlasmaAir3Zheleznyak");

  // Pressure form input script
  pp.get("pressure",           m_p);
  pp.get("quenching_pressure", m_pq);
  pp.get("temperature",        m_T);
  pp.get("frac_N2",            m_N2frac);
  pp.get("frac_O2",            m_O2frac);
  pp.get("photoi_factor",      m_factor);

  m_p  *= Units::atm2pascal;
  m_pq *= Units::atm2pascal;
  m_N = m_p*Units::Na/(m_T*Units::R);
}

void CdrPlasmaAir3Zheleznyak::parseElectronMobility(){
  m_e_mobility = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Zheleznyak::s_bolsig_mobility, "");
  m_e_mobility.sort();
  m_e_mobility.scale<0>(m_N*Units::Td);
  m_e_mobility.scale<1>(1./m_N); 
  m_e_mobility.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir3Zheleznyak::parseElectronDiffusionCoefficient(){
  m_e_diffco = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Zheleznyak::s_bolsig_diffco, "");
  m_e_diffco.sort();
  m_e_diffco.scale<0>(m_N*Units::Td);
  m_e_diffco.scale<1>(1./m_N); 
  m_e_diffco.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir3Zheleznyak::parseAlpha(){
  m_e_alpha = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Zheleznyak::s_bolsig_alpha, "");
  m_e_alpha.sort();
  m_e_alpha.scale<0>(m_N*Units::Td);
  m_e_alpha.scale<1>(m_N); 
  m_e_alpha.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir3Zheleznyak::parseEta(){
  m_e_eta = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Zheleznyak::s_bolsig_eta, "");
  m_e_eta.sort();
  m_e_eta.scale<0>(m_N*Units::Td);
  m_e_eta.scale<1>(m_N);
  m_e_eta.makeUniform(m_uniform_entries);
}

void CdrPlasmaAir3Zheleznyak::parseSEE(){
  ParmParse pp("CdrPlasmaAir3Zheleznyak");
  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void CdrPlasmaAir3Zheleznyak::initRNG(){
  ParmParse pp("CdrPlasmaAir3Zheleznyak");
  pp.get("rng_seed", m_rng_seed);

  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  if(m_rng_seed < 0) {
    m_rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  m_udist11 = new std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_rng     = new std::mt19937_64(m_rng_seed);
}

void CdrPlasmaAir3Zheleznyak::initSpecies(){
  m_numCdrSpecies = 3;
  m_numRtSpecies = 1;

  m_elec_idx = 0;
  m_plus_idx = 1;
  m_minu_idx = 2;
  m_phot_idx = 0;


  m_cdrSpecies.resize(m_numCdrSpecies);
  m_cdrSpecies[m_elec_idx]  = RefCountedPtr<CdrSpecies>      (new CdrPlasmaAir3Zheleznyak::Electron());
  m_cdrSpecies[m_plus_idx]  = RefCountedPtr<CdrSpecies>      (new CdrPlasmaAir3Zheleznyak::MPlus());
  m_cdrSpecies[m_minu_idx]  = RefCountedPtr<CdrSpecies>      (new CdrPlasmaAir3Zheleznyak::MMinus());

  m_rtSpecies.resize(m_numRtSpecies);
  m_rtSpecies[m_phot_idx] = RefCountedPtr<RtSpecies> (new CdrPlasmaAir3Zheleznyak::uv_Photon());
}

void CdrPlasmaAir3Zheleznyak::parse_initial_particles(){



  // Get some parameters from the input script
  Vector<Real> vec;
  Real weight, uniform_pairs, rad_pairs, gaussian_pairs;
  RealVect center_pairs;

  ParmParse pp("CdrPlasmaAir3Zheleznyak");
  pp.get("particle_weight",            weight);
  pp.get("uniform_pairs",              uniform_pairs);
  pp.get("gaussian_pairs",             gaussian_pairs);
  pp.get("gaussian_pairs_radius",      rad_pairs);
  pp.getarr("gaussian_pairs_center",   vec, 0, SpaceDim); center_pairs   = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  
  List<Particle> Electron_ion_pairs;

  // Add various types of particles
  addUniformParticles(Electron_ion_pairs,  round(uniform_pairs),    weight);
  addGaussianParticles(Electron_ion_pairs, round(gaussian_pairs),   weight, rad_pairs,   center_pairs);
  
  // Set initial particles
  auto& initialElectrons = m_cdrSpecies[m_elec_idx]->getInitialParticles();
  auto& initialIons      = m_cdrSpecies[m_elec_idx]->getInitialParticles();

  initialElectrons = Electron_ion_pairs;
  initialIons      = Electron_ion_pairs;  

  // Set the deposition scheme
  std::string str;
  DepositionType deposition;
  pp.get("particle_deposition", str);
  if(str == "cic"){
    deposition = DepositionType::CIC;
  }
  else if(str == "ngp"){
    deposition = DepositionType::NGP;
  }
  else if(str == "tsc"){
    deposition = DepositionType::TSC;
  }
  else{
    MayDay::Abort("CdrPlasmaAir3Zheleznyak::parse_initial_particles - unknown deposition type requested");
  }
}

void CdrPlasmaAir3Zheleznyak::addUniformParticles(List<Particle>& a_particles, const int a_num, const Real a_weight){

  // Get Lo/Hi sides of domain
  RealVect lo, hi;
  Vector<Real> vec(SpaceDim);
  ParmParse pp("AmrMesh");
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

void CdrPlasmaAir3Zheleznyak::addGaussianParticles(List<Particle>& a_particles,
						   const int       a_num,
						   const Real      a_weight,
						   const Real      a_rad,
						   const RealVect  a_center){
  m_gauss = std::normal_distribution<Real>(0., a_rad);

  for (int i = 0; i < a_num; i++){
    RealVect pos = a_center + randomGaussian();
    a_particles.add(Particle(a_weight, pos));
  }
}

RealVect CdrPlasmaAir3Zheleznyak::randomGaussian(){

  const Real rad = m_gauss(*m_rng);
  return rad*randomDirection();
}

RealVect CdrPlasmaAir3Zheleznyak::randomDirection(){
#if CH_SPACEDIM == 2
  return randomDirection2D();
#else
  return randomDirection3D();
#endif
}

#if CH_SPACEDIM == 2
RealVect CdrPlasmaAir3Zheleznyak::randomDirection2D(){
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
RealVect CdrPlasmaAir3Zheleznyak::randomDirection3D(){
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

void CdrPlasmaAir3Zheleznyak::parseDomainBc(){

  ParmParse pp("CdrPlasmaAir3Zheleznyak");
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

void CdrPlasmaAir3Zheleznyak::advanceReactionNetwork(Vector<Real>&          a_particle_sources,
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

void CdrPlasmaAir3Zheleznyak::advanceChemistryEuler(Vector<Real>&          a_particle_sources,
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
  // R1: e + M -> e + e + M+  alpha*Xe
  // R2: e + M -> M-+         eta*Xe
  // R3: e + M -> c4v0        alpha*Xe*exc_eff(c4v0)
  // R4: e + M -> c4v1        alpha*Xe*exc_eff(c4v1)
  // R5: e + M -> b1v1        alpha*Xe*exc_eff(b1v1)
  const Real volume = pow(a_dx, SpaceDim);
  const Real E      = a_E.vectorLength();
  const Real ve     = E*m_e_mobility.getEntry<1>(E);
  
  // Ionization and attachment coefficients
  Real alpha  = m_e_alpha.getEntry<1>(E);
  Real eta    = m_e_eta.getEntry<1>(E);

  // Modify alpha
  if(m_alpha_corr){
    const RealVect Eunit = a_E/a_E.vectorLength();
    const Real De        = m_e_diffco.getEntry<1>(E);
    const RealVect gNe   = a_particle_gradients[m_elec_idx];

    Real fcorr = 1.0;
    fcorr = 1.0 + PolyGeom::dot(Eunit, De*gNe)/(1.0+a_particle_densities[m_elec_idx]*ve);
    fcorr = Min(fcorr, 1.0);
    fcorr = Max(0.0, fcorr);

    alpha = alpha*fcorr;
  }
  

  const Real R1 = alpha*ve*a_particle_densities[m_elec_idx];
  const Real R2 = eta*ve*a_particle_densities[m_elec_idx];

  Real& Se = a_particle_sources[m_elec_idx];
  Real& Sp = a_particle_sources[m_plus_idx];
  Real& Sm = a_particle_sources[m_minu_idx];

  Se = 0.0;
  Sp = 0.0;
  Sm = 0.0;

  // e + M => e + e + M+
  Se += R1;
  Sp += R1;

  // e + M => M-
  Se -= R2;
  Sm += R2;

  // Photoionization, M + y => e + M+
  for (int i = 0; i < a_Photon_densities.size(); i++){
    Se += a_Photon_densities[i]/a_dt;
    Sp += a_Photon_densities[i]/a_dt;
  }

  const Real quench  = m_pq/(m_p + m_pq);
  const Real xfactor = quench*excitationRates(E)*sergeyFactor(m_O2frac)*m_factor;
  const Real Rgamma  = R1*xfactor;
  const int num_phot = poissonReaction(Rgamma*volume, a_dt);
  a_Photon_sources[m_phot_idx] = 1.0*num_phot;


  return;
}

int CdrPlasmaAir3Zheleznyak::poissonReaction(const Real a_propensity, const Real a_dt) const{
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


Vector<Real> CdrPlasmaAir3Zheleznyak::computeCdrDiffusionCoefficients(const Real         a_time,
								      const RealVect     a_pos,
								      const RealVect     a_E,
								      const Vector<Real> a_cdr_densities) const {

  Vector<Real> dco(m_numCdrSpecies, 0.0);
  dco[m_elec_idx] = m_e_diffco.getEntry<1>(a_E.vectorLength());
  dco[m_plus_idx] = m_ion_diffusion;
  dco[m_minu_idx] = m_ion_diffusion;
  
  return dco;

}
  
Vector<RealVect> CdrPlasmaAir3Zheleznyak::computeCdrDriftVelocities(const Real         a_time,
								    const RealVect     a_pos,
								    const RealVect     a_E,
								    const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_numCdrSpecies, RealVect::Zero);

  vel[m_elec_idx] = -a_E*m_e_mobility.getEntry<1>(a_E.vectorLength());
  vel[m_plus_idx] =  a_E*m_ion_mobility;
  vel[m_minu_idx] = -a_E*m_ion_mobility;
  
  return vel;
}
  
Vector<Real> CdrPlasmaAir3Zheleznyak::computeCdrDomainFluxes(const Real           a_time,
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
  
Vector<Real> CdrPlasmaAir3Zheleznyak::computeCdrElectrodeFluxes(const Real         a_time,
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

Vector<Real> CdrPlasmaAir3Zheleznyak::computeCdrDielectricFluxes(const Real         a_time,
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

Vector<Real> CdrPlasmaAir3Zheleznyak::computeCdrFluxes(const Real         a_time,
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

  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;

  // Switch for setting drift flux to zero for charge species
  Vector<Real> aj(m_numCdrSpecies, 0.0);
  for (int i = 0; i < m_numCdrSpecies; i++){
    if(DataOps::sgn(m_cdrSpecies[i]->getChargeNumber())*PolyGeom::dot(a_E, a_normal) < 0){
      aj[i] = 1.0;
    }
    else {
      aj[i] = 0.0;
    }
  }

  // Drift outflow for now
  for (int i = 0; i < m_numCdrSpecies; i++){
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  return fluxes;
}

Real CdrPlasmaAir3Zheleznyak::initialSigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

Real CdrPlasmaAir3Zheleznyak::computeAlpha(const Real a_E, const RealVect a_pos) const{
  const Real alpha = m_e_alpha.getEntry<1>(a_E);

  return alpha;
}

Real CdrPlasmaAir3Zheleznyak::excitationRates(const Real a_E) const{
  const Real Etd = a_E/(m_N*Units::Td);

  Real y = 1.0;
  if(Etd > 100){
    y = 0.1*exp(233/Etd);
  }

  return y;
}

Real CdrPlasmaAir3Zheleznyak::sergeyFactor(const Real a_O2frac) const{
  return 3E-2 + 0.4*pow(a_O2frac, 0.6);
}

#include <CD_NamespaceFooter.H>
