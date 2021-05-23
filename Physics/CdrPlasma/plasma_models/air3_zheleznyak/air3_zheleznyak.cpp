/*!
  @file   air3_zheleznyak.H
  @brief  3-species and 8-Photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air3_zheleznyak.H"
#include "air3_zheleznyak_species.H"
#include <CD_DataOps.H>
#include <CD_Units.H> 

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include <ParmParse.H>
#include <PolyGeom.H>

#include "CD_NamespaceHeader.H"
using namespace Physics::CdrPlasma;

std::string air3_zheleznyak::s_bolsig_mobility = "# Electron mobility (E/N, mu*N)";
std::string air3_zheleznyak::s_bolsig_diffco   = "# Electron diffusion coefficient (E/N, D*N)";
std::string air3_zheleznyak::s_bolsig_alpha    = "# Townsend ionization coeff (E/N, alpha/N)";
std::string air3_zheleznyak::s_bolsig_eta      = "# Townsend attachment coeff (E/N, eta/N)";

air3_zheleznyak::air3_zheleznyak() {

  instantiate_species();      

  parse_transport_file(); 
  parse_transport();
  parse_chemistry();
  parse_gas_params();
  parse_electron_mobility();
  parse_electron_diffco();
  parse_alpha();
  parse_eta();
  parse_see();
  parseDomainBc();

  init_rng();                 // Initialize random number generators
  
  parse_initial_particles();  // Parse initial particles
}

air3_zheleznyak::~air3_zheleznyak() {

}

void air3_zheleznyak::read_file_entries(LookupTable& a_table, const std::string a_string){
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

void air3_zheleznyak::parse_chemistry(){
  ParmParse pp("air3_zheleznyak");

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

void air3_zheleznyak::parse_transport_file(){
  ParmParse pp("air3_zheleznyak");
  pp.get("transport_file",  m_transport_file);
  pp.get("uniform_tables",  m_uniform_entries);
  std::ifstream infile(m_transport_file);
  if(!infile.good()){
    MayDay::Abort("air3_zheleznyak::parse_transport_file - could not find transport data");
  }
  else{
    infile.close();
  }
}

void air3_zheleznyak::parse_transport(){
  ParmParse pp("air3_zheleznyak");

  std::string str;

  pp.get("use_alpha_corr", str);      m_alpha_corr          = (str == "true") ? true : false;
  pp.get("mobile_electrons", str);    m_isMobile_electrons    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_isDiffusive_electrons = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);      m_isDiffusive_ions      = (str == "true") ? true : false;
  pp.get("mobile_ions", str);         m_isMobile_ions         = (str == "true") ? true : false;
  
  pp.get("ion_mobility", m_ion_mobility);

  m_ion_diffusion = m_ion_mobility*(Units::kb*m_T)/Units::Qe;
}

void air3_zheleznyak::parse_gas_params(){
  ParmParse pp("air3_zheleznyak");

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

void air3_zheleznyak::parse_electron_mobility(){
  ParmParse pp("air3_zheleznyak");

  read_file_entries(m_e_mobility, air3_zheleznyak::s_bolsig_mobility);
  m_e_mobility.scaleX(m_N*Units::Td);
  m_e_mobility.scaleY(1./m_N); 
  m_e_mobility.makeUniform(m_uniform_entries);
}

void air3_zheleznyak::parse_electron_diffco(){
  ParmParse pp("air3_zheleznyak");
  
  read_file_entries(m_e_diffco, air3_zheleznyak::s_bolsig_diffco);
  m_e_diffco.scaleX(m_N*Units::Td);
  m_e_diffco.scaleY(1./m_N); 
  m_e_diffco.makeUniform(m_uniform_entries);
}

void air3_zheleznyak::parse_alpha(){
  ParmParse pp("air3_zheleznyak");
  read_file_entries(m_e_alpha, air3_zheleznyak::s_bolsig_alpha);
  m_e_alpha.scaleX(m_N*Units::Td);
  m_e_alpha.scaleY(m_N); 
  m_e_alpha.makeUniform(m_uniform_entries);
}

void air3_zheleznyak::parse_eta(){
  ParmParse pp("air3_zheleznyak");
  read_file_entries(m_e_eta, air3_zheleznyak::s_bolsig_eta);
  m_e_eta.scaleX(m_N*Units::Td);
  m_e_eta.scaleY(m_N);
  m_e_eta.makeUniform(m_uniform_entries);
}

void air3_zheleznyak::parse_see(){
  ParmParse pp("air3_zheleznyak");
  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void air3_zheleznyak::init_rng(){
  ParmParse pp("air3_zheleznyak");
  pp.get("rng_seed", m_rng_seed);

  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  if(m_rng_seed < 0) {
    m_rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  m_udist11 = new std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_rng     = new std::mt19937_64(m_rng_seed);
}

void air3_zheleznyak::instantiate_species(){
  m_numCdrSpecies = 3;
  m_numRtSpecies = 1;

  m_elec_idx = 0;
  m_plus_idx = 1;
  m_minu_idx = 2;
  m_phot_idx = 0;


  m_CdrSpecies.resize(m_numCdrSpecies);
  m_CdrSpecies[m_elec_idx]  = RefCountedPtr<CdrSpecies>      (new air3_zheleznyak::electron());
  m_CdrSpecies[m_plus_idx]  = RefCountedPtr<CdrSpecies>      (new air3_zheleznyak::M_plus());
  m_CdrSpecies[m_minu_idx]  = RefCountedPtr<CdrSpecies>      (new air3_zheleznyak::M_minus());

  m_RtSpecies.resize(m_numRtSpecies);
  m_RtSpecies[m_phot_idx] = RefCountedPtr<RtSpecies> (new air3_zheleznyak::uv_Photon());
}

void air3_zheleznyak::parse_initial_particles(){



  // Get some parameters from the input script
  Vector<Real> vec;
  Real weight, uniform_pairs, rad_pairs, gaussian_pairs;
  RealVect center_pairs;

  ParmParse pp("air3_zheleznyak");
  pp.get("particle_weight",            weight);
  pp.get("uniform_pairs",              uniform_pairs);
  pp.get("gaussian_pairs",             gaussian_pairs);
  pp.get("gaussian_pairs_radius",      rad_pairs);
  pp.getarr("gaussian_pairs_center",   vec, 0, SpaceDim); center_pairs   = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  
  List<Particle> electron_ion_pairs;

  // Add various types of particles
  add_uniform_particles(electron_ion_pairs,  round(uniform_pairs),    weight);
  add_gaussian_particles(electron_ion_pairs, round(gaussian_pairs),   weight, rad_pairs,   center_pairs);
  
  // Set initial particles
  m_CdrSpecies[m_elec_idx]->getInitialParticles() = electron_ion_pairs;
  m_CdrSpecies[m_plus_idx]->getInitialParticles() = electron_ion_pairs;

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
    MayDay::Abort("air3_zheleznyak::parse_initial_particles - unknown deposition type requested");
  }
  
  for (int i = 0; i < m_numCdrSpecies; i++){
    m_CdrSpecies[i]->getDeposition() = deposition;
  }
}

void air3_zheleznyak::add_uniform_particles(List<Particle>& a_particles, const int a_num, const Real a_weight){

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

void air3_zheleznyak::add_gaussian_particles(List<Particle>& a_particles,
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

RealVect air3_zheleznyak::randomGaussian(){

  const Real rad = m_gauss(*m_rng);
  return rad*randomDirection();
}

RealVect air3_zheleznyak::randomDirection(){
#if CH_SPACEDIM == 2
  return randomDirection2D();
#else
  return randomDirection3D();
#endif
}

#if CH_SPACEDIM == 2
RealVect air3_zheleznyak::randomDirection2D(){
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
RealVect air3_zheleznyak::randomDirection3D(){
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

void air3_zheleznyak::parseDomainBc(){

  ParmParse pp("air3_zheleznyak");
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

void air3_zheleznyak::advanceReactionNetwork(Vector<Real>&          a_particle_sources,
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
    if(m_chemistryAlgorithm == chemistryAlgorithm::euler){
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Increment
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi[i] + cdr_src[i]*dt;
      }

      // Add Photons produced in the substep
      for (int i = 0; i < m_numRtSpecies; i++){
	a_Photon_sources[i] += rte_src[i];
      }
    }
    else if(m_chemistryAlgorithm == chemistryAlgorithm::rk2){

      // Compute slope at k
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
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
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);

      // Funky notation, but checks out
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi[i] + dt*(0.5*cdr_src[i] - 0.5*k1[i]);
      }
    }
    else if(m_chemistryAlgorithm == chemistryAlgorithm::rk4){
      const Vector<Real> cdr_phi0 = cdr_phi;
      
      // Compute k1 slope
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k1 = cdr_src;

      // Only Euler update for Photons.
      for (int i = 0; i < m_numRtSpecies; i++){
	a_Photon_sources[i] += rte_src[i];
      }

      // Compute k2 slope
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k1[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k2 = cdr_src;

      // Compute k3 slope
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + 0.5*dt*k2[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
      const Vector<Real> k3 = cdr_src;

      // Compute k4 slope
      for (int i = 0; i < m_numCdrSpecies; i++){
	cdr_phi[i] = cdr_phi0[i] + dt*k3[i];
      }
      advance_chemistry_euler(cdr_src, rte_src, cdr_phi, a_particle_gradients, rte_phi, a_E, a_pos, a_dx, dt, time, a_kappa);
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

void air3_zheleznyak::advance_chemistry_euler(Vector<Real>&          a_particle_sources,
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
  const Real ve     = E*m_e_mobility.getEntry(E);
  
  // Ionization and attachment coefficients
  Real alpha  = m_e_alpha.getEntry(E);
  Real eta    = m_e_eta.getEntry(E);

  // Modify alpha
  if(m_alpha_corr){
    const RealVect Eunit = a_E/a_E.vectorLength();
    const Real De        = m_e_diffco.getEntry(E);
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
  const Real xfactor = quench*excitation_rates(E)*sergey_factor(m_O2frac)*m_factor;
  const Real Rgamma  = R1*xfactor;
  const int num_phot = poisson_reaction(Rgamma*volume, a_dt);
  a_Photon_sources[m_phot_idx] = 1.0*num_phot;


  return;
}

int air3_zheleznyak::poisson_reaction(const Real a_propensity, const Real a_dt) const{
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


Vector<Real> air3_zheleznyak::computeCdrDiffusionCoefficients(const Real         a_time,
								 const RealVect     a_pos,
								 const RealVect     a_E,
								 const Vector<Real> a_cdr_densities) const {

  Vector<Real> dco(m_numCdrSpecies, 0.0);
  dco[m_elec_idx] = m_e_diffco.getEntry(a_E.vectorLength());
  dco[m_plus_idx] = m_ion_diffusion;
  dco[m_minu_idx] = m_ion_diffusion;
  
  return dco;

}
  
Vector<RealVect> air3_zheleznyak::computeCdrDriftVelocities(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_E,
							 const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_numCdrSpecies, RealVect::Zero);

  vel[m_elec_idx] = -a_E*m_e_mobility.getEntry(a_E.vectorLength());
  vel[m_plus_idx] =  a_E*m_ion_mobility;
  vel[m_minu_idx] = -a_E*m_ion_mobility;
  
  return vel;
}
  
Vector<Real> air3_zheleznyak::computeCdrDomainFluxes(const Real           a_time,
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
  
Vector<Real> air3_zheleznyak::computeCdrElectrodeFluxes(const Real         a_time,
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

Vector<Real> air3_zheleznyak::computeCdrDielectricFluxes(const Real         a_time,
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

Vector<Real> air3_zheleznyak::computeCdrFluxes(const Real         a_time,
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
    if(DataOps::sgn(m_CdrSpecies[i]->getChargeNumber())*PolyGeom::dot(a_E, a_normal) < 0){
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

Real air3_zheleznyak::initialSigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

Real air3_zheleznyak::computeAlpha(const RealVect a_E) const{
  const Real E     = a_E.vectorLength();
  const Real alpha = m_e_alpha.getEntry(E);
  const Real eta   = m_e_eta.getEntry(E);

  return alpha;
}

Real air3_zheleznyak::excitation_rates(const Real a_E) const{
  const Real Etd = a_E/(m_N*Units::Td);

  Real y = 1.0;
  if(Etd > 100){
    y = 0.1*exp(233/Etd);
  }

  return y;
}

Real air3_zheleznyak::sergey_factor(const Real a_O2frac) const{
  return 3E-2 + 0.4*pow(a_O2frac, 0.6);
}
#include "CD_NamespaceFooter.H"
