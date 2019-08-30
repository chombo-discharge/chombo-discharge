/*!
  @file   morrow_network.cpp
  @brief  Implementation of morrow_network.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Really, really need to revise these functions since we've changed the scaling for the rte equations
*/

#include "morrow_network.H"
#include "units.H"
#include "data_ops.H"

#include <ParmParse.H>
#include <PolyGeom.H>

#include <chrono>

morrow_network::morrow_network(){
  m_num_species = 3;
  m_num_photons = 1;

  m_species.resize(m_num_species);
  m_photons.resize(m_num_photons);

  m_nelec_idx   = 0;
  m_nplus_idx   = 1;
  m_nminu_idx   = 2;
  m_photon1_idx = 0;

  // Instantiate species
  m_species[m_nelec_idx]    = RefCountedPtr<species>      (new morrow_network::electron());
  m_species[m_nplus_idx]    = RefCountedPtr<species>      (new morrow_network::positive_species());
  m_species[m_nminu_idx]    = RefCountedPtr<species>      (new morrow_network::negative_species());
  m_photons[m_photon1_idx]  = RefCountedPtr<photon_group> (new morrow_network::uv_photon());

  // Parse some basic settings
  parse_gas_params();
  parse_see();
  parse_domain_bc();
  parse_reaction_settings();

  // Init RNG
  if(m_seed < 0) {
    m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_rng = new std::mt19937_64(m_seed);

  // Parse initial data
  parse_initial_particles();

  // Convert to correct units and compute necessary things
  m_p  *= units::s_atm2pascal;
  m_pq *= units::s_atm2pascal;
  m_N   = m_p*units::s_Na/(m_T*units::s_R);
}

morrow_network::~morrow_network(){


}

void morrow_network::advance_reaction_network(Vector<Real>&          a_particle_sources,
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

  // Six reactions for this plasma model:
  // ===================================
  // 1. e + M   => 2e + M+ 
  // 2. e + M   => M-
  // 3. e + M+  => 0
  // 4. M- + M+ => 0
  // 5  y + M   => e + M+
  // 6. e + M   => e + M + y
  //
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

void morrow_network::network_tau(Vector<Real>&          a_particle_sources,
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
  const Real volume = pow(a_dx, SpaceDim);

  Vector<int> particle_numbers(m_num_species, 0);
  Vector<int> photon_numbers(m_num_photons, 0);
  
  for (int i = 0; i < m_num_species; i++){
    a_particle_sources[i] = 0.0;
    particle_numbers[i] = floor(a_particle_densities[i]*volume);
  }

  for (int i = 0; i < m_num_photons; i++){
    a_photon_sources[i] = 0.0;
    photon_numbers[i] = floor(a_photon_densities[i]);
  }



  int se = 0;
  int sp = 0;
  int sm = 0;
  int sy = 0;

  // Get some aux stuff for propensity functions
  const RealVect Ve = compute_ve(a_E);
  const Real ve     = Ve.vectorLength();
  const Real alpha  = compute_alpha(a_E);
  const Real eta    = compute_eta(a_E);
  const Real beta   = compute_beta(a_E);
  
  // Reaction 1: e + M => 2e + M+
  const Real a1 = Max(0.0, particle_numbers[m_nelec_idx]*alpha*ve);
  const int S1  = poisson_reaction(a1, a_dt);
  se += S1;
  sp += S1;

  // Reaction 2: e + M   => M-
  const Real a2 = Max(0.0, particle_numbers[m_nelec_idx]*eta*ve);
  const int S2  = poisson_reaction(a2, a_dt);
  se -= S2;
  sm += S2;

  // Reaction 3: e + M+  => 0
  const Real a3 = Max(0.0, particle_numbers[m_nelec_idx]*particle_numbers[m_nplus_idx]*beta/volume);
  const int S3  = poisson_reaction(a3, a_dt);
  se -= S3;
  sp -= S3;

  // Reaction 3: M+ + M-  => 0
  const Real a4 = Max(0.0, particle_numbers[m_nplus_idx]*particle_numbers[m_nminu_idx]*beta/volume);
  const int  S4 = poisson_reaction(a4, a_dt);
  sp -= S4;
  sm -= S4;

  // Reaction 5: y + M   => e + M+
  se += photon_numbers[0];
  sp += photon_numbers[0];

  // Reaction 6: e + M => e + M + y
  const Real a6 = a1*m_exc_eff*m_pq/(m_p+m_pq);
  const int S6  = poisson_reaction(a6, a_dt);
  a_photon_sources[0] = 1.0*S6;

#if 1
  const int res = -se + sp - sm;
  if(res != 0) MayDay::Abort("nope");
#endif

  // Do some scaling
  Real& Se = a_particle_sources[m_nelec_idx];
  Real& Sp = a_particle_sources[m_nplus_idx];
  Real& Sm = a_particle_sources[m_nminu_idx];

  Se = 1.0*se;
  Sp = 1.0*sp;
  Sm = 1.0*sm;

  
  // Do the proper scaling
  const Real factor = 1./(volume*a_dt);
  for (int i = 0; i < a_particle_sources.size(); i++){
    a_particle_sources[i] *= factor;
  }

  


  return;
}

// Tau leaping method
void morrow_network::network_rre(Vector<Real>&          a_particle_sources,
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
  MayDay::Abort("morrow_network::network_rre - not implemented");
}

// Tau leaping method
void morrow_network::network_ssa(Vector<Real>&          a_particle_sources,
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
  MayDay::Abort("morrow_network::network_ssa - not implemented");
}

Vector<RealVect> morrow_network::compute_cdr_velocities(const Real         a_time,
						    const RealVect     a_pos,
						    const RealVect     a_E,
						    const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> velocities(m_num_species);
  
  velocities[m_nelec_idx] = this->compute_ve(a_E);
  velocities[m_nplus_idx] = this->compute_vp(a_E);
  velocities[m_nminu_idx] = this->compute_vn(a_E);

  return velocities;
}

RealVect morrow_network::compute_ve(const RealVect a_E) const{
  RealVect ve = RealVect::Zero;

  const RealVect E = a_E*1.E-2;          // Morrow-Lowke wants E in V/cm
  const Real Emag  = E.vectorLength();   //
  const Real N     = m_N*1.E-6;          // Morrow-Lowke wants N in cm^3
  const Real EbyN  = Emag/N;             //

  const Real lim0 = 2.6E-17;
  const Real lim1 = 1.E-16;
  const Real lim2 = 2.0E-15;
  const Real E_SI = a_E.vectorLength();
  if(EbyN <= lim0){
    ve = -E/Emag*(6.87E22*EbyN + 3.38E4);
  }
  else if(EbyN > lim0 && EbyN <= lim1){
    ve = -E/Emag*(7.293E21*EbyN + 1.63E6);
  }
  else if(EbyN > lim1 && EbyN <= lim2){
    ve = -E/Emag*(1.03E22*EbyN + 1.3E6);
  }
  else if(EbyN > lim2){
    ve = -E/Emag*(7.4E21*EbyN + 7.1E6);
  }

  ve *= 0.01; // Morrow-Lowke expressions are in cm/s
  return ve;
}

RealVect morrow_network::compute_vp(const RealVect a_E) const{
  const RealVect E = a_E*1.E-2;           // E in V/cm
  RealVect vp = 2.34*E*m_p/units::s_atm2pascal;  // Morrow-Lowke wants V/cm
  vp *= 0.01;                             // Morrow-Lowke expression is in cm/s

  return vp;  
}

RealVect morrow_network::compute_vn(const RealVect a_E) const{
  RealVect vn = RealVect::Zero;

  const RealVect E = a_E*1.E-2;       // Morrow-Lowke wants E in V/cm
  const Real Emag = E.vectorLength(); // 
  const Real N    = m_N*1.E-6;        // Morrow-Lowke weants N in cm^3
  const Real EbyN = Emag/N;           //

  const Real lim0 = 5.0E-16;
  if(EbyN <= lim0){
    vn = -2.7*E*m_p/units::s_atm2pascal;
  }
  else{
    vn = -1.86*E*m_p/units::s_atm2pascal;
  }

  vn *= 0.01; // Morrow-Lowke expression is in cm/s
  return vn;
}

Real morrow_network::compute_alpha(const RealVect a_E) const{
  Real alpha    = 0.;
  Real alphabyN = 0.;

  const RealVect E = a_E*1.E-2;        // Morrow-Lowke wants E in V/cm
  const Real Emag  = E.vectorLength(); // 
  const Real N     = m_N*1.E-6;        // Morrow-Lowke wants N in cm^3
  const Real EbyN  = Emag/N;           // This is now V/cm^2

  const Real lim = 1.5E-15;
  if(EbyN > lim){
    alphabyN = 2.0E-16*exp(-7.248E-15/EbyN);
  }
  else{
    alphabyN = 6.619E-17*exp(-5.593E-15/EbyN);
  }

  alphabyN *= 1.E-4; // Morrow-Lowke expression is in cm^2
  alpha = alphabyN*m_N;

  return alpha;
}

Real morrow_network::compute_eta(const RealVect a_E) const{

  const Real eta2 = this->compute_eta2(a_E); 
  const Real eta3 = this->compute_eta3(a_E);
  const Real eta  = eta2 + eta3;

  return eta;
}

Real morrow_network::compute_eta2(const RealVect a_E) const{
  Real eta2    = 0.;
  Real eta2byN = 0.;

  //
  const RealVect E = a_E*1.E-2;        // Morrow-Lowke wants E in V/cm
  const Real Emag  = E.vectorLength(); //
  const Real N     = m_N*1.E-6;        // Morrow-Lowke weants N in cm^3
  const Real EbyN  = Emag/N;           // This is now V/cm^2

  const Real lim = 1.05E-15;
  if(EbyN > lim){
    eta2byN = 8.889E-5*EbyN + 2.567E-19;
  }
  else{
    eta2byN = 6.089E-4*EbyN - 2.893E-19;
  }

  eta2byN *= 1.E-4;   // Morrow-Lowke expression is in cm^2, make it m^2
  eta2 = eta2byN*m_N; //

  return eta2;
}

Real morrow_network::compute_eta3(const RealVect a_E) const{
  const RealVect E = a_E*1.E-2;         // Morrow-Lowke wants E in V/cm
  const Real Emag  = E.vectorLength();  //
  const Real N     = m_N*1.E-6;         // Morrow-Lowke weants N in cm^3
  const Real EbyN  = Emag/N;            // This is now V/cm^2

  //
  Real eta3    = 0.;
  Real eta3byN = 4.7778E-59*pow(EbyN, -1.2749);

  //
  eta3byN *= 1.E-10; // Morrow-Lowke expression is in cm^5. Make it m^5
  eta3 = eta3byN*m_N*m_N;

  return eta3;
}

Real morrow_network::compute_beta(const RealVect a_E) const{
  Real beta = 2.0E-7;
  beta *= 1.E-6; // Morrow-Lowke expression is in cm^3. Make it m^3
  return beta;
}

Real morrow_network::compute_De(const RealVect a_E) const{
  const RealVect E  = a_E*1.E-2;                 // Morrow-Lowke wants E in V/cm
  const Real Emag   = E.vectorLength();          //
  const Real N      = m_N*1.E-6;                 // Morrow-Lowke weants N in cm^3
  const Real EbyN   = Emag/N;                    //

  //
  const RealVect Ve = this->compute_ve(a_E);     // Does it's own conversion, comes out in m/s (aka. mentally sane units)
  const Real ve     = Ve.vectorLength()*1.E-2;   // Make it cm/s
  Real De = 0.3341E9*pow(EbyN, 0.54069)*ve/Emag;
  
  De *= 1.E-4; // Morrow-Lowke expression is in cm^2/s. Make it m^2/s

  return De;
}


Vector<Real> morrow_network::compute_cdr_diffusion_coefficients(const Real         a_time,
								const RealVect     a_pos,
								const RealVect     a_E,
								const Vector<Real> a_cdr_densities) const{
  Vector<Real> diffCo(m_num_species, 0.0);
  diffCo[m_nelec_idx] = this->compute_De(a_E);
  diffCo[m_nplus_idx] = 0.;
  diffCo[m_nminu_idx] = 0.;
  
  return diffCo;
}

Vector<Real> morrow_network::compute_cdr_fluxes(const Real         a_time,
						const RealVect     a_pos,
						const RealVect     a_normal,
						const RealVect     a_E,
						const Vector<Real> a_cdr_densities,
						const Vector<Real> a_cdr_velocities,
						const Vector<Real> a_cdr_gradients,
						const Vector<Real> a_rte_fluxes,
						const Vector<Real> a_extrap_cdr_fluxes,
						const Real         a_townsend2,
						const Real         a_quantum_efficiency) const {
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

Vector<Real> morrow_network::compute_cdr_domain_fluxes(const Real           a_time,
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

  int idx;
  int sgn;
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
    MayDay::Abort("morrow_lowke::compute_cdr_domain_fluxes - uknown domain bc requested");
  }

  
  return fluxes;
}

Vector<Real> morrow_network::compute_cdr_electrode_fluxes(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_normal,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities,
						      const Vector<Real> a_cdr_velocities,
						      const Vector<Real> a_cdr_gradients,
						      const Vector<Real> a_rte_fluxes,
						      const Vector<Real> a_extrap_cdr_fluxes) const {

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
				  a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
}

Vector<Real> morrow_network::compute_cdr_dielectric_fluxes(const Real         a_time,
						       const RealVect     a_pos,
						       const RealVect     a_normal,
						       const RealVect     a_E,
						       const Vector<Real> a_cdr_densities,
						       const Vector<Real> a_cdr_velocities,
						       const Vector<Real> a_cdr_gradients,
						       const Vector<Real> a_rte_fluxes,
						       const Vector<Real> a_extrap_cdr_fluxes) const {

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
				  a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

int morrow_network::poisson_reaction(const Real a_propensity, const Real a_dt) const{
  int value = 0.0;
  const Real mean = a_propensity*a_dt;
  if(mean < m_cutoff_poisson){
    value = round(a_propensity*a_dt);
  }
  else{
    if(mean < m_poiss_exp_swap){
      std::poisson_distribution<int> dist(mean);
      value = dist(*m_rng);
    }
    else{
      std::normal_distribution<double> dist(mean, sqrt(mean));
      value = dist(*m_rng);
    }
  }

  return value;
}

Real morrow_network::initial_sigma(const Real a_time, const RealVect a_pos) const{
  return 0.;
}

morrow_network::electron::electron(){
  m_name      = "electron";
  m_charge    = -1;
  m_diffusive = true;
  m_mobile    = true;
  m_unit      = "m-3";

  Vector<Real> pos(SpaceDim);
  std::string str;

  ParmParse pp("morrow_network");
  pp.get("uniform_density",     m_uniform_density);
  pp.get("seed_density",        m_seed_density);
  pp.get("seed_radius",         m_seed_radius);
  pp.get("diffusive_electrons", str); m_diffusive = (str == "true") ? true : false;

  pp.getarr("seed_position", pos, 0, SpaceDim); m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
}

morrow_network::positive_species::positive_species(){
  m_name      = "positive_species";
  m_charge    = 1;
  m_diffusive = false;
  m_mobile    = false;
  m_unit      = "m-3";

  Vector<Real> pos(SpaceDim);
  std::string str;

  ParmParse pp("morrow_network");
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density",    m_seed_density);
  pp.get("seed_radius",     m_seed_radius);
  pp.get("mobile_ions",     str); m_mobile    = (str == "true") ? true : false;

  pp.getarr("seed_position", pos, 0, SpaceDim); m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
}

morrow_network::negative_species::negative_species(){
  m_name      = "negative_species";
  m_charge    = -1;
  m_diffusive = false;
  m_mobile    = false;
  m_unit      = "m-3";

  std::string str;

  ParmParse pp("morrow_network");
  pp.get("mobile_ions",    str); m_mobile    = (str == "true") ? true : false;
}

morrow_network::electron::~electron(){
  
}

morrow_network::positive_species::~positive_species(){
  
}

morrow_network::negative_species::~negative_species(){
  
}

Real morrow_network::electron::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);

  return m_uniform_density + seed;
}

Real morrow_network::positive_species::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);
  
  return m_uniform_density + seed;
}

Real morrow_network::negative_species::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.;
}

morrow_network::uv_photon::uv_photon(){
  m_name   = "uv_photon";

  Real pressure, O2_frac;
  
  ParmParse pp("morrow_network");
  pp.get("photoi_f1",      m_f1);
  pp.get("photoi_f2",      m_f2);
  pp.get("photoi_K1",      m_K1);
  pp.get("photoi_K2",      m_K2);

  pp.get("gas_O2_frac",  O2_frac);
  pp.get("gas_pressure", pressure);
  pp.get("seed",         m_seed);

  // Convert units
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
  m_K1  = m_K1*m_pO2;
  m_K2  = m_K2*m_pO2;

  // Seed the RNG
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
}

morrow_network::uv_photon::~uv_photon(){
  
}

Real morrow_network::uv_photon::get_kappa(const RealVect a_pos) const {
  MayDay::Abort("morrow_network::uv_photon::get_kappa - should not be called. morrow_network is used with the mc_photo module");
}

Real morrow_network::uv_photon::get_random_kappa() const {
  const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
  return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
}

void morrow_network::parse_gas_params(){
  ParmParse pp("morrow_network");
  std::string str;
  pp.get("gas_temperature",            m_T);
  pp.get("gas_N2_frac",                m_fracN2);
  pp.get("gas_O2_frac",                m_fracO2);
  pp.get("gas_pressure",               m_p);
  pp.get("gas_quenching_pressure",     m_pq);
  pp.get("excitation_efficiency",      m_exc_eff);
  pp.get("photoionization_efficiency", m_photoi_eff);
}

void morrow_network::parse_see(){

  ParmParse pp("morrow_network");

  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void morrow_network::parse_domain_bc(){

  ParmParse pp("morrow_network");
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

void morrow_network::parse_reaction_settings(){
  ParmParse pp("morrow_network");
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
    MayDay::Abort("morrow_network::parse_reaction_settings - stop!");
  }

  pp.get("seed",           m_seed);
  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  pp.get("cutoff_poisson", m_cutoff_poisson);
}

void morrow_network::parse_initial_particles(){

  List<Particle> p;

  // Get types of particles
  add_uniform_particles(p);
  add_gaussian_particles(p);
  
  // Copy initial particles to various species
  m_species[m_nelec_idx]->get_initial_particles() = p;
  m_species[m_nplus_idx]->get_initial_particles() = p;

  // Get the initial deposition scheme
  ParmParse pp("morrow_network");
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
    MayDay::Abort("morrow_network::parse_initial_particles - unknown deposition type requested");
  }
  
  m_species[m_nelec_idx]->get_deposition() = deposition;
  m_species[m_nplus_idx]->get_deposition() = deposition;
}

void morrow_network::add_uniform_particles(List<Particle>& a_particles){
  
  // Get lo/hi sides
  RealVect lo, hi;
  Vector<Real> vec(SpaceDim);
  {
    ParmParse pp1("physical_domain");

    pp1.getarr("lo_corner", vec, 0, SpaceDim); lo = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    pp1.getarr("hi_corner", vec, 0, SpaceDim); hi = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  }


  // Make RNGs
  auto rngX = new std::uniform_real_distribution<Real>(lo[0], hi[0]);
  auto rngY = new std::uniform_real_distribution<Real>(lo[1], hi[1]);
#if CH_SPACEDIM==3
  auto rngZ = new std::uniform_real_distribution<Real>(lo[2], hi[2]);
#endif

  int num_uniform_particles;
  Real num_particles;
  
  // Create uniform particles
  ParmParse pp("morrow_network");
  pp.get("uniform_particles",  num_particles);
  num_uniform_particles = round(num_particles);

  for (int i = 0; i < num_uniform_particles; i++){
    const Real x = (*rngX)(*m_rng);
    const Real y = (*rngY)(*m_rng);
#if CH_SPACEDIM==3
    const Real z = (*rngZ)(*m_rng);
#endif
    RealVect pos = RealVect(D_DECL(x, y, z));
    a_particles.add(Particle(1.0, pos));
  }
}

void morrow_network::add_gaussian_particles(List<Particle>& a_particles){

  // Create uniform particles
  ParmParse pp("morrow_network");

  int num_gaussian_particles;
  Real num_particles;
  Real gaussian_radius;
  RealVect gaussian_center;
  Vector<Real> vec(SpaceDim);
  // Create Gaussian seed particles
  pp.get("gaussian_particles", num_particles);
  pp.get("gaussian_radius",    gaussian_radius);
  pp.getarr("gaussian_center", vec, 0, SpaceDim); gaussian_center = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  num_gaussian_particles = round(num_particles);

  auto rngGX  = new std::normal_distribution<Real>(gaussian_center[0], gaussian_radius);
  auto rngGY  = new std::normal_distribution<Real>(gaussian_center[1], gaussian_radius);
#if CH_SPACEDIM==3
  auto rngGZ  = new std::normal_distribution<Real>(gaussian_center[2], gaussian_radius);
#endif
  for (int i = 0; i < num_gaussian_particles; i++){
    const Real x = (*rngGX)(*m_rng);
    const Real y = (*rngGY)(*m_rng);
#if CH_SPACEDIM==3
    const Real z = (*rngGZ)(*m_rng);
#endif
    RealVect pos = RealVect(D_DECL(x, y, z));
    a_particles.add(Particle(1.0, pos));
  }
}
