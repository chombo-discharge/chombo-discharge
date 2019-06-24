/*!
  @file   morrow_fhd.cpp
  @brief  Implementation of morrow_fhd.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Really, really need to revise these functions since we've changed the scaling for the rte equations
*/

#include "morrow_fhd.H"
#include "units.H"
#include "data_ops.H"

#include <ParmParse.H>
#include <PolyGeom.H>

#include <chrono>

morrow_fhd::morrow_fhd(){

  m_num_species = 3;
  m_num_photons = 1;

  m_species.resize(m_num_species);
  m_photons.resize(m_num_photons);

  m_nelec_idx   = 0;
  m_nplus_idx   = 1;
  m_nminu_idx   = 2;
  m_photon1_idx = 0;
  
  m_species[m_nelec_idx]    = RefCountedPtr<species>      (new morrow_fhd::electron());
  m_species[m_nplus_idx]    = RefCountedPtr<species>      (new morrow_fhd::positive_species());
  m_species[m_nminu_idx]    = RefCountedPtr<species>      (new morrow_fhd::negative_species());
  m_photons[m_photon1_idx]  = RefCountedPtr<photon_group> (new morrow_fhd::uv_photon());

  ParmParse pp("morrow_fhd");
  std::string str;
  pp.get("gas_temperature",            m_T);
  pp.get("gas_N2_frac",                m_fracN2);
  pp.get("gas_O2_frac",                m_fracO2);
  pp.get("gas_pressure",               m_p);
  pp.get("gas_quenching_pressure",     m_pq);
  pp.get("excitation_efficiency",      m_exc_eff);
  pp.get("photoionization_efficiency", m_photoi_eff);

  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);

  pp.get("use_fhd", str);   m_fhd = (str == "true") ? true : false;
  pp.get("seed",           m_seed);
  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  pp.get("cutoff_poisson", m_cutoff_poisson);

  // Boundary conditions at domain walls
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

  // Convert to correct units and compute necessary things
  m_p  *= units::s_atm2pascal;
  m_pq *= units::s_atm2pascal;
  m_N   = m_p*units::s_Na/(m_T*units::s_R);

  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);
}

morrow_fhd::~morrow_fhd(){


}

Vector<RealVect> morrow_fhd::compute_cdr_velocities(const Real&         a_time,
						    const RealVect&     a_pos,
						    const RealVect&     a_E,
						    const Vector<Real>& a_cdr_densities) const{
  Vector<RealVect> velocities(m_num_species);
  
  velocities[m_nelec_idx] = this->compute_ve(a_E);
  velocities[m_nplus_idx] = this->compute_vp(a_E);
  velocities[m_nminu_idx] = this->compute_vn(a_E);

  return velocities;
}

RealVect morrow_fhd::compute_ve(const RealVect& a_E) const{
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

RealVect morrow_fhd::compute_vp(const RealVect& a_E) const{
  const RealVect E = a_E*1.E-2;           // E in V/cm
  RealVect vp = 2.34*E*m_p/units::s_atm2pascal;  // Morrow-Lowke wants V/cm
  vp *= 0.01;                             // Morrow-Lowke expression is in cm/s

  return vp;  
}

RealVect morrow_fhd::compute_vn(const RealVect& a_E) const{
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

Vector<Real> morrow_fhd::compute_cdr_source_terms(const Real              a_time,
						  const Real              a_kappa,
						  const Real              a_dx,
						  const RealVect&         a_pos,
						  const RealVect&         a_E,
						  const RealVect&         a_gradE,
						  const Vector<Real>&     a_cdr_densities,
						  const Vector<Real>&     a_rte_densities,
						  const Vector<RealVect>& a_grad_cdr) const {
  Vector<Real> source(m_num_species, 0.0);

  const Real Ne  = a_cdr_densities[m_nelec_idx]; 
  const Real Np  = a_cdr_densities[m_nplus_idx];
  const Real Nn  = a_cdr_densities[m_nminu_idx];
  
  const Real Ve    = (compute_ve(a_E)).vectorLength(); // Electron velocity
  const Real alpha = this->compute_alpha(a_E);         // Ionization coefficient
  const Real eta   = this->compute_eta(a_E);           // Attachment coefficient
  const Real beta  = this->compute_beta(a_E);          // Recombination coefficient
  const Real vol   = pow(a_dx, SpaceDim);
  
  Real products, p;

  // Impact ionization
  p = alpha*Ne*Ve;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_nelec_idx] += products;
  source[m_nplus_idx] += products;

  // Attachment
  p = eta*Ne*Ve;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_nelec_idx] -= products;
  source[m_nminu_idx] += products;

  // Electron-ion recombination
  p = beta*Ne*Np;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_nelec_idx] -= products;
  source[m_nplus_idx] -= products;

  // Ion-ion recombination
  p = beta*Np*Nn;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_nplus_idx] -= products;
  source[m_nminu_idx] -= products;

  // Photoionization
  p = m_photoi_eff*a_rte_densities[m_photon1_idx]/m_dt;
  products = p;//m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_nelec_idx] += products;
  source[m_nplus_idx] += products;

  return source;
}

Real morrow_fhd::compute_alpha(const RealVect& a_E) const{
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


Real morrow_fhd::compute_eta(const RealVect& a_E) const{

  const Real eta2 = this->compute_eta2(a_E); 
  const Real eta3 = this->compute_eta3(a_E);
  const Real eta  = eta2 + eta3;

  return eta;
}

Real morrow_fhd::compute_eta2(const RealVect& a_E) const{
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

Real morrow_fhd::compute_eta3(const RealVect& a_E) const{
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

Real morrow_fhd::compute_beta(const RealVect& a_E) const{
  Real beta = 2.0E-7;
  beta *= 1.E-6; // Morrow-Lowke expression is in cm^3. Make it m^3
  return beta;
}

Real morrow_fhd::compute_De(const RealVect& a_E) const{
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


Vector<Real> morrow_fhd::compute_cdr_diffusion_coefficients(const Real&         a_time,
							    const RealVect&     a_pos,
							    const RealVect&     a_E,
							    const Vector<Real>& a_cdr_densities) const{
  Vector<Real> diffCo(m_num_species, 0.0);
  diffCo[m_nelec_idx] = this->compute_De(a_E);
  diffCo[m_nplus_idx] = 0.;
  diffCo[m_nminu_idx] = 0.;
  
  return diffCo;
}

Vector<Real> morrow_fhd::compute_cdr_fluxes(const Real&         a_time,
					    const RealVect&     a_pos,
					    const RealVect&     a_normal,
					    const RealVect&     a_E,
					    const Vector<Real>& a_cdr_densities,
					    const Vector<Real>& a_cdr_velocities,
					    const Vector<Real>& a_cdr_gradients,
					    const Vector<Real>& a_rte_fluxes,
					    const Vector<Real>& a_extrap_cdr_fluxes,
					    const Real&         a_townsend2,
					    const Real&         a_quantum_efficiency) const {

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

Vector<Real> morrow_fhd::compute_cdr_domain_fluxes(const Real&           a_time,
						   const RealVect&       a_pos,
						   const int&            a_dir,
						   const Side::LoHiSide& a_side,
						   const RealVect&       a_E,
						   const Vector<Real>&   a_cdr_densities,
						   const Vector<Real>&   a_cdr_velocities,
						   const Vector<Real>&   a_cdr_gradients,
						   const Vector<Real>&   a_rte_fluxes,
						   const Vector<Real>&   a_extrap_cdr_fluxes) const{
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

Vector<Real> morrow_fhd::compute_cdr_electrode_fluxes(const Real&         a_time,
						      const RealVect&     a_pos,
						      const RealVect&     a_normal,
						      const RealVect&     a_E,
						      const Vector<Real>& a_cdr_densities,
						      const Vector<Real>& a_cdr_velocities,
						      const Vector<Real>& a_cdr_gradients,
						      const Vector<Real>& a_rte_fluxes,
						      const Vector<Real>& a_extrap_cdr_fluxes) const {

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
				  a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
}

Vector<Real> morrow_fhd::compute_cdr_dielectric_fluxes(const Real&         a_time,
						       const RealVect&     a_pos,
						       const RealVect&     a_normal,
						       const RealVect&     a_E,
						       const Vector<Real>& a_cdr_densities,
						       const Vector<Real>& a_cdr_velocities,
						       const Vector<Real>& a_cdr_gradients,
						       const Vector<Real>& a_rte_fluxes,
						       const Vector<Real>& a_extrap_cdr_fluxes) const {

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
				  a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

Vector<Real> morrow_fhd::compute_rte_source_terms(const Real&         a_time,
						  const Real&         a_kappa,
						  const Real&         a_dx,
						  const RealVect&     a_pos,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities) const{
  Vector<Real> ret(m_num_photons);

  const Real alpha           = this->compute_alpha(a_E);           // Compute ionization coefficient
  const Real Ne              = a_cdr_densities[m_nelec_idx];       // Electron density
  const Real ve              = compute_ve(a_E).vectorLength();     // Electron velocity
  const Real Se              = Max(0., alpha*Ne*ve);               // Excitations = alpha*Ne*ve

  ret[m_photon1_idx] = Se*m_exc_eff*(m_pq/(m_pq + m_p));

  return ret;
}

Real morrow_fhd::stochastic_reaction(const Real a_S, const Real a_vol, const Real a_dt) const{
  Real value = 0.0;
  const Real mean = a_S*a_vol*a_dt;
  if(mean < m_cutoff_poisson){
    value = a_S;
  }
  else{
    if(mean < m_poiss_exp_swap){
      std::poisson_distribution<int> dist(mean);
      value = Max(0.0, 1.0*dist(*m_rng)/(a_vol * a_dt));
    }
    else{
      std::normal_distribution<double> dist(mean, sqrt(mean));
      value = Max(0.0, 1.0*dist(*m_rng)/(a_vol * a_dt));
    }
  }

  return value;
}

Real morrow_fhd::initial_sigma(const Real a_time, const RealVect& a_pos) const{
  return 0.;
}

morrow_fhd::electron::electron(){
  m_name      = "electron";
  m_charge    = -1;
  m_diffusive = true;
  m_mobile    = true;
  m_unit      = "m-3";

  Vector<Real> pos(SpaceDim);
  std::string str;

  ParmParse pp("morrow_fhd");
  pp.get("uniform_density",     m_uniform_density);
  pp.get("seed_density",        m_seed_density);
  pp.get("seed_radius",         m_seed_radius);
  pp.get("diffusive_electrons", str); m_diffusive = (str == "true") ? true : false;

  pp.getarr("seed_position", pos, 0, SpaceDim); m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
}

morrow_fhd::positive_species::positive_species(){
  m_name      = "positive_species";
  m_charge    = 1;
  m_diffusive = false;
  m_mobile    = false;
  m_unit      = "m-3";

  Vector<Real> pos(SpaceDim);
  std::string str;

  ParmParse pp("morrow_fhd");
  pp.get("uniform_density",     m_uniform_density);
  pp.get("seed_density",        m_seed_density);
  pp.get("seed_radius",         m_seed_radius);
  
  pp.get("mobile_ions",    str); m_mobile    = (str == "true") ? true : false;

  pp.getarr("seed_position", pos, 0, SpaceDim); m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
}

morrow_fhd::negative_species::negative_species(){
  m_name      = "negative_species";
  m_charge    = -1;
  m_diffusive = false;
  m_mobile    = false;
  m_unit      = "m-3";

  std::string str;

  ParmParse pp("morrow_fhd");
  pp.get("mobile_ions",    str); m_mobile    = (str == "true") ? true : false;
}

morrow_fhd::electron::~electron(){
  
}

morrow_fhd::positive_species::~positive_species(){
  
}

morrow_fhd::negative_species::~negative_species(){
  
}

Real morrow_fhd::electron::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);

  return m_uniform_density + seed;
}

Real morrow_fhd::positive_species::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);
  
  return m_uniform_density + seed;
}

Real morrow_fhd::negative_species::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.;
}

morrow_fhd::uv_photon::uv_photon(){
  m_name   = "uv_photon";

  Real pressure, O2_frac;
  
  ParmParse pp("morrow_fhd");
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

morrow_fhd::uv_photon::~uv_photon(){
  
}

Real morrow_fhd::uv_photon::get_kappa(const RealVect a_pos) const {
  MayDay::Abort("morrow_fhd::uv_photon::get_kappa - should not be called. morrow_fhd is used with the mc_photo module");
}

Real morrow_fhd::uv_photon::get_random_kappa() const {
  const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
  return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
}
