/*!
  @file   morrow_lowke.cpp
  @brief  Implementation of morrow_lowke.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Really, really need to revise these functions since we've changed the scaling for the rte equations
*/

#include "morrow_lowke.H"
#include "units.H"

#include <ParmParse.H>
#include <PolyGeom.H>

morrow_lowke::morrow_lowke(){
  CH_TIME("morrow_lowke::morrow_lowke");

  // Number of cdr and rte equations
  m_num_species = 3;
  m_num_photons = 3;


  m_species.resize(m_num_species);
  m_photons.resize(m_num_photons);

  m_nelec_idx   = 0;
  m_nplus_idx   = 1;
  m_nminu_idx   = 2;
  m_photon1_idx = 0;
  m_photon2_idx = 1;
  m_photon3_idx = 2;
  
  m_species[m_nelec_idx]    = RefCountedPtr<species>      (new morrow_lowke::electron());
  m_species[m_nplus_idx]    = RefCountedPtr<species>      (new morrow_lowke::positive_species());
  m_species[m_nminu_idx]    = RefCountedPtr<species>      (new morrow_lowke::negative_species());
  m_photons[m_photon1_idx]  = RefCountedPtr<photon_group> (new morrow_lowke::photon_one());
  m_photons[m_photon2_idx]  = RefCountedPtr<photon_group> (new morrow_lowke::photon_two());
  m_photons[m_photon3_idx]  = RefCountedPtr<photon_group> (new morrow_lowke::photon_three());


  // Default parameters. All of these can be changed through the command line or an input script.
  m_temp      = 300.;      // Gas temperature
  m_fracN2    = 0.8;       // Gas composition
  m_fracO2    = 0.2;       // Gas composition
  m_p         = 1.0;       // Gas pressure (in atm)
  m_pq        = 0.03947;   // Quenching pressure (in atm)
  m_exc_eff   = 0.6;       // Excitation efficient (excitations per collisions)
  m_photo_eff = 0.1;       // Photo-efficiency (emissions per excitation)
  
  m_townsend2_conductor  = 1.E-4; // Second Townsend coefficient on conductor surfaces
  m_townsend2_dielectric = 1.E-4; // Second Townsend coefficient on dielectric surfaces
  m_electrode_yield      = 1.E-6; // Photo-emission yield on conductor surfaces
  m_dielectric_yield     = 1.E-6; // Photo-emission yield on dielectric surfaces
  m_dielectric_work      = 3.0;   // Dielectric work function (in eV)

  m_noise_amp     = 0.0;
  m_noise_freq    = 0.0*RealVect::Unit;
  m_noise_persist = 0.5;
  m_noise_octaves = 1;

  { // Gas composition
    ParmParse pp("morrow_lowke");
    pp.query("gas_temperature",            m_temp);
    pp.query("gas_N2_frac",                m_fracN2);
    pp.query("gas_O2_frac",                m_fracO2);
    pp.query("gas_pressure",               m_p);
    pp.query("gas_quenching_pressure",     m_pq);
    pp.query("excitation_efficiency",      m_exc_eff);
    pp.query("photoionization_efficiency", m_photo_eff);
  }

  { // Electrode things
    ParmParse pp("morrow_lowke");
    pp.query("electrode_townsend2",           m_townsend2_conductor);
    pp.query("electrode_quantum_efficiency",  m_electrode_yield);
    pp.query("dielectric_townsend2",          m_townsend2_dielectric);
    pp.query("dielectric_quantum_efficiency", m_dielectric_yield);
    pp.query("dielectric_work_function",      m_dielectric_work);
  }

  { // Noise parameters for initial dat
    ParmParse pp("morrow_lowke");
    pp.query("noise_amplitude",   m_noise_amp);
    pp.query("noise_octaves",     m_noise_octaves);
    pp.query("noise_persistence", m_noise_persist);
    if(pp.contains("noise_frequency")){
      Vector<Real> freq(SpaceDim);
      pp.queryarr("noise_frequency", freq, 0, SpaceDim);
      m_noise_freq = RealVect(D_DECL(freq[0], freq[1], freq[2]));
    }
  }

  // Initiate noise function and give this to electron and positive species
  m_perlin = RefCountedPtr<perlin_if> (new perlin_if(1.0, m_noise_freq, m_noise_persist, m_noise_octaves));
  morrow_lowke::electron* electron    = static_cast<morrow_lowke::electron*> (&(*m_species[m_nelec_idx]));
  morrow_lowke::positive_species* pos = static_cast<morrow_lowke::positive_species*> (&(*m_species[m_nplus_idx]));
  electron->set_noise(m_perlin);
  pos->set_noise(m_perlin);

  // Convert to correct units and compute necessary things
  m_p  *= units::s_atm2pascal;
  m_pq *= units::s_atm2pascal;

  m_N   = m_p*units::s_Na/(m_temp*units::s_R);
}

morrow_lowke::~morrow_lowke(){


}

Vector<RealVect> morrow_lowke::compute_velocities(const RealVect& a_E) const{
  Vector<RealVect> velocities(m_num_species);
  
  velocities[m_nelec_idx] = this->compute_ve(a_E);
  velocities[m_nplus_idx] = this->compute_vp(a_E);
  velocities[m_nminu_idx] = this->compute_vn(a_E);

  return velocities;
}

RealVect morrow_lowke::compute_ve(const RealVect& a_E) const{
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

RealVect morrow_lowke::compute_vp(const RealVect& a_E) const{
  const RealVect E = a_E*1.E-2;           // E in V/cm
  RealVect vp = 2.34*E*m_p/units::s_atm2pascal;  // Morrow-Lowke wants V/cm
  vp *= 0.01;                             // Morrow-Lowke expression is in cm/s

  return vp;  
}

RealVect morrow_lowke::compute_vn(const RealVect& a_E) const{
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

Vector<Real> morrow_lowke::compute_source_terms(const Vector<Real>& a_species_densities,
						const Vector<Real>& a_photon_densities,
						const RealVect&     a_E) const {
  const Vector<RealVect> vel = this->compute_velocities(a_E); // Does it's own conversion

  const Real alpha  = this->compute_alpha(a_E); // Ionization coefficient
  const Real eta    = this->compute_eta(a_E);   // Attachment coefficient
  const Real beta   = this->compute_beta(a_E);  // Recombination coefficient

  // Cast so we can get A-coefficients
  const morrow_lowke::photon_one*   photon1 = static_cast<morrow_lowke::photon_one*>   (&(*m_photons[m_photon1_idx]));
  const morrow_lowke::photon_two*   photon2 = static_cast<morrow_lowke::photon_two*>   (&(*m_photons[m_photon2_idx]));
  const morrow_lowke::photon_three* photon3 = static_cast<morrow_lowke::photon_three*> (&(*m_photons[m_photon3_idx]));
  
  // Densities and velocities
  const Real Ne  = a_species_densities[m_nelec_idx]; 
  const Real Np  = a_species_densities[m_nplus_idx];
  const Real Nn  = a_species_densities[m_nminu_idx];
  const Real Ve  = vel[m_nelec_idx].vectorLength();
  const Real Vp  = vel[m_nplus_idx].vectorLength();
  const Real Vn  = vel[m_nminu_idx].vectorLength();
  const Real Sph = m_photo_eff*units::s_c0*m_fracO2*m_p*(photon1->get_A()*a_photon_densities[m_photon1_idx]
							 + photon2->get_A()*a_photon_densities[m_photon2_idx]
							 + photon3->get_A()*a_photon_densities[m_photon3_idx]);

  Vector<Real> source(m_num_species, 0.0); 
  Real& Se = source[m_nelec_idx];
  Real& Sp = source[m_nplus_idx];
  Real& Sn = source[m_nminu_idx];

  Se = alpha*Ne*Ve - eta*Ne*Ve   - beta*Ne*Np + Sph;
  Sp = alpha*Ne*Ve - beta*Np*Nn  - beta*Ne*Np + Sph;
  Sn = eta*Ne*Ve   - beta*Np*Nn;

#if 0 // Debug
  Se = 0.0;
  Sp = 0.0;
  Sn = 0.0;
#endif

  return source;
}

Real morrow_lowke::compute_alpha(const RealVect& a_E) const{
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


Real morrow_lowke::compute_eta(const RealVect& a_E) const{

  const Real eta2 = this->compute_eta2(a_E); 
  const Real eta3 = this->compute_eta3(a_E);
  const Real eta  = eta2 + eta3;

  return eta;
}

Real morrow_lowke::compute_eta2(const RealVect& a_E) const{
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
    eta2byN = 6.089E-4*EbyN - 2.983E-19;
  }

  eta2byN *= 1.E-4;   // Morrow-Lowke expression is in cm^2, make it m^2
  eta2 = eta2byN*m_N; //

  return eta2;
}


Real morrow_lowke::compute_eta3(const RealVect& a_E) const{
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

Real morrow_lowke::compute_beta(const RealVect& a_E) const{
  Real beta = 2.0E-7;
  beta *= 1.E-6; // Morrow-Lowke expression is in cm^3. Make it m^3
  return beta;
}

Real morrow_lowke::compute_De(const RealVect& a_E) const{
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

Vector<Real> morrow_lowke::compute_diffusion_coefficients(const RealVect& a_E) const {

  Vector<Real> diffCo(m_num_species, 0.0);
  diffCo[m_nelec_idx] = this->compute_De(a_E);
  diffCo[m_nplus_idx] = 0.;
  diffCo[m_nminu_idx] = 0.;
  
  return diffCo;
}

Vector<Real> morrow_lowke::compute_dielectric_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						     const Vector<Real>& a_ion_densities,
						     const Vector<Real>& a_ion_velocities,
						     const Vector<Real>& a_photon_fluxes,
						     const RealVect&     a_E,
						     const RealVect&     a_pos,
						     const RealVect&     a_normal,
						     const Real&         a_time) const{
  // Outflux of species
  Vector<Real> fluxes(m_num_species, 0.0); 

  if(PolyGeom::dot(a_E, a_normal) > 0.0){ // Field points into gas phase
    fluxes[m_nelec_idx] = Max(0.0, a_extrapolated_fluxes[m_nelec_idx]); // Outflow for electrons
    fluxes[m_nminu_idx] = Max(0.0, a_extrapolated_fluxes[m_nminu_idx]); // Outflow for negative species
  }
  else if(PolyGeom::dot(a_E, a_normal) < 0.0){ // Field points into dielectric
    fluxes[m_nplus_idx] = Max(0.0, a_extrapolated_fluxes[m_nplus_idx]); // Outflow for positive species
  }
  
  // Add in photoelectric effect and ion bombardment for electrons by positive ions
  if(PolyGeom::dot(a_E, a_normal) <= 0.){
    fluxes[m_nelec_idx] += -a_photon_fluxes[m_photon1_idx]*m_dielectric_yield;
    fluxes[m_nelec_idx] += -a_photon_fluxes[m_photon2_idx]*m_dielectric_yield;
    fluxes[m_nelec_idx] += -a_photon_fluxes[m_photon3_idx]*m_dielectric_yield;
    fluxes[m_nelec_idx] += -Max(0.0, a_extrapolated_fluxes[m_nplus_idx])*m_townsend2_dielectric;
  }


#if 0 // This currently does not work.. why?
  // Also add in Schottky emission
  if(PolyGeom::dot(a_E, a_normal) <= 0.){
    const Real W  = m_dielectric_work*units::s_eV;
    const Real dW = sqrt(units::s_Qe*units::s_Qe*units::s_Qe*a_E.vectorLength()/(4.0*units::s_pi*units::s_eps0));
    const Real T  = m_temp;
    const Real A  = 1200000;
    fluxes[m_nelec_idx] += -A*T*T*exp(-(W-dW)/(units::s_kb*T))/units::s_Qe;
  }
#endif


  return fluxes;
}

Vector<Real> morrow_lowke::compute_conductor_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						    const Vector<Real>& a_ion_densities,
						    const Vector<Real>& a_ion_velocities,
						    const Vector<Real>& a_photon_fluxes,
						    const RealVect&     a_E,
						    const RealVect&     a_pos,
						    const RealVect&     a_normal,
						    const Real&         a_time) const{
  Vector<Real> fluxes(m_num_species, 0.0);

  // Treat anode and cathode differently
  const bool is_cathode = PolyGeom::dot(a_E, a_normal) < 0.;
  const bool is_anode   = PolyGeom::dot(a_E, a_normal) > 0.;
  if(is_cathode){
    fluxes = this->compute_cathode_flux(a_extrapolated_fluxes,
					a_ion_densities,
					a_ion_velocities,
					a_photon_fluxes,
					a_E,
					a_pos,
					a_normal,
					a_time);
  }
  else if(is_anode){
    fluxes = this->compute_anode_flux(a_extrapolated_fluxes,
				      a_ion_densities,
				      a_ion_velocities,
				      a_photon_fluxes,
				      a_E,
				      a_pos,
				      a_normal,
				      a_time);
  }

  return fluxes;
}


Vector<Real> morrow_lowke::compute_cathode_flux(const Vector<Real>& a_extrapolated_fluxes,
						const Vector<Real>& a_ion_densities,
						const Vector<Real>& a_ion_velocities,
						const Vector<Real>& a_photon_fluxes,
						const RealVect&     a_E,
						const RealVect&     a_pos,
						const RealVect&     a_normal,
						const Real&         a_time) const{
  Vector<Real> fluxes(m_num_species);

  // Set everything to outflow
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = Max(0., a_extrapolated_fluxes[i]);
  }

  // For electrons, we add ion bombardment of positive ions and the photoelectric effect
  fluxes[m_nelec_idx] = 0.;
  fluxes[m_nelec_idx] += -Max(0., a_extrapolated_fluxes[m_nplus_idx])*m_townsend2_conductor;

  // Photoelectric effect
  fluxes[m_nelec_idx] += -a_photon_fluxes[m_photon1_idx]*m_electrode_yield;
  fluxes[m_nelec_idx] += -a_photon_fluxes[m_photon2_idx]*m_electrode_yield;
  fluxes[m_nelec_idx] += -a_photon_fluxes[m_photon3_idx]*m_electrode_yield;

  return fluxes;
}

Vector<Real> morrow_lowke::compute_anode_flux(const Vector<Real>& a_extrapolated_fluxes,
					      const Vector<Real>& a_ion_densities,
					      const Vector<Real>& a_ion_velocities,
					      const Vector<Real>& a_photon_fluxes,
					      const RealVect&     a_E,
					      const RealVect&     a_pos,
					      const RealVect&     a_normal,
					      const Real&         a_time) const{
  Vector<Real> fluxes(m_num_species);

  // Set to outflux
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = Max(0., a_extrapolated_fluxes[i]);
  }

  return fluxes;
}

Vector<Real> morrow_lowke::compute_rte_source_terms(const Vector<Real>& a_densities, const RealVect& a_E) const{
  Vector<Real> ret(m_num_photons);

  const Vector<RealVect> vel = this->compute_velocities(a_E);      // Compute velocities
  const Real alpha           = this->compute_alpha(a_E);           // Compute ionization coefficient
  const Real Ne              = a_densities[m_nelec_idx];           // Electron density
  const Real ve              = vel[m_nelec_idx].vectorLength();    // Electron velocity
  const Real Se              = Max(0., alpha*Ne*ve);               // Excitations = alpha*Ne*ve

  // Photo emissions = electron excitations * efficiency * quenching
  ret[m_photon1_idx] = Se*m_exc_eff*(m_pq/(m_pq + m_p));
  ret[m_photon2_idx] = Se*m_exc_eff*(m_pq/(m_pq + m_p));
  ret[m_photon3_idx] = Se*m_exc_eff*(m_pq/(m_pq + m_p));

  return ret;
}

Real morrow_lowke::initial_sigma(const Real a_time, const RealVect& a_pos) const{
  return 0.;
}

morrow_lowke::electron::electron(){
  m_name      = "electron";
  m_charge    = -1;
  m_diffusive = true;
  m_unit      = "m-3";

  m_uniform_density = 1.E10;
  m_seed_density    = 0.0;
  m_seed_radius     = 1.0;
  m_noise_density   = 0.0;
  m_seed_pos        = RealVect::Zero;

  { // Get from input script or command line
    ParmParse pp("morrow_lowke");
    std::string str;
    pp.query("uniform_density",    m_uniform_density);
    pp.query("seed_density",       m_seed_density);
    pp.query("seed_radius",        m_seed_radius);
    pp.query("noise_amplitude",    m_noise_density);
    pp.query("electron_diffusion", str);
    if(pp.contains("seed_position")){
      if(str == "true"){
	m_diffusive = true;
      }
      else if(str == "false"){
	m_diffusive = false;
      }
      Vector<Real> pos(SpaceDim);
      pp.queryarr("seed_position", pos, 0, SpaceDim);
      m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
    }
  }
}

morrow_lowke::electron::~electron(){
}

Real morrow_lowke::electron::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);
  const Real noise  = pow(m_perlin->value(a_pos),10)*m_noise_density;;

  return seed + m_uniform_density + noise;
}

void morrow_lowke::electron::set_noise(RefCountedPtr<perlin_if> a_perlin){
  m_perlin = a_perlin;
}

morrow_lowke::positive_species::positive_species(){
  m_name      = "positive_species";
  m_charge    = 1;
  m_diffusive = false;
  m_unit      = "m-3";

  m_uniform_density = 1.E10;
  m_seed_density    = 0.0;
  m_seed_radius     = 1.0;
  m_noise_density   = 0.0;
  m_seed_pos        = RealVect::Zero;
  
  { // Get from input script or command line
    ParmParse pp("morrow_lowke");
    pp.query("uniform_density",  m_uniform_density);
    pp.query("seed_density",     m_seed_density);
    pp.query("seed_radius",      m_seed_radius);
    pp.query("noise_amplitude",  m_noise_density);
    if(pp.contains("seed_position")){
      Vector<Real> pos(SpaceDim);
      pp.queryarr("seed_position", pos, 0, SpaceDim);
      m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
    }
  }
}

morrow_lowke::positive_species::~positive_species(){
}

Real morrow_lowke::positive_species::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);
  const Real noise  = pow(m_perlin->value(a_pos),10)*m_noise_density;;
  
  return seed + m_uniform_density + noise;
}

void morrow_lowke::positive_species::set_noise(RefCountedPtr<perlin_if> a_perlin){
  m_perlin = a_perlin;
}

morrow_lowke::negative_species::negative_species(){
  m_name      = "negative_species";
  m_charge    = -1;
  m_diffusive = false;
  m_unit      = "m-3";
}

morrow_lowke::negative_species::~negative_species(){
}

Real morrow_lowke::negative_species::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.;
}

void morrow_lowke::negative_species::set_noise(RefCountedPtr<perlin_if> a_perlin){
  m_perlin = a_perlin;
}

morrow_lowke::photon_one::photon_one(){
  m_name   = "photon_one";

  m_A      = 1.12E-4;
  m_lambda = 4.15E-2;

  { // Parameters
    ParmParse pp("morrow_lowke");
    pp.query("photon1_A_coeff",      m_A);
    pp.query("photon1_lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.2;
  Real pressure = 1.0;
  {
    ParmParse pp("morrow_lowke");
    pp.query("gas_O2_frac",  O2_frac);
    pp.query("gas_pressure", pressure);
  }
  
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
}

morrow_lowke::photon_one::~photon_one(){
  
}

Real morrow_lowke::photon_one::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.
}

morrow_lowke::photon_two::photon_two(){
  m_name   = "photon_two";

  m_A      = 2.88E-3;
  m_lambda = 1.09E-1;

  { // Parameters
    ParmParse pp("morrow_lowke");
    pp.query("photon2_A_coeff",      m_A);
    pp.query("photon2_lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.2;
  Real pressure = 1.0;
  {
    ParmParse pp("morrow_lowke");
    pp.query("gas_O2_frac",  O2_frac);
    pp.query("gas_pressure", pressure);
  }
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
}

morrow_lowke::photon_two::~photon_two(){
}

Real morrow_lowke::photon_two::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.
}

morrow_lowke::photon_three::photon_three(){
  m_name   = "photon_three";

  m_A      = 2.76E-1;
  m_lambda = 6.69E-1;

  { // Parameters
    ParmParse pp("morrow_lowke");
    pp.query("photon3_A_coeff",      m_A);
    pp.query("photon3_lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.2;
  Real pressure = 1.0;
  {
    ParmParse pp("morrow_lowke");
    pp.query("gas_O2_frac",  O2_frac);
    pp.query("gas_pressure", pressure);
  }
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;  
}

morrow_lowke::photon_three::~photon_three(){
}

Real morrow_lowke::photon_three::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.

}
