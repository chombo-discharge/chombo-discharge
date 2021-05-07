/*!
  @file   air_bolsig.cpp
  @brief  Implementation of air_bolsig.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "air_bolsig.H"
#include "units.H"
#include "perlin_if.H"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ParmParse.H>
#include <PolyGeom.H>

// These are things that we use to identify lines in BOLSIG+ output files. If you want to
// read transport data in a different format, you're on your own. 
std::string air_bolsig::s_skip_fit     = "Fit coefficients y=exp(A+B*ln(x)+C/x+D/x^2+E/x^3)";
std::string air_bolsig::s_bolsig_alpha = "E/N (Td)	Townsend ioniz. coef. alpha/N (m2)";
std::string air_bolsig::s_bolsig_mob   = "E/N (Td)	Mobility *N (1/m/V/s)";
std::string air_bolsig::s_bolsig_diff  = "E/N (Td)	Diffusion coefficient *N (1/m/s)";
std::string air_bolsig::s_bolsig_eta   = "E/N (Td)	Townsend attach. coef. eta/N (m2)";

air_bolsig::air_bolsig(){
  
  // Default parameters
  m_gas_temp = 300.;
  m_frac_O2  = 0.2;
  m_frac_N2  = 0.8;
  m_p        = 1.0;
  m_pq       = 0.03947;

  m_background_rate               = 0.0;
  m_electron_recombination        = 5.E-14;
  m_electron_detachment           = 1.E-18;
  m_ion_recombination             = 2.07E-12;
  m_positive_species_mobility     = 2.E-4;
  m_negative_species_mobility     = 2.E-4;
  m_excitation_efficiency         = 0.6;
  m_photoionization_efficiency    = 0.1;
  m_townsend2_electrode           = 1.E-6;
  m_townsend2_dielectric          = 1.E-6;
  m_electrode_quantum_efficiency  = 1.E-6;
  m_dielectric_quantum_efficiency = 1.E-6;

  m_noise_amplitude   = 0.0;
  m_noise_octaves     = 1;
  m_noise_persistence = 0.5;
  m_noise_frequency   = RealVect::Unit;
      
  { // Get path to BOLSIG and database file. 
    ParmParse pp("air_bolsig");
    pp.get("transport_file",     m_transport_file);
    pp.get("gas_temperature",    m_gas_temp);
    pp.get("gas_O2_frac",        m_frac_O2);
    pp.get("gas_N2_frac",        m_frac_N2);
  }

  { // Get gas parameters
    ParmParse pp("air_bolsig");

    pp.query("gas_pressure",           m_p);
    pp.query("gas_quenching_pressure", m_pq);
  }

  { // Kinetic things not covered by BOLSIG
    ParmParse pp("air_bolsig");
    pp.query("ionization_rate",               m_background_rate);
    pp.query("electron_recombination",        m_electron_recombination);
    pp.query("electron_detachment",           m_electron_detachment);
    pp.query("ion_recombination",             m_ion_recombination);
    pp.query("positive_species_mobility",     m_positive_species_mobility);
    pp.query("negative_species_mobility",     m_negative_species_mobility);
    pp.query("excitation_efficiency",         m_excitation_efficiency);
    pp.query("photoionization_efficiency",    m_photoionization_efficiency);
    pp.query("electrode_townsend2"       ,    m_townsend2_electrode);
    pp.query("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
    pp.query("dielectric_townsend2"       ,   m_townsend2_dielectric);
    pp.query("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
  }

  { // Noise things. Noise function is passed by pointer to species since it must be unique
    ParmParse pp("air_bolsig");
    pp.query("noise_amplitude",   m_noise_amplitude);
    pp.query("noise_octaves",     m_noise_octaves);
    pp.query("noise_persistence", m_noise_persistence);
    if(pp.contains("noise_frequency")){
      Vector<Real> freq(SpaceDim);
      pp.queryarr("noise_frequency", freq, 0, SpaceDim);
      m_noise_frequency = RealVect(D_DECL(freq[0], freq[1], freq[2]));
    }
  }

  { // Transform to SI units
    m_p  *= units::s_atm2pascal;
    m_pq *= units::s_atm2pascal;
    m_N   = m_p*units::s_Na/(m_gas_temp*units::s_R);
  }

  m_num_species = 3;
  m_species.resize(m_num_species);
  m_nelec_idx   = 0;
  m_nplus_idx   = 1;
  m_nminu_idx   = 2;
  m_species[m_nelec_idx]    = RefCountedPtr<species>      (new air_bolsig::electron());
  m_species[m_nplus_idx]    = RefCountedPtr<species>      (new air_bolsig::positive_species());
  m_species[m_nminu_idx]    = RefCountedPtr<species>      (new air_bolsig::negative_species());

  m_num_photons = 3;
  m_photons.resize(m_num_photons);
  m_photon1_idx = 0;
  m_photon2_idx = 1;
  m_photon3_idx = 2;
  m_photons[m_photon1_idx]  = RefCountedPtr<photon_group> (new air_bolsig::photon_one());
  m_photons[m_photon2_idx]  = RefCountedPtr<photon_group> (new air_bolsig::photon_two());
  m_photons[m_photon3_idx]  = RefCountedPtr<photon_group> (new air_bolsig::photon_three());

  // Compute transport data for electrons by calling BOLSIG
  this->compute_transport_coefficients();

  // Normalize things that come out of BOLSIG.
  m_alpha.scale_y(m_N);
  m_eta.scale_y(m_N);
  m_electron_mobility.scale_y(1./m_N);
  m_electron_diffusion.scale_y(1./m_N);

  // Instantiate noise function and pass it down to electron and positive species classes
  m_perlin = RefCountedPtr<perlin_if> (new perlin_if(1.0, m_noise_frequency, m_noise_persistence, m_noise_octaves));
  air_bolsig::electron* electron    = static_cast<air_bolsig::electron*> (&(*m_species[m_nelec_idx]));
  air_bolsig::positive_species* pos = static_cast<air_bolsig::positive_species*> (&(*m_species[m_nplus_idx]));
  electron->set_noise(m_perlin);
  pos->set_noise(m_perlin);
}

air_bolsig::~air_bolsig(){

}

Vector<RealVect> air_bolsig::compute_velocities(const RealVect& a_E) const {
  Vector<RealVect> velocities(m_num_species);

  // Set ion and electron velocityes 
  velocities[m_nplus_idx] = +a_E*m_positive_species_mobility;
  velocities[m_nminu_idx] = -a_E*m_negative_species_mobility;

  // Set the electron velocity
  const Real ET       = a_E.vectorLength()/(units::s_Td*m_N);
  const Real mobility = m_electron_mobility.get_entry(ET);
  velocities[m_nelec_idx] = -a_E*mobility;

  return velocities;
}


Vector<Real> air_bolsig::compute_cdr_source_terms(const Real              a_time,
						  const RealVect&         a_pos,
						  const RealVect&         a_E,
						  const RealVect&         a_gradE,
						  const Vector<Real>&     a_cdr_densities,
						  const Vector<Real>&     a_rte_densities,
						  const Vector<RealVect>& a_grad_cdr) const {
  const Vector<RealVect> vel = this->compute_velocities(a_E);
  const Vector<Real> diffCo  = this->compute_diffusion_coefficients(a_E);

  // Compute ionization and attachment coeffs
  const Real ET    = a_E.vectorLength()/(units::s_Td*m_N);
  const Real alpha = m_alpha.get_entry(ET);
  const Real eta   = m_eta.get_entry(ET);
  const Real zero  = 0.0;

  // Cast so we can get A-coefficients from photons
  const air_bolsig::photon_one*   photon1 = static_cast<air_bolsig::photon_one*>   (&(*m_photons[m_photon1_idx]));
  const air_bolsig::photon_two*   photon2 = static_cast<air_bolsig::photon_two*>   (&(*m_photons[m_photon2_idx]));
  const air_bolsig::photon_three* photon3 = static_cast<air_bolsig::photon_three*> (&(*m_photons[m_photon3_idx]));
  
  // Densities and velocities
  const Real Ne  = a_cdr_densities[m_nelec_idx];
  const Real Np  = a_cdr_densities[m_nplus_idx];
  const Real Nn  = a_cdr_densities[m_nminu_idx];
  const Real ve  = vel[m_nelec_idx].vectorLength();
  const Real De  = diffCo[m_nelec_idx];
  const Real Sph = m_photoionization_efficiency*units::s_c0*m_frac_O2*m_p*(photon1->get_A()*a_rte_densities[m_photon1_idx]
									   + photon2->get_A()*a_rte_densities[m_photon2_idx]
									   + photon3->get_A()*a_rte_densities[m_photon3_idx]);
  
  Vector<Real> source(m_num_species, 0.0); 
  Real& Se = source[m_nelec_idx];
  Real& Sp = source[m_nplus_idx];
  Real& Sn = source[m_nminu_idx];
  
  // Here are the source terms.
  const Real& bep  = m_electron_recombination;
  const Real& bpn  = m_ion_recombination;
  const Real& kdet = m_electron_detachment;

  const Real factor = PolyGeom::dot(a_E,De*a_grad_cdr[m_nelec_idx])/((1.0 + Ne)*PolyGeom::dot(vel[m_nelec_idx], a_E));
  //  const Real alpha_corr = alpha*(1 - Max(factor, 0.0));
  const Real alpha_corr = alpha;
  const Real eta_corr   = eta;
  Se = alpha_corr*Ne*ve - eta_corr*Ne*ve - bep*Ne*Np + m_background_rate + Sph + kdet*Nn*m_N;
  Sp = alpha_corr*Ne*ve - bep*Ne*Np - bpn*Np*Nn + m_background_rate + Sph;
  Sn = eta_corr*Ne*ve   - bpn*Np*Nn - kdet*Nn*m_N;

#if 0 // debug code
  Se = 0.0;
  Sp = 0.0;
  Sn = 0.0;
#endif
  
  return source;
}


Vector<Real> air_bolsig::compute_rte_source_terms(const Vector<Real>& a_cdr_densities, const RealVect& a_E) const {
  Vector<Real> ret(m_num_photons);

  const Real zero            = 0.0;
  const Real ET              = a_E.vectorLength()/(units::s_Td*m_N); // Field in Townsend
  const Vector<RealVect> vel = this->compute_velocities(a_E);        // Compute velocities
  const Real alpha           = m_alpha.get_entry(ET);                // Compute ionization coefficient
  const Real Ne              = a_cdr_densities[m_nelec_idx];         // Electron density
  const Real ve              = vel[m_nelec_idx].vectorLength();      // Electron velocity
  const Real Se              = Max(zero, alpha*Ne*ve);                 // Excitations = alpha*Ne*ve

  // Photosource = electron excitations * efficiency * quenching
  ret[m_photon1_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_photon2_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_photon3_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));

  return ret;
}

Vector<Real> air_bolsig::compute_diffusion_coefficients(const RealVect& a_E) const {
  Vector<Real> diffCo(m_num_species);

  const Real ET    = a_E.vectorLength()/(units::s_Td*m_N);
  
  diffCo[m_nelec_idx] = m_electron_diffusion.get_entry(ET);
  diffCo[m_nplus_idx] = 0.0;
  diffCo[m_nminu_idx] = 0.0;


  return diffCo;
}

Vector<Real> air_bolsig::compute_conductor_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						  const Vector<Real>& a_cdr_densities,
						  const Vector<Real>& a_cdr_velocities,
						  const Vector<Real>& a_rte_fluxes,
						  const RealVect&     a_E,
						  const RealVect&     a_pos,
						  const RealVect&     a_normal,
						  const Real&         a_time) const {
  Vector<Real> fluxes(m_num_species, 0.0); 

  // Treat anode and cathode differently
  const bool is_cathode = PolyGeom::dot(a_E, a_normal) < 0.;
  const bool is_anode   = PolyGeom::dot(a_E, a_normal) > 0.;
  if(is_cathode){
    fluxes = this->compute_cathode_flux(a_extrapolated_fluxes,
					a_cdr_densities,
					a_cdr_velocities,
					a_rte_fluxes,
					a_E,
					a_pos,
					a_normal,
					a_time);
  }
  else if(is_anode){
    fluxes = this->compute_anode_flux(a_extrapolated_fluxes,
				      a_cdr_densities,
				      a_cdr_velocities,
				      a_rte_fluxes,
				      a_E,
				      a_pos,
				      a_normal,
				      a_time);
  }

  return fluxes;
}

Vector<Real> air_bolsig::compute_cathode_flux(const Vector<Real>& a_extrapolated_fluxes,
					      const Vector<Real>& a_cdr_densities,
					      const Vector<Real>& a_cdr_velocities,
					      const Vector<Real>& a_rte_fluxes,
					      const RealVect&     a_E,
					      const RealVect&     a_pos,
					      const RealVect&     a_normal,
					      const Real&         a_time) const {
  Vector<Real> fluxes(m_num_species);

  const Real zero = 0.0;
  // Set everything to outflow
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = Max(zero, a_extrapolated_fluxes[i]);
  }

  // For electrons, we add ion bombardment and the photoelectric effect

  // Ion bombardment
  fluxes[m_nelec_idx] = zero;
  fluxes[m_nelec_idx] += -Max(zero, a_extrapolated_fluxes[m_nplus_idx])*m_townsend2_electrode;

  // Photoelectric effect
  fluxes[m_nelec_idx] += -a_rte_fluxes[m_photon1_idx]*m_electrode_quantum_efficiency;
  fluxes[m_nelec_idx] += -a_rte_fluxes[m_photon2_idx]*m_electrode_quantum_efficiency;
  fluxes[m_nelec_idx] += -a_rte_fluxes[m_photon3_idx]*m_electrode_quantum_efficiency;

  return fluxes;
}

Vector<Real> air_bolsig::compute_anode_flux(const Vector<Real>& a_extrapolated_fluxes,
					    const Vector<Real>& a_cdr_densities,
					    const Vector<Real>& a_cdr_velocities,
					    const Vector<Real>& a_rte_fluxes,
					    const RealVect&     a_E,
					    const RealVect&     a_pos,
					    const RealVect&     a_normal,
					    const Real&         a_time) const {
  Vector<Real> fluxes(m_num_species);
  const Real zero = 0.0;
  // Set to outflux
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = Max(zero, a_extrapolated_fluxes[i]);
  }

  return fluxes;
}


Vector<Real> air_bolsig::compute_dielectric_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						   const Vector<Real>& a_cdr_densities,
						   const Vector<Real>& a_cdr_velocities,
						   const Vector<Real>& a_rte_fluxes,
						   const RealVect&     a_E,
						   const RealVect&     a_pos,
						   const RealVect&     a_normal,
						   const Real&         a_time) const {

  // Outflux of species
  Vector<Real> fluxes(m_num_species, 0.0);
  const Real zero = 0.0;
  
  if(PolyGeom::dot(a_E, a_normal) > 0.0){ // Field points into gas phase
    fluxes[m_nelec_idx] = Max(zero, a_extrapolated_fluxes[m_nelec_idx]); // Outflow for electrons
    fluxes[m_nminu_idx] = Max(zero, a_extrapolated_fluxes[m_nminu_idx]); // Outflow for negative species
  }
  else if(PolyGeom::dot(a_E, a_normal) < 0.0){ // Field points into dielectric
    fluxes[m_nplus_idx] = Max(zero, a_extrapolated_fluxes[m_nplus_idx]); // Outflow for positive species
  }
  
  // Add in photoelectric effect and ion bombardment for electrons by positive ions
  if(PolyGeom::dot(a_E, a_normal) < 0.){ // Field points into dielectric
    fluxes[m_nelec_idx] += -a_rte_fluxes[m_photon1_idx]*m_dielectric_quantum_efficiency;
    fluxes[m_nelec_idx] += -a_rte_fluxes[m_photon2_idx]*m_dielectric_quantum_efficiency;
    fluxes[m_nelec_idx] += -a_rte_fluxes[m_photon3_idx]*m_dielectric_quantum_efficiency;
    fluxes[m_nelec_idx] += -Max(zero, a_extrapolated_fluxes[m_nplus_idx])*m_townsend2_dielectric;
  }


  return fluxes;
}

void air_bolsig::compute_transport_coefficients(){
  Real EbyN;
  Real entry;
  bool readLine = false;
  lookup_table* which_table = NULL;
  std::ifstream infile(m_transport_file);
  std::string line;
  while (std::getline(infile, line)){

    // Right trim string
    line.erase(line.find_last_not_of(" \n\r\t")+1);
      
    // Start reading here
    if(line == air_bolsig::s_bolsig_alpha){
      which_table = &m_alpha;
      readLine   = true;
      continue;
    }
    else if(line == air_bolsig::s_bolsig_eta){
      which_table = &m_eta;
      readLine   = true;
      continue;
    }
    else if(line == air_bolsig::s_bolsig_diff){
      which_table = &m_electron_diffusion;
      readLine   = true;
      continue;
    }
    else if(line == air_bolsig::s_bolsig_mob){
      which_table = &m_electron_mobility;
      readLine   = true;
      continue;
    }

    // Stop when we encounter an empty line
    if((line == "" || line == s_skip_fit) && readLine){
      readLine = false;
      continue;
    }

    if(readLine){
      std::istringstream iss(line);
      if (!(iss >> EbyN >> entry)) {
    	continue;
      }
      which_table->add_entry(EbyN, entry);
    }
  }
  infile.close();
}

void air_bolsig::electron::set_noise(RefCountedPtr<perlin_if> a_perlin){
  m_perlin = a_perlin;
}

void air_bolsig::positive_species::set_noise(RefCountedPtr<perlin_if> a_perlin){
  m_perlin = a_perlin;
}

void air_bolsig::negative_species::set_noise(RefCountedPtr<perlin_if> a_perlin){
  m_perlin = a_perlin;
}

Real air_bolsig::initial_sigma(const Real a_time, const RealVect& a_pos) const {
  return 0.0;
}

Real air_bolsig::electron::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);
  const Real noise  = pow(m_perlin->value(a_pos),10)*m_noise_density;;

  return seed + m_uniform_density + noise;
}

Real air_bolsig::positive_species::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);
  const Real noise  = pow(m_perlin->value(a_pos),10)*m_noise_density;;
  
  return seed + m_uniform_density + noise;
}

Real air_bolsig::negative_species::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.;
}

Real air_bolsig::photon_one::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.
}

Real air_bolsig::photon_two::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.
}

Real air_bolsig::photon_three::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.

}

air_bolsig::electron::electron(){
  m_name      = "electron density";
  m_unit      = "m-3";
  m_charge    = -1;
  m_diffusive = false;

  m_uniform_density = 1.0;
  m_seed_density    = 0.0;
  m_seed_radius     = 1.0;
  m_noise_density   = 0.0;
  m_seed_pos        = RealVect::Zero;

  { // Get from input script or command line
    ParmParse pp("air_bolsig");
    std::string str;
    pp.query("uniform_density", m_uniform_density);
    pp.query("seed_density",    m_seed_density);
    pp.query("seed_radius",     m_seed_radius);
    pp.query("noise_amplitude", m_noise_density);
    pp.query("electron_diffusion", str);
    if(str == "true"){
      m_diffusive = true;
    }
    else if(str == "false"){
      m_diffusive = false;
    }
    if(pp.contains("seed_position")){
      Vector<Real> pos(SpaceDim);
      pp.queryarr("seed_position", pos, 0, SpaceDim);
      m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
    }
  }
}

air_bolsig::electron::~electron(){
  
}

air_bolsig::positive_species::positive_species(){
  m_name      = "positive_species density";
  m_unit      = "m-3";
  m_charge    = 1;
  m_diffusive = false;
  m_mobile    = false;

  m_uniform_density = 1.0;
  m_seed_density    = 0.0;
  m_seed_radius     = 1.0;
  m_noise_density   = 0.0;
  m_seed_pos        = RealVect::Zero;
  
  { // Get from input script or command line
    ParmParse pp("air_bolsig");
    pp.query("uniform_density", m_uniform_density);
    pp.query("seed_density",    m_seed_density);
    pp.query("seed_radius",     m_seed_radius);
    pp.query("noise_amplitude",     m_noise_density);
    if(pp.contains("seed_position")){
      Vector<Real> pos(SpaceDim);
      pp.queryarr("seed_position", pos, 0, SpaceDim);
      m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
    }
  }
}

air_bolsig::positive_species::~positive_species(){
  
}

air_bolsig::negative_species::negative_species(){
  m_name      = "negative_species density";
  m_unit      = "m-3";
  m_charge    = -1;
  m_diffusive = false;
  m_mobile    = false;
}

air_bolsig::negative_species::~negative_species(){
  
}

air_bolsig::photon_one::photon_one(){
  m_name   = "photon_one";

  m_A      = 1.12E-4;
  m_lambda = 4.15E-2;

  { // Parameters
    ParmParse pp("air_bolsig");
    pp.query("photon1_A_coeff",      m_A);
    pp.query("photon1_lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.2;
  Real pressure = 1.0;
  {
    ParmParse pp("air_bolsig");
    pp.query("gas_O2_frac",  O2_frac);
    pp.query("gas_pressure", pressure);
  }
  
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
}

air_bolsig::photon_one::~photon_one(){
  
}

air_bolsig::photon_two::photon_two(){
  m_name   = "photon_two";

  m_A      = 2.88E-3;
  m_lambda = 1.09E-1;

  { // Parameters
    ParmParse pp("air_bolsig");
    pp.query("photon2_A_coeff",      m_A);
    pp.query("photon2_lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.2;
  Real pressure = 1.0;
  {
    ParmParse pp("air_bolsig");
    pp.query("gas_O2_frac",  O2_frac);
    pp.query("gas_pressure", pressure);
  }
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
}

air_bolsig::photon_two::~photon_two(){
  
}

air_bolsig::photon_three::photon_three(){
  m_name   = "photon_three";

  m_A      = 2.76E-1;
  m_lambda = 6.69E-1;

  { // Parameters
    ParmParse pp("air_bolsig");
    pp.query("photon3_A_coeff",      m_A);
    pp.query("photon3_lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.2;
  Real pressure = 1.0;
  {
    ParmParse pp("air_bolsig");
    pp.query("gas_O2_frac",  O2_frac);
    pp.query("gas_pressure", pressure);
  }
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;  
}

air_bolsig::photon_three::~photon_three(){
  
#include "CD_NamespaceFooter.H"


