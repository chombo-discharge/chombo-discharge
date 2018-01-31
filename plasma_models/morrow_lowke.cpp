/*!
  @file   morrow_lowke.H
  @brief  Implementation of morrow_lowke.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "morrow_lowke.H"
#include "units.H"

#include <ParmParse.H>

morrow_lowke::morrow_lowke(){
  CH_TIME("morrow_lowke::morrow_lowke");

  // Number of cdr and rte equations
  const int num_species = 3;
  const int num_photons = 3;

  m_species.resize(num_species);
  m_photons.resize(num_photons);

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
  m_fracN2    = 0.2;       // Gas composition
  m_fracO2    = 0.8;       // Gas composition
  m_p         = 1.0;       // Gas pressure (in atm)
  m_pq        = 0.03947;   // Quenching pressure (in atm)
  m_exc_eff   = 0.6;       // Excitation efficient (excitations per collisions)
  m_photo_eff = 0.1;       // Photo-efficiency (emissions per excitation)
  
  m_townsend2_conductor  = 1.E-4; // Second Townsend coefficient on conductor surfaces
  m_townsend2_dielectric = 1.E-4; // Second Townsend coefficient on dielectric surfaces
  m_electrode_yield      = 1.E-3; // Photo-emission yield on conductor surfaces
  m_dielectric_yield     = 1.E-6; // Photo-emission yield on dielectric surfaces
  m_dielectric_work      = 3.0;   // Dielectric work function (in eV)

  m_noise_amp     = 0.0;
  m_noise_freq    = 0.0*RealVect::Unit;
  m_noise_persist = 0.5;
  m_noise_octaves = 1;

  { // Gas composition
    ParmParse pp("gas");
    pp.query("temperature",        m_temp);
    pp.query("frac_N2",            m_fracN2);
    pp.query("frac_O2",            m_fracO2);
    pp.query("pressure",           m_p);
    pp.query("quenching_pressure", m_pq);
    pp.query("excitation_eff",     m_exc_eff);
    pp.query("photo_eff",          m_photo_eff);
  }

  { // Electrode things
    ParmParse pp("electrode");
    pp.query("townsend_coeff", m_townsend2_conductor);
    pp.query("quantum_eff",    m_electrode_yield);
  }

  { // Dielectric things
    ParmParse pp("dielectric");
    pp.query("townsend_coeff", m_townsend2_dielectric);
    pp.query("quantum_eff",    m_dielectric_yield);
    pp.query("work_function",  m_dielectric_work);
  }

  { // Noise parameters for initial dat
    ParmParse pp("init_data");
    pp.query("noise_amplitude", m_noise_amp);
    pp.query("noise_octaves",   m_noise_octaves);
    pp.query("noise_persist",   m_noise_persist);
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

Real morrow_lowke::initial_sigma(const RealVect& a_pos) const{
  return 0.;
}

morrow_lowke::electron::electron(){
  m_name      = "electron";
  m_charge    = -1;
  m_diffusive = true;

  m_uniform_density = 1.0;
  m_seed_density    = 0.0;
  m_seed_radius     = 1.0;
  m_noise_density   = 0.0;
  m_seed_pos        = RealVect::Zero;

  { // Get from input script or command line
    ParmParse pp("init_data");
    pp.query("uniform_density", m_uniform_density);
    pp.query("seed_density",    m_seed_density);
    pp.query("seed_radius",     m_seed_radius);
    pp.query("noise_density",     m_noise_density);
    if(pp.contains("seed_position")){
      Vector<Real> pos(SpaceDim);
      pp.queryarr("seed_position", pos, 0, SpaceDim);
      m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
    }
  }
}

morrow_lowke::electron::~electron(){
}

const Real morrow_lowke::electron::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/(m_seed_radius*m_seed_radius);
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

  m_uniform_density = 1.0;
  m_seed_density    = 0.0;
  m_seed_radius     = 1.0;
  m_noise_density   = 0.0;
  m_seed_pos        = RealVect::Zero;
  
  { // Get from input script or command line
    ParmParse pp("init_data");
    pp.query("uniform_density", m_uniform_density);
    pp.query("seed_density",    m_seed_density);
    pp.query("seed_radius",     m_seed_radius);
    pp.query("noise_density",     m_noise_density);
    if(pp.contains("seed_position")){
      Vector<Real> pos(SpaceDim);
      pp.queryarr("seed_position", pos, 0, SpaceDim);
      m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
    }
  }
}

morrow_lowke::positive_species::~positive_species(){
}

const Real morrow_lowke::positive_species::initial_data(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/(m_seed_radius*m_seed_radius);
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
}

morrow_lowke::negative_species::~negative_species(){
}

const Real morrow_lowke::negative_species::initial_data(const RealVect a_pos, const Real a_time) const {
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
    ParmParse pp("photon_one");
    pp.query("A_coeff",      m_A);
    pp.query("lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.8;
  Real pressure = 1.0;
  {
    ParmParse pp("gas");
    pp.query("frac_O2",  O2_frac);
    pp.query("pressure", pressure);
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
    ParmParse pp("photon_two");
    pp.query("A_coeff",      m_A);
    pp.query("lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.8;
  Real pressure = 1.0;
  {
    ParmParse pp("gas");
    pp.query("frac_O2",  O2_frac);
    pp.query("pressure", pressure);
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
    ParmParse pp("photon_three");
    pp.query("A_coeff",      m_A);
    pp.query("lambda_coeff", m_lambda);
  }

  // Find pressure. Need gas state for this. 
  Real O2_frac  = 0.8;
  Real pressure = 1.0;
  {
    ParmParse pp("gas");
    pp.query("frac_O2",  O2_frac);
    pp.query("pressure", pressure);
  }
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;  
}

morrow_lowke::photon_three::~photon_three(){
}

Real morrow_lowke::photon_three::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.

}
