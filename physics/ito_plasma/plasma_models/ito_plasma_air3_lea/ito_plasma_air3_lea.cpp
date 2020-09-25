/*!
  @file   ito_plasma_air3_lea.cpp
  @brief  Implementation of ito_plasma_air3_lea.H
  @author Robert Marskar
  @date   Aug. 2020
*/

#include "ito_plasma_air3_lea.H"
#include "units.H"

#include <ParmParse.H>

using namespace physics::ito_plasma;

ito_plasma_air3_lea::ito_plasma_air3_lea(){
  m_num_ito_species = 3;
  m_num_rte_species = 1;

  m_coupling = ito_plasma_physics::coupling::LEA;

  ParmParse pp("ito_plasma_air3_lea");
  Vector<Real> v;
  
  // Stuff for initial particles
  pp.get   ("seed",            m_seed);
  pp.get   ("blob_radius",     m_blob_radius);
  pp.get   ("num_particles",   m_num_particles);
  pp.get   ("particle_weight", m_particle_weight);
  pp.getarr("blob_center",     v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));

  // Photostuff
  pp.get("quenching_pressure", m_pq);
  pp.get("photoi_factor",      m_photoi_factor);

  // Algorithm stuff
  std::string str;
  pp.get("react_ppc",      m_ppc);
  pp.get("poisson_switch", m_poisson_switch);
  pp.get("Ncrit",          m_Ncrit);
  pp.get("prop_eps",       m_eps);
  pp.get("NSSA",           m_NSSA);
  pp.get("SSAlim",         m_SSAlim);
  
  // Get algorithm
  pp.get("algorithm",      str);
  if(str == "hybrid"){
    m_algorithm = algorithm::hybrid;
  }
  else if(str == "tau"){
    m_algorithm = algorithm::tau;
  }
  else if(str == "ssa"){
    m_algorithm = algorithm::ssa;
  }
  else{
    MayDay::Abort("ito_plasma_air3_lea::ito_plasma_air3_lea - unknown algorithm requested");
  }

  // Random seed
  if(m_seed < 0) {
    m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }

  // Standard air. 
  m_p = 1.0;
  m_T = 300;
  m_N2frac = 0.8;
  m_O2frac = 0.2;

  // Convert to SI units
  m_p  = m_p*units::s_atm2pascal;
  m_pq = m_pq*units::s_atm2pascal;
  m_N  = m_p*units::s_Na/(m_T*units::s_R);

  // Set up species
  m_ito_species.resize(m_num_ito_species);
  m_rte_species.resize(m_num_rte_species);

  m_electron_idx = 0;
  m_positive_idx = 1;
  m_negative_idx = 2;
  m_photonZ_idx  = 0;

  // Read transport tables
  this->read_tables();

  // Initiate species
  m_ito_species[m_electron_idx] = RefCountedPtr<ito_species> (new electron(m_tables.at("mobility"), m_tables.at("diffco")));
  m_ito_species[m_positive_idx] = RefCountedPtr<ito_species> (new positive());
  m_ito_species[m_negative_idx] = RefCountedPtr<ito_species> (new negative());
  m_rte_species[m_photonZ_idx]  = RefCountedPtr<rte_species> (new photonZ());

  // To avoid that MPI ranks draw the same particle positions, increment the seed for each rank
  m_seed += procID();
  m_rng   = std::mt19937_64(m_seed);

  List<ito_particle>& electrons = m_ito_species[m_electron_idx]->get_initial_particles();
  List<ito_particle>& positives = m_ito_species[m_positive_idx]->get_initial_particles();
  List<ito_particle>& negatives = m_ito_species[m_negative_idx]->get_initial_particles();

  electrons.clear();
  positives.clear();
  negatives.clear();
  
  this->draw_sphere_particles(electrons, positives, m_num_particles, m_blob_center, m_blob_radius, m_particle_weight, 2.0, 0.0);

  // Electron loss function
  std::pair<int, Real> impact_loss     = std::make_pair(m_electron_idx, -13.0);  //  14eV per reaction of this type.
  std::pair<int, Real> friction_loss   = std::make_pair(m_electron_idx, -1.0);   // -12eV per reaction of this type.
  std::pair<int, Real> photo_loss      = std::make_pair(m_electron_idx,  -15.0); //  15 eV per photoexcitation
  std::pair<int, Real> photo_gain      = std::make_pair(m_electron_idx,  0.0);   //  Energy of appearing photoelectrons

  // Particle-particle reactions
  m_reactions.emplace("impact_ionization",      ito_reaction({m_electron_idx}, {m_electron_idx, m_electron_idx, m_positive_idx}, {impact_loss}));
  // m_reactions.emplace("electron_attachment",    ito_reaction({m_electron_idx}, {m_negative_idx}));
  // m_reactions.emplace("electron_recombination", ito_reaction({m_electron_idx, m_positive_idx}, {}));
  // m_reactions.emplace("ion_recombination",      ito_reaction({m_positive_idx, m_negative_idx}, {}));
  m_reactions.emplace("photo_excitation",       ito_reaction({m_electron_idx}, {m_electron_idx}, {m_photonZ_idx}));
  m_reactions.emplace("electron_scattering",    ito_reaction({m_electron_idx}, {m_electron_idx}, {friction_loss}));

  // Photo-reactions
  m_photo_reactions.emplace("zheleznyak",  photo_reaction({m_photonZ_idx}, {m_electron_idx, m_positive_idx}, {photo_gain}));
}

ito_plasma_air3_lea::~ito_plasma_air3_lea(){

}

void ito_plasma_air3_lea::read_tables(){

  this->add_table("mobility",  "mobility.dat");
  this->add_table("diffco",    "diffusion.dat");
  this->add_table("alpha",     "alpha.dat");
  this->add_table("eta",       "eta.dat");
  this->add_table("alpha_lfa", "alpha_lfa.dat");
  this->add_table("collision", "collision_rate.dat");

  // Need to scale energy tables. First column is in eV
  m_tables["mobility"].scale_y(1./m_N);
  m_tables["diffco"].scale_y(1./m_N);
  m_tables["alpha"].scale_y(m_N);
  m_tables["eta"].scale_y(m_N);
  m_tables["collision"].scale_y(m_N);

  // LFA table for Townsed coefficient
  m_tables["alpha_lfa"].scale_x(m_N*units::s_Td);
  m_tables["alpha_lfa"].scale_y(m_N);
}

Real ito_plasma_air3_lea::compute_alpha(const RealVect a_E) const {
  const Real E = a_E.vectorLength();

  return m_tables.at("alpha_lfa").get_entry(E);
}

void ito_plasma_air3_lea::update_reaction_rates_lea(const RealVect a_E, const Vector<Real> a_mean_energies, const Real a_dx, const Real a_kappa) const {

  //  const Real electron_energy = a_mean_energies[m_electron_idx];
  //  const Real electron_energy = 12; // Use 12eV for everything for now.
  const Real electron_energy = Min(a_mean_energies[m_electron_idx], 30.0);

  // Compute the reaction rates.
  const Real dV      = pow(a_dx, SpaceDim);//*a_kappa;
  const Real E       = a_E.vectorLength();
  const Real alpha   = m_tables.at("alpha").get_entry(electron_energy);
  const Real eta     = m_tables.at("eta").get_entry(electron_energy);
  const Real mu      = m_tables.at("mobility").get_entry(electron_energy);
  const Real scat    = m_tables.at("collision").get_entry(electron_energy);
  const Real velo    = mu*a_E.vectorLength();
  const Real Te      = 2.0*a_mean_energies[m_electron_idx]/(3.0*units::s_kb);
  const Real xfactor = (m_pq/(m_p + m_pq))*excitation_rates(E)*sergey_factor(m_O2frac)*m_photoi_factor;
  const Real bpn     = 2E-13*sqrt(300/m_T)/dV;
  const Real bpe     = 1.138E-11*pow(Te, -0.7)/dV;
  
  m_reactions.at("impact_ionization").rate()      = alpha*velo;
  // m_reactions.at("electron_attachment").rate()    = eta*velo;
  // m_reactions.at("electron_recombination").rate() = bpe;
  // m_reactions.at("ion_recombination").rate()      = bpn;
  m_reactions.at("photo_excitation").rate()       = alpha*velo*xfactor;
  m_reactions.at("electron_scattering").rate()    = scat;
}

Real ito_plasma_air3_lea::excitation_rates(const Real a_E) const{
  const Real Etd = a_E/(m_N*units::s_Td);

  Real y = 1.0;
  if(Etd > 100){
    y = 0.1*exp(233/Etd);
  }

  return y;
}

Real ito_plasma_air3_lea::sergey_factor(const Real a_O2frac) const{
  return 3E-2 + 0.4*pow(a_O2frac, 0.6);
}

Real ito_plasma_air3_lea::photo_rate(const Real a_E) const {
  return excitation_rates(a_E)*sergey_factor(m_O2frac)*m_photoi_factor;
}

ito_plasma_air3_lea::electron::electron(const lookup_table& a_mobility, const lookup_table& a_diffusion) : m_mobility(a_mobility),
													   m_diffusion(a_diffusion) {
  m_mobile    = true;
  m_diffusive = true;
  m_name      = "electron";
  m_charge    = -1;
}

ito_plasma_air3_lea::electron::~electron(){

}

Real ito_plasma_air3_lea::electron::mobility(const Real a_energy) const {
  return m_mobility.get_entry(a_energy);
}

Real ito_plasma_air3_lea::electron::diffusion(const Real a_energy) const {
  return m_diffusion.get_entry(a_energy);
}

ito_plasma_air3_lea::positive::positive(){
  m_mobile    = false;
  m_diffusive = false;
  m_name      = "positive";
  m_charge    = 1;
}  

ito_plasma_air3_lea::positive::~positive(){

}

ito_plasma_air3_lea::negative::negative(){
  m_mobile    = false;
  m_diffusive = false;
  m_name      = "negative";
  m_charge    = -1;
}  

ito_plasma_air3_lea::negative::~negative(){
}

ito_plasma_air3_lea::photonZ::photonZ(){
  m_name   = "photonZ";

  const Real O2_frac  = 0.2;
  const Real pressure = 1.0;
  
  ParmParse pp("ito_plasma_air3_lea");
  
  pp.get("photoi_f1",   m_f1);
  pp.get("photoi_f2",   m_f2);
  pp.get("photoi_K1",   m_K1);
  pp.get("photoi_K2",   m_K2);
  pp.get("photoi_seed", m_seed);

  // Convert units
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
  m_K1  = m_K1*m_pO2;
  m_K2  = m_K2*m_pO2;

  // Seed the RNG
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
}

ito_plasma_air3_lea::photonZ::~photonZ(){

}

Real ito_plasma_air3_lea::photonZ::get_kappa(const RealVect a_pos) const {
  const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
  return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
}
