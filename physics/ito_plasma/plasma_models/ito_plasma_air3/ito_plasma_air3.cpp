/*!
  @file   ito_plasma_air3.cpp
  @brief  Implementation of ito_plasma_air3.H
  @author Robert Marskar
  @date   Aug. 2020
*/

#include "ito_plasma_air3.H"
#include "units.H"

#include <ParmParse.H>

using namespace physics::ito_plasma;

ito_plasma_air3::ito_plasma_air3(){
  m_num_ito_species = 3;
  m_num_rte_species = 1;

  ParmParse pp("ito_plasma_air3");
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
    MayDay::Abort("ito_plasma_air3::ito_plasma_air3 - unknown algorithm requested");
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

  m_ito_species[m_electron_idx] = RefCountedPtr<ito_species> (new electron());
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
  
  this->draw_sphere_particles(electrons, positives, m_num_particles, m_blob_center, m_blob_radius, m_particle_weight);


  // Particle-particle reactions
  m_reactions.emplace("impact_ionization",      ito_reaction({m_electron_idx}, {m_electron_idx, m_electron_idx, m_positive_idx}, 
  m_reactions.emplace("electron_attachment",    ito_reaction({m_electron_idx}, {m_negative_idx}));
  m_reactions.emplace("electron_recombination", ito_reaction({m_electron_idx, m_positive_idx}, {}));
  m_reactions.emplace("ion_recombination",      ito_reaction({m_positive_idx, m_negative_idx}, {}));
  m_reactions.emplace("photo_excitation",       ito_reaction({m_electron_idx}, {m_electron_idx}, {m_photonZ_idx}));

  // Photo-reactions
  m_photo_reactions.emplace("zheleznyak",  photo_reaction({m_photonZ_idx}, {m_electron_idx, m_positive_idx}));

  // Read reaction tables. This also makes them "uniform".
  this->read_tables();
}

ito_plasma_air3::~ito_plasma_air3(){

}

void ito_plasma_air3::read_tables(){

  this->add_table("mobility", "mobility.dat");
  this->add_table("diffco",   "diffusion.dat");
  this->add_table("alpha",    "alpha.dat");
  this->add_table("eta",      "eta.dat");
  this->add_table("Te",       "energy.dat");

  // Need to scale tables. First column is in Td and second in /N
  m_tables["mobility"].scale_x(m_N*units::s_Td);
  m_tables["mobility"].scale_y(1./m_N);

  m_tables["diffco"].scale_x(m_N*units::s_Td);
  m_tables["diffco"].scale_y(1./m_N);

  m_tables["alpha"].scale_x(m_N*units::s_Td);
  m_tables["alpha"].scale_y(m_N);

  m_tables["eta"].scale_x(m_N*units::s_Td);
  m_tables["eta"].scale_y(m_N);

  m_tables["Te"].swap_xy();
  m_tables["Te"].scale_x(m_N*units::s_Td);
  m_tables["Te"].scale_y(2.0*units::s_Qe/(3.0*units::s_kb));
}

Real ito_plasma_air3::compute_alpha(const RealVect a_E) const {
  const Real E = a_E.vectorLength();

  return m_tables.at("alpha").get_entry(E);
}

Vector<Real> ito_plasma_air3::compute_ito_mobilities_lfa(const Real a_time, const RealVect a_pos, const RealVect a_E) const {
  Vector<Real> mobilities(m_num_ito_species, 0.0);
  mobilities[m_electron_idx] = m_tables.at("mobility").get_entry(a_E.vectorLength());

  return mobilities;
}

RealVect ito_plasma_air3::compute_electron_velocity(const RealVect a_E) const {
  return -m_tables.at("mobility").get_entry(a_E.vectorLength())*a_E;
}

Vector<Real> ito_plasma_air3::compute_ito_diffusion_lfa(const Real         a_time,
						    const RealVect     a_pos,
						    const RealVect     a_E,
						    const Vector<Real> a_cdr_densities) const {
  Vector<Real> D(m_num_ito_species, 0.0);
  D[m_electron_idx] = m_tables.at("diffco").get_entry(a_E.vectorLength());
  
  return D;
}

void ito_plasma_air3::update_reaction_rates_lfa(const RealVect a_E, const Real a_dx, const Real a_kappa) const{

  // Compute the reaction rates.
  const Real dV      = pow(a_dx, SpaceDim);//*a_kappa;
  const Real E       = a_E.vectorLength();
  const Real alpha   = m_tables.at("alpha").get_entry(E);
  const Real eta     = m_tables.at("eta").get_entry(E);
  const Real Te      = m_tables.at("Te").get_entry(E);
  const Real velo    = this->compute_electron_velocity(a_E).vectorLength();
  const Real xfactor = (m_pq/(m_p + m_pq))*excitation_rates(E)*sergey_factor(m_O2frac)*m_photoi_factor;
  const Real bpn     = 2E-13*sqrt(300/m_T)/dV;
  const Real bpe     = 1.138E-11*pow(Te, -0.7)/dV;
  
  m_reactions.at("impact_ionization").rate()      = alpha*velo;
  m_reactions.at("electron_attachment").rate()    = eta*velo;
  m_reactions.at("electron_recombination").rate() = bpe;
  m_reactions.at("ion_recombination").rate()      = bpn;
  m_reactions.at("photo_excitation").rate()       = alpha*velo*xfactor;
}

Real ito_plasma_air3::excitation_rates(const Real a_E) const{
  const Real Etd = a_E/(m_N*units::s_Td);

  Real y = 1.0;
  if(Etd > 100){
    y = 0.1*exp(233/Etd);
  }

  return y;
}

Real ito_plasma_air3::sergey_factor(const Real a_O2frac) const{
  return 3E-2 + 0.4*pow(a_O2frac, 0.6);
}

Real ito_plasma_air3::photo_rate(const Real a_E) const {
  return excitation_rates(a_E)*sergey_factor(m_O2frac)*m_photoi_factor;
}

ito_plasma_air3::electron::electron(){
  m_mobile    = true;
  m_diffusive = true;
  m_name      = "electron";
  m_charge    = -1;
}

ito_plasma_air3::electron::~electron(){

}

ito_plasma_air3::positive::positive(){
  m_mobile    = false;
  m_diffusive = false;
  m_name      = "positive";
  m_charge    = 1;
}  

ito_plasma_air3::positive::~positive(){

}

ito_plasma_air3::negative::negative(){
  m_mobile    = false;
  m_diffusive = false;
  m_name      = "negative";
  m_charge    = -1;
}  

ito_plasma_air3::negative::~negative(){
}  

ito_plasma_air3::photonZ::photonZ(){
  m_name   = "photonZ";

  const Real O2_frac  = 0.2;
  const Real pressure = 1.0;
  
  ParmParse pp("ito_plasma_air3");
  
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

ito_plasma_air3::photonZ::~photonZ(){

}

Real ito_plasma_air3::photonZ::get_kappa(const RealVect a_pos) const {
  const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
  return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
}
