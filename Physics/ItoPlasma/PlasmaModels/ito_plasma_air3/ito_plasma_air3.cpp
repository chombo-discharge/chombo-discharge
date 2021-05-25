/*!
  @file   ito_plasma_air3.cpp
  @brief  Implementation of ito_plasma_air3.H
  @author Robert Marskar
  @date   Aug. 2020
*/

#include "ito_plasma_air3.H"
#include <CD_Units.H>

#include <ParmParse.H>

#include "CD_NamespaceHeader.H"
using namespace Physics::ItoPlasma;

ito_plasma_air3::ito_plasma_air3(){
  m_num_ItoSpecies = 3;
  m_numRtSpecies = 1;

  m_coupling == ItoPlasmaPhysics::coupling::LFA;

  ParmParse pp("ito_plasma_air3");
  Vector<Real> v;

  // For controlling time step
  pp.get("dX", m_deltaX);
  
  // Stuff for initial particles
  pp.get   ("seed",            m_seed);
  pp.get   ("blob_radius",     m_blob_radius);
  pp.get   ("num_particles",   m_num_particles);
  pp.get   ("particle_weight", m_particle_weight);
  pp.getarr("blob_center",     v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));

  // Mobile ions or not
  pp.get("mobile_ions",  m_isMobile_ions);
  pp.get("ion_mobility", m_ion_mu);

  // Photostuff
  pp.get("quenching_pressure", m_pq);
  pp.get("photoi_factor",      m_photoi_factor);

  // Algorithm stuff
  std::string str;
  pp.get("react_ppc",      m_ppc);
  pp.get("poisson_switch", m_fieldSolver_switch);
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
  m_p  = m_p*Units::atm2pascal;
  m_pq = m_pq*Units::atm2pascal;
  m_N  = m_p*Units::Na/(m_T*Units::R);

  // Set up species
  m_ItoSpecies.resize(m_num_ItoSpecies);
  m_RtSpecies.resize(m_numRtSpecies);

  m_electron_idx = 0;
  m_positive_idx = 1;
  m_negative_idx = 2;
  m_PhotonZ_idx  = 0;

  m_ItoSpecies[m_electron_idx] = RefCountedPtr<ItoSpecies> (new electron());
  m_ItoSpecies[m_positive_idx] = RefCountedPtr<ItoSpecies> (new positive());
  m_ItoSpecies[m_negative_idx] = RefCountedPtr<ItoSpecies> (new negative());
  m_RtSpecies[m_PhotonZ_idx]  = RefCountedPtr<RtSpecies> (new PhotonZ());

  // To avoid that MPI ranks draw the same particle positions, increment the seed for each rank
  m_seed += procID();
  m_rng   = std::mt19937_64(m_seed);

  List<ItoParticle>& electrons = m_ItoSpecies[m_electron_idx]->getInitialParticles();
  List<ItoParticle>& positives = m_ItoSpecies[m_positive_idx]->getInitialParticles();
  List<ItoParticle>& negatives = m_ItoSpecies[m_negative_idx]->getInitialParticles();

  electrons.clear();
  positives.clear();
  negatives.clear();
  
  this->drawSphereParticles(electrons, positives, m_num_particles, m_blob_center, m_blob_radius, m_particle_weight, 0.0, 0.0);

  // Particle-particle reactions
  m_reactions.emplace("impact_ionization",      ItoPlasmaReaction({m_electron_idx}, {m_electron_idx, m_electron_idx, m_positive_idx}));
  m_reactions.emplace("electron_attachment",    ItoPlasmaReaction({m_electron_idx}, {m_negative_idx}));
  m_reactions.emplace("electron_recombination", ItoPlasmaReaction({m_electron_idx, m_positive_idx}, {}));
  m_reactions.emplace("ion_recombination",      ItoPlasmaReaction({m_positive_idx, m_negative_idx}, {}));
  m_reactions.emplace("photo_excitation",       ItoPlasmaReaction({m_electron_idx}, {m_electron_idx}, {m_PhotonZ_idx}));

  // Photo-reactions
  m_ItoPlasmaPhotoReactions.emplace("zheleznyak",  ItoPlasmaPhotoReaction({m_PhotonZ_idx}, {m_electron_idx, m_positive_idx}));

  // Set the ions diffusion coefficient
  m_ion_D = m_ion_mu*Units::kb*m_T/Units::Qe;

  this->read_tables();
}

ito_plasma_air3::~ito_plasma_air3(){

}

void ito_plasma_air3::read_tables(){

  this->addTable("mobility", "mobility.dat");
  this->addTable("diffco",   "diffusion.dat");
  this->addTable("alpha",    "alpha.dat");
  this->addTable("eta",      "eta.dat");
  this->addTable("Te",       "energy.dat");

  // Need to scale tables. First column is in Td and second in /N
  m_tables["mobility"].scaleX(m_N*Units::Td);
  m_tables["mobility"].scaleY(1./m_N);

  m_tables["diffco"].scaleX(m_N*Units::Td);
  m_tables["diffco"].scaleY(1./m_N);

  m_tables["alpha"].scaleX(m_N*Units::Td);
  m_tables["alpha"].scaleY(m_N);

  m_tables["eta"].scaleX(m_N*Units::Td);
  m_tables["eta"].scaleY(m_N);

  m_tables["Te"].swapXY();
  m_tables["Te"].scaleX(m_N*Units::Td);
  m_tables["Te"].scaleY(2.0*Units::Qe/(3.0*Units::kb));
}

Real ito_plasma_air3::computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_cdr_densities) const {
  const Real alpha = this->computeAlpha(a_E);
  const Real velo  = this->compute_electron_velocity(a_E).vectorLength();

  const Real k_alpha = alpha*velo;

  return log(m_deltaX)/k_alpha;
}

Real ito_plasma_air3::computeAlpha(const RealVect a_E) const {
  const Real E = a_E.vectorLength();

  return m_tables.at("alpha").getEntry(E);
}

Vector<Real> ito_plasma_air3::computeItoMobilitiesLFA(const Real a_time, const RealVect a_pos, const RealVect a_E) const {
  Vector<Real> mobilities(m_num_ItoSpecies, m_ion_mu);
  mobilities[m_electron_idx] = m_tables.at("mobility").getEntry(a_E.vectorLength());

  return mobilities;
}

RealVect ito_plasma_air3::compute_electron_velocity(const RealVect a_E) const {
  return -m_tables.at("mobility").getEntry(a_E.vectorLength())*a_E;
}

Vector<Real> ito_plasma_air3::computeItoDiffusionLFA(const Real         a_time,
							const RealVect     a_pos,
							const RealVect     a_E,
							const Vector<Real> a_cdr_densities) const {
  Vector<Real> D(m_num_ItoSpecies, m_ion_D);
  D[m_electron_idx] = m_tables.at("diffco").getEntry(a_E.vectorLength());
  
  return D;
}

void ito_plasma_air3::updateReactionRatesLFA(const RealVect a_E, const Real a_dx, const Real a_kappa) const{

  // Compute the reaction rates.
  const Real dV      = pow(a_dx, SpaceDim);//*a_kappa;
  const Real E       = a_E.vectorLength();
  const Real alpha   = m_tables.at("alpha").getEntry(E);
  const Real eta     = m_tables.at("eta").getEntry(E);
  const Real Te      = m_tables.at("Te").getEntry(E);
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
  const Real Etd = a_E/(m_N*Units::Td);

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
  m_isMobile    = true;
  m_isDiffusive = true;
  m_name      = "electron";
  m_chargeNumber    = -1;
}

ito_plasma_air3::electron::~electron(){

}

ito_plasma_air3::positive::positive(){
  m_isMobile    = false;
  m_isDiffusive = false;
  m_name      = "positive";
  m_chargeNumber    = 1;

  ParmParse pp("ito_plasma_air3");
  pp.get("mobile_ions", m_isMobile);
  pp.get("mobile_ions", m_isDiffusive);
}  

ito_plasma_air3::positive::~positive(){

}

ito_plasma_air3::negative::negative(){
  m_isMobile    = false;
  m_isDiffusive = false;
  m_name      = "negative";
  m_chargeNumber    = -1;

  
  ParmParse pp("ito_plasma_air3");
  pp.get("mobile_ions", m_isMobile);
  pp.get("mobile_ions", m_isDiffusive);
}  

ito_plasma_air3::negative::~negative(){
}

ito_plasma_air3::PhotonZ::PhotonZ(){
  m_name   = "PhotonZ";

  const Real O2_frac  = 0.2;
  const Real pressure = 1.0;
  
  ParmParse pp("ito_plasma_air3");
  
  pp.get("photoi_f1",   m_f1);
  pp.get("photoi_f2",   m_f2);
  pp.get("photoi_K1",   m_K1);
  pp.get("photoi_K2",   m_K2);
  pp.get("photoi_seed", m_seed);

  // Convert units
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
  m_K1  = m_K1*m_pO2;
  m_K2  = m_K2*m_pO2;

  // Seed the RNG
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
}

ito_plasma_air3::PhotonZ::~PhotonZ(){

}

Real ito_plasma_air3::PhotonZ::getKappa(const RealVect a_pos) const {
  const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
  return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
}
#include "CD_NamespaceFooter.H"
