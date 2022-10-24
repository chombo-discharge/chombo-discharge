/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaAir3LFA.cpp
  @brief  Implementation of CD_ItoPlasmaAir3LFA.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoPlasmaAir3LFA.H>
#include <CD_Units.H>
#include <CD_ParticleOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaAir3LFA::ItoPlasmaAir3LFA()
{
  m_num_ItoSpecies = 3;
  m_numRtSpecies   = 1;

  m_coupling == ItoPlasmaPhysics::coupling::LFA;

  ParmParse    pp("ItoPlasmaAir3LFA");
  Vector<Real> v;

  // For controlling time step
  pp.get("dX", m_deltaX);

  // Stuff for initial particles
  pp.get("seed", m_seed);
  pp.get("blob_radius", m_blob_radius);
  pp.get("num_particles", m_num_particles);
  pp.get("particle_weight", m_particle_weight);
  pp.getarr("blob_center", v, 0, SpaceDim);
  m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));

  // Mobile ions or not
  pp.get("mobile_ions", m_isMobile_ions);
  pp.get("ion_mobility", m_ion_mu);

  // Photostuff
  pp.get("quenching_pressure", m_pq);
  pp.get("photoi_factor", m_photoi_factor);

  // Algorithm stuff
  std::string str;
  pp.get("react_ppc", m_ppc);
  pp.get("poisson_switch", m_fieldSolver_switch);
  pp.get("Ncrit", m_Ncrit);
  pp.get("prop_eps", m_eps);
  pp.get("NSSA", m_NSSA);
  pp.get("SSAlim", m_SSAlim);

  // Get algorithm
  pp.get("algorithm", str);
  if (str == "hybrid") {
    m_algorithm = algorithm::hybrid;
  }
  else if (str == "tau") {
    m_algorithm = algorithm::tau;
  }
  else if (str == "ssa") {
    m_algorithm = algorithm::ssa;
  }
  else {
    MayDay::Abort("ItoPlasmaAir3LFA::ItoPlasmaAir3LFA - unknown algorithm requested");
  }

  // Random seed
  if (m_seed < 0) {
    m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }

  // Standard air.
  m_p      = 1.0;
  m_T      = 300;
  m_N2frac = 0.8;
  m_O2frac = 0.2;

  // Convert to SI units
  m_p  = m_p * Units::atm2pascal;
  m_pq = m_pq * Units::atm2pascal;
  m_N  = m_p * Units::Na / (m_T * Units::R);

  // Set up species
  m_ItoSpecies.resize(m_num_ItoSpecies);
  m_rtSpecies.resize(m_numRtSpecies);

  m_ElectronIdx = 0;
  m_PositiveIdx = 1;
  m_NegativeIdx = 2;
  m_PhotonZ_idx = 0;

  m_ItoSpecies[m_ElectronIdx] = RefCountedPtr<ItoSpecies>(new Electron());
  m_ItoSpecies[m_PositiveIdx] = RefCountedPtr<ItoSpecies>(new Positive());
  m_ItoSpecies[m_NegativeIdx] = RefCountedPtr<ItoSpecies>(new Negative());
  m_rtSpecies[m_PhotonZ_idx]  = RefCountedPtr<RtSpecies>(new PhotonZ());

  // To avoid that MPI ranks draw the same particle positions, increment the seed for each rank
  m_seed += procID();
  m_rng = std::mt19937_64(m_seed);

  List<ItoParticle>& Electrons = m_ItoSpecies[m_ElectronIdx]->getInitialParticles();
  List<ItoParticle>& Positives = m_ItoSpecies[m_PositiveIdx]->getInitialParticles();
  List<ItoParticle>& Negatives = m_ItoSpecies[m_NegativeIdx]->getInitialParticles();

  Electrons.clear();
  Positives.clear();
  Negatives.clear();

  // Draw some initial electrons and add corresponding positive ions. 
  ParticleOps::drawSphereParticles(Electrons, m_num_particles, m_blob_center, m_blob_radius);

  for (ListIterator<ItoParticle> lit(Electrons); lit.ok(); ++lit) {
    const RealVect& x = lit().position();

    lit().weight() = m_particle_weight;
    
    Positives.add(ItoParticle(m_particle_weight, x));    
  }

  // Particle-particle reactions
  m_reactions.emplace("impact_ionization",
                      ItoPlasmaReaction({m_ElectronIdx}, {m_ElectronIdx, m_ElectronIdx, m_PositiveIdx}));
  m_reactions.emplace("Electron_attachment", ItoPlasmaReaction({m_ElectronIdx}, {m_NegativeIdx}));
  m_reactions.emplace("Electron_recombination", ItoPlasmaReaction({m_ElectronIdx, m_PositiveIdx}, {}));
  m_reactions.emplace("ion_recombination", ItoPlasmaReaction({m_PositiveIdx, m_NegativeIdx}, {}));
  m_reactions.emplace("photo_excitation", ItoPlasmaReaction({m_ElectronIdx}, {m_ElectronIdx}, {m_PhotonZ_idx}));

  // Photo-reactions
  m_photoReactions.emplace("zheleznyak", ItoPlasmaPhotoReaction({m_PhotonZ_idx}, {m_ElectronIdx, m_PositiveIdx}));

  // Set the ions diffusion coefficient
  m_ion_D = m_ion_mu * Units::kb * m_T / Units::Qe;

  this->readTables();
}

ItoPlasmaAir3LFA::~ItoPlasmaAir3LFA() {}

void
ItoPlasmaAir3LFA::readTables()
{

  this->addTable("mobility", "mobility.dat");
  this->addTable("diffco", "diffusion.dat");
  this->addTable("alpha", "alpha.dat");
  this->addTable("eta", "eta.dat");
  this->addTable("Te", "energy.dat");

  // Need to scale tables. First column is in Td and second in /N
  m_tables["mobility"].scale<0>(m_N * Units::Td);
  m_tables["mobility"].scale<1>(1. / m_N);

  m_tables["diffco"].scale<0>(m_N * Units::Td);
  m_tables["diffco"].scale<1>(1. / m_N);

  m_tables["alpha"].scale<0>(m_N * Units::Td);
  m_tables["alpha"].scale<1>(m_N);

  m_tables["eta"].scale<0>(m_N * Units::Td);
  m_tables["eta"].scale<1>(m_N);

  m_tables["Te"].swap(0, 1);
  m_tables["Te"].scale<0>(m_N * Units::Td);
  m_tables["Te"].scale<1>(2.0 * Units::Qe / (3.0 * Units::kb));
}

Real
ItoPlasmaAir3LFA::computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_cdr_densities) const
{
  const Real alpha = this->computeAlpha(a_E);
  const Real velo  = this->computeElectronDriftVelocity(a_E).vectorLength();

  const Real k_alpha = alpha * velo;

  return log(m_deltaX) / k_alpha;
}

Real
ItoPlasmaAir3LFA::computeAlpha(const RealVect a_E) const
{
  const Real E = a_E.vectorLength();

  return m_tables.at("alpha").getEntry<1>(E);
}

Vector<Real>
ItoPlasmaAir3LFA::computeItoMobilitiesLFA(const Real a_time, const RealVect a_pos, const RealVect a_E) const
{
  Vector<Real> mobilities(m_num_ItoSpecies, m_ion_mu);
  mobilities[m_ElectronIdx] = m_tables.at("mobility").getEntry<1>(a_E.vectorLength());

  return mobilities;
}

RealVect
ItoPlasmaAir3LFA::computeElectronDriftVelocity(const RealVect a_E) const
{
  return -m_tables.at("mobility").getEntry<1>(a_E.vectorLength()) * a_E;
}

Vector<Real>
ItoPlasmaAir3LFA::computeItoDiffusionLFA(const Real         a_time,
                                         const RealVect     a_pos,
                                         const RealVect     a_E,
                                         const Vector<Real> a_cdr_densities) const
{
  Vector<Real> D(m_num_ItoSpecies, m_ion_D);
  D[m_ElectronIdx] = m_tables.at("diffco").getEntry<1>(a_E.vectorLength());

  return D;
}

void
ItoPlasmaAir3LFA::updateReactionRatesLFA(const RealVect a_E, const Real a_dx, const Real a_kappa) const
{

  // Compute the reaction rates.
  const Real dV      = pow(a_dx, SpaceDim); //*a_kappa;
  const Real E       = a_E.vectorLength();
  const Real alpha   = m_tables.at("alpha").getEntry<1>(E);
  const Real eta     = m_tables.at("eta").getEntry<1>(E);
  const Real Te      = m_tables.at("Te").getEntry<1>(E);
  const Real velo    = this->computeElectronDriftVelocity(a_E).vectorLength();
  const Real xfactor = (m_pq / (m_p + m_pq)) * excitationRates(E) * sergeyFactor(m_O2frac) * m_photoi_factor;
  const Real bpn     = 2E-13 * sqrt(300 / m_T) / dV;
  const Real bpe     = 1.138E-11 * pow(Te, -0.7) / dV;

  m_reactions.at("impact_ionization").rate()      = alpha * velo;
  m_reactions.at("Electron_attachment").rate()    = eta * velo;
  m_reactions.at("Electron_recombination").rate() = bpe;
  m_reactions.at("ion_recombination").rate()      = bpn;
  m_reactions.at("photo_excitation").rate()       = alpha * velo * xfactor;
}

Real
ItoPlasmaAir3LFA::excitationRates(const Real a_E) const
{
  const Real Etd = a_E / (m_N * Units::Td);

  Real y = 1.0;
  if (Etd > 100) {
    y = 0.1 * exp(233 / Etd);
  }

  return y;
}

Real
ItoPlasmaAir3LFA::sergeyFactor(const Real a_O2frac) const
{
  return 3E-2 + 0.4 * pow(a_O2frac, 0.6);
}

Real
ItoPlasmaAir3LFA::photoionizationRate(const Real a_E) const
{
  return excitationRates(a_E) * sergeyFactor(m_O2frac) * m_photoi_factor;
}

ItoPlasmaAir3LFA::Electron::Electron()
{
  m_isMobile     = true;
  m_isDiffusive  = true;
  m_name         = "Electron";
  m_chargeNumber = -1;
}

ItoPlasmaAir3LFA::Electron::~Electron() {}

ItoPlasmaAir3LFA::Positive::Positive()
{
  m_isMobile     = false;
  m_isDiffusive  = false;
  m_name         = "Positive";
  m_chargeNumber = 1;

  ParmParse pp("ItoPlasmaAir3LFA");
  pp.get("mobile_ions", m_isMobile);
  pp.get("mobile_ions", m_isDiffusive);
}

ItoPlasmaAir3LFA::Positive::~Positive() {}

ItoPlasmaAir3LFA::Negative::Negative()
{
  m_isMobile     = false;
  m_isDiffusive  = false;
  m_name         = "Negative";
  m_chargeNumber = -1;

  ParmParse pp("ItoPlasmaAir3LFA");
  pp.get("mobile_ions", m_isMobile);
  pp.get("mobile_ions", m_isDiffusive);
}

ItoPlasmaAir3LFA::Negative::~Negative() {}

ItoPlasmaAir3LFA::PhotonZ::PhotonZ()
{
  m_name = "PhotonZ";

  const Real O2_frac  = 0.2;
  const Real pressure = 1.0;

  ParmParse pp("ItoPlasmaAir3LFA");

  pp.get("photoi_f1", m_f1);
  pp.get("photoi_f2", m_f2);
  pp.get("photoi_K1", m_K1);
  pp.get("photoi_K2", m_K2);
  pp.get("photoi_seed", m_seed);

  // Convert units
  m_pO2 = pressure * O2_frac * Units::atm2pascal;
  m_K1  = m_K1 * m_pO2;
  m_K2  = m_K2 * m_pO2;

  // Seed the RNG
  if (m_seed < 0)
    m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng     = new std::mt19937_64(m_seed);
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
}

ItoPlasmaAir3LFA::PhotonZ::~PhotonZ() {}

Real
ItoPlasmaAir3LFA::PhotonZ::getAbsorptionCoefficient(const RealVect a_pos) const
{
  const Real f = m_f1 + (*m_udist01)(*m_rng) * (m_f2 - m_f1);
  return m_K1 * pow(m_K2 / m_K1, (f - m_f1) / (m_f2 - m_f1));
}

#include <CD_NamespaceFooter.H>
