/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaAir3LEA.cpp
  @brief  Implementation of CD_ItoPlasmaAir3LEA.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoPlasmaAir3LEA.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaAir3LEA::ItoPlasmaAir3LEA()
{
  m_num_ItoSpecies = 3;
  m_numRtSpecies   = 1;

  m_coupling = ItoPlasmaPhysics::coupling::LEA;

  ParmParse    pp("ItoPlasmaAir3LEA");
  Vector<Real> v;

  // Stuff for initial particles
  pp.get("seed", m_seed);
  pp.get("blob_radius", m_blob_radius);
  pp.get("num_particles", m_num_particles);
  pp.get("particle_weight", m_particle_weight);
  pp.getarr("blob_center", v, 0, SpaceDim);
  m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));

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
    MayDay::Abort("ItoPlasmaAir3LEA::ItoPlasmaAir3LEA - unknown algorithm requested");
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

  m_electronIdx = 0;
  m_positiveIdx = 1;
  m_negativeIdx = 2;
  m_PhotonZ_idx = 0;

  // Read transport tables
  this->readTables();

  // Initiate species
  m_ItoSpecies[m_electronIdx] = RefCountedPtr<ItoSpecies>(new Electron(m_tables.at("mobility"), m_tables.at("diffco")));
  m_ItoSpecies[m_positiveIdx] = RefCountedPtr<ItoSpecies>(new Positive());
  m_ItoSpecies[m_negativeIdx] = RefCountedPtr<ItoSpecies>(new Negative());
  m_rtSpecies[m_PhotonZ_idx]  = RefCountedPtr<RtSpecies>(new PhotonZ());

  // To avoid that MPI ranks draw the same particle positions, increment the seed for each rank
  m_seed += procID();
  m_rng = std::mt19937_64(m_seed);

  List<ItoParticle>& Electrons = m_ItoSpecies[m_electronIdx]->getInitialParticles();
  List<ItoParticle>& Positives = m_ItoSpecies[m_positiveIdx]->getInitialParticles();
  List<ItoParticle>& Negatives = m_ItoSpecies[m_negativeIdx]->getInitialParticles();

  Electrons.clear();
  Positives.clear();
  Negatives.clear();

  this->drawSphereParticles(Electrons,
                            Positives,
                            m_num_particles,
                            m_blob_center,
                            m_blob_radius,
                            m_particle_weight,
                            2.0,
                            0.0);

  // Electron loss function
  std::pair<int, Real> impact_loss   = std::make_pair(m_electronIdx, -13.0); //  14eV per reaction of this type.
  std::pair<int, Real> friction_loss = std::make_pair(m_electronIdx, -1.0);  // -1eV per "friction" collision.
  std::pair<int, Real> photo_loss    = std::make_pair(m_electronIdx, -15.0); //  15 eV per photoexcitation
  std::pair<int, Real> photo_gain    = std::make_pair(m_electronIdx, 2.0);   //  Energy of appearing photoElectrons

  // Particle-particle reactions
  m_reactions.emplace("impact_ionization",
                      ItoPlasmaReaction({m_electronIdx}, {m_electronIdx, m_electronIdx, m_positiveIdx}, {impact_loss}));
  // m_reactions.emplace("Electron_attachment",    ItoPlasmaReaction({m_electronIdx}, {m_negativeIdx}));
  // m_reactions.emplace("Electron_recombination", ItoPlasmaReaction({m_electronIdx, m_positiveIdx}, {}));
  // m_reactions.emplace("ion_recombination",      ItoPlasmaReaction({m_positiveIdx, m_negativeIdx}, {}));
  // m_reactions.emplace("photo_excitation",       ItoPlasmaReaction({m_electronIdx}, {m_electronIdx}, {m_PhotonZ_idx}, {photo_loss}));
  m_reactions.emplace("Electron_scattering", ItoPlasmaReaction({m_electronIdx}, {m_electronIdx}, {friction_loss}));

  // Photo-reactions
  //  m_photoReactions.emplace("zheleznyak",  ItoPlasmaPhotoReaction({m_PhotonZ_idx}, {m_electronIdx, m_positiveIdx}, {photo_gain}));
}

ItoPlasmaAir3LEA::~ItoPlasmaAir3LEA() {}

void
ItoPlasmaAir3LEA::readTables()
{

  this->addTable("mobility", "mobility.dat");
  this->addTable("diffco", "diffusion.dat");
  this->addTable("alpha", "alpha.dat");
  this->addTable("eta", "eta.dat");
  this->addTable("alpha_lfa", "alpha_lfa.dat");
  this->addTable("collision", "collision_rate.dat");

  // Need to scale energy tables. First column is in eV
  m_tables["mobility"].scale<1>(1. / m_N);
  m_tables["diffco"].scale<1>(1. / m_N);
  m_tables["alpha"].scale<1>(m_N);
  m_tables["eta"].scale<1>(m_N);
  m_tables["collision"].scale<1>(m_N);

  // LFA table for Townsed coefficient
  m_tables["alpha_lfa"].scale<0>(m_N * Units::Td);
  m_tables["alpha_lfa"].scale<1>(m_N);
}

Real
ItoPlasmaAir3LEA::computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_cdr_densities) const
{
  return 1.E99;
}

Real
ItoPlasmaAir3LEA::computeAlpha(const RealVect a_E) const
{
  const Real E = a_E.vectorLength();

  return m_tables.at("alpha_lfa").getEntry<1>(E);
}

void
ItoPlasmaAir3LEA::updateReactionRatesLEA(const RealVect     a_E,
                                         const Vector<Real> a_mean_energies,
                                         const Real         a_dx,
                                         const Real         a_kappa) const
{

  //  const Real Electron_energy = a_mean_energies[m_electronIdx];
  //  const Real Electron_energy = 12; // Use 12eV for everything for now.
  const Real Electron_energy = Min(a_mean_energies[m_electronIdx], 30.0);

  // Compute the reaction rates.
  const Real dV      = pow(a_dx, SpaceDim); //*a_kappa;
  const Real E       = a_E.vectorLength();
  const Real alpha   = m_tables.at("alpha").getEntry<1>(Electron_energy);
  const Real eta     = m_tables.at("eta").getEntry<1>(Electron_energy);
  const Real mu      = m_tables.at("mobility").getEntry<1>(Electron_energy);
  const Real scat    = m_tables.at("collision").getEntry<1>(Electron_energy);
  const Real velo    = mu * a_E.vectorLength();
  const Real Te      = 2.0 * a_mean_energies[m_electronIdx] / (3.0 * Units::kb);
  const Real xfactor = (m_pq / (m_p + m_pq)) * excitationRates(E) * sergeyFactor(m_O2frac) * m_photoi_factor;
  const Real bpn     = 2E-13 * sqrt(300 / m_T) / dV;
  const Real bpe     = 1.138E-11 * pow(Te, -0.7) / dV;

  m_reactions.at("impact_ionization").rate() = alpha * velo;
  // m_reactions.at("Electron_attachment").rate()    = eta*velo;
  // m_reactions.at("Electron_recombination").rate() = bpe;
  // m_reactions.at("ion_recombination").rate()      = bpn;
  // m_reactions.at("photo_excitation").rate()       = alpha*velo*xfactor;
  m_reactions.at("Electron_scattering").rate() = scat;
}

Real
ItoPlasmaAir3LEA::excitationRates(const Real a_E) const
{
  const Real Etd = a_E / (m_N * Units::Td);

  Real y = 1.0;
  if (Etd > 100) {
    y = 0.1 * exp(233 / Etd);
  }

  return y;
}

Real
ItoPlasmaAir3LEA::sergeyFactor(const Real a_O2frac) const
{
  return 3E-2 + 0.4 * pow(a_O2frac, 0.6);
}

Real
ItoPlasmaAir3LEA::photoionizationRate(const Real a_E) const
{
  return excitationRates(a_E) * sergeyFactor(m_O2frac) * m_photoi_factor;
}

ItoPlasmaAir3LEA::Electron::Electron(const LookupTable1D<2>& a_mobility, const LookupTable1D<2>& a_diffusion)
  : m_mobility(a_mobility), m_diffusion(a_diffusion)
{
  m_isMobile     = true;
  m_isDiffusive  = true;
  m_name         = "Electron";
  m_chargeNumber = -1;
}

ItoPlasmaAir3LEA::Electron::~Electron() {}

Real
ItoPlasmaAir3LEA::Electron::mobility(const Real a_energy) const
{
  const Real energy = Max(a_energy, 30.0);
  return m_mobility.getEntry<1>(energy);
}

Real
ItoPlasmaAir3LEA::Electron::diffusion(const Real a_energy) const
{
  const Real energy = Max(a_energy, 30.0);
  return m_diffusion.getEntry<1>(energy);
}

ItoPlasmaAir3LEA::Positive::Positive()
{
  m_isMobile     = false;
  m_isDiffusive  = false;
  m_name         = "Positive";
  m_chargeNumber = 1;
}

ItoPlasmaAir3LEA::Positive::~Positive() {}

ItoPlasmaAir3LEA::Negative::Negative()
{
  m_isMobile     = false;
  m_isDiffusive  = false;
  m_name         = "Negative";
  m_chargeNumber = -1;
}

ItoPlasmaAir3LEA::Negative::~Negative() {}

ItoPlasmaAir3LEA::PhotonZ::PhotonZ()
{
  m_name = "PhotonZ";

  const Real O2_frac  = 0.2;
  const Real pressure = 1.0;

  ParmParse pp("ItoPlasmaAir3LEA");

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

ItoPlasmaAir3LEA::PhotonZ::~PhotonZ() {}

Real
ItoPlasmaAir3LEA::PhotonZ::getAbsorptionCoefficient(const RealVect a_pos) const
{
  const Real f = m_f1 + (*m_udist01)(*m_rng) * (m_f2 - m_f1);
  return m_K1 * pow(m_K2 / m_K1, (f - m_f1) / (m_f2 - m_f1));
}

#include <CD_NamespaceFooter.H>
