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
#include <CD_DataParser.H>
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaAir3LFA::ItoPlasmaAir3LFA()
{
  m_className        = "ItoPlasmaAir3LFA";
  m_numPlasmaSpecies = 3;
  m_numRtSpecies     = 1;

  // Default parameter for lookup tables
  m_table_entries = 1000;

  ParmParse    pp("ItoPlasmaAir3LFA");
  Vector<Real> v;

  // For controlling time step
  pp.get("dX", m_deltaX);

  // Stuff for initial particles
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

  this->parsePPC();
  this->parseAlgorithm();
  this->parseDebug();

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
  m_plasmaSpecies.resize(m_numPlasmaSpecies);
  m_rtSpecies.resize(m_numRtSpecies);

  m_ElectronIdx = 0;
  m_PositiveIdx = 1;
  m_NegativeIdx = 2;
  m_PhotonZ_idx = 0;

  m_plasmaSpecies[m_ElectronIdx] = RefCountedPtr<ItoSpecies>(new Electron());
  m_plasmaSpecies[m_PositiveIdx] = RefCountedPtr<ItoSpecies>(new Positive());
  m_plasmaSpecies[m_NegativeIdx] = RefCountedPtr<ItoSpecies>(new Negative());
  m_rtSpecies[m_PhotonZ_idx]     = RefCountedPtr<RtSpecies>(new PhotonZ());

  List<ItoParticle>& Electrons = m_plasmaSpecies[m_ElectronIdx]->getInitialParticles();
  List<ItoParticle>& Positives = m_plasmaSpecies[m_PositiveIdx]->getInitialParticles();
  List<ItoParticle>& Negatives = m_plasmaSpecies[m_NegativeIdx]->getInitialParticles();

  Electrons.clear();
  Positives.clear();
  Negatives.clear();

  // Draw some initial electrons and add corresponding positive ions.
  ParticleManagement::drawSphereParticles(Electrons, m_num_particles, m_blob_center, m_blob_radius);

  for (ListIterator<ItoParticle> lit(Electrons); lit.ok(); ++lit) {
    const RealVect& x = lit().position();

    lit().weight() = m_particle_weight;

    Positives.add(ItoParticle(m_particle_weight, x));
  }

  auto r1 = std::make_shared<KMCReaction>(std::list<size_t>{0}, std::list<size_t>{0, 0, 1}, std::list<size_t>{});
  auto r2 = std::make_shared<KMCReaction>(std::list<size_t>{0}, std::list<size_t>{0, 2}, std::list<size_t>{});
  auto r3 = std::make_shared<KMCReaction>(std::list<size_t>{0, 1}, std::list<size_t>{}, std::list<size_t>{});
  auto r4 = std::make_shared<KMCReaction>(std::list<size_t>{1, 2}, std::list<size_t>{}, std::list<size_t>{});
  auto r5 = std::make_shared<KMCReaction>(std::list<size_t>{1}, std::list<size_t>{1}, std::list<size_t>{0});

  m_kmcReactions.emplace_back(r1);
  m_kmcReactions.emplace_back(r2);
  m_kmcReactions.emplace_back(r3);
  m_kmcReactions.emplace_back(r4);
  m_kmcReactions.emplace_back(r5);

  // Photo-reactions
  m_photoReactions.emplace("zheleznyak", ItoPlasmaPhotoReaction({m_PhotonZ_idx}, {m_ElectronIdx, m_PositiveIdx}));

  // Set the ions diffusion coefficient
  m_ion_D = m_ion_mu * Units::kb * m_T / Units::Qe;

  this->defineKMC();

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
ItoPlasmaAir3LFA::computeItoMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept
{
  Vector<Real> mobilities(m_numPlasmaSpecies, m_ion_mu);
  mobilities[m_ElectronIdx] = m_tables.at("mobility").getEntry<1>(a_E.vectorLength());

  return mobilities;
}

RealVect
ItoPlasmaAir3LFA::computeElectronDriftVelocity(const RealVect a_E) const
{
  return -m_tables.at("mobility").getEntry<1>(a_E.vectorLength()) * a_E;
}

inline void
ItoPlasmaAir3LFA::addTable(const std::string a_table_name, const std::string a_file)
{

  LookupTable1D<2> table = DataParser::simpleFileReadASCII(a_file);

  table.sort();
  table.makeUniform(m_table_entries);    // Make table into a unifom table
  m_tables.emplace(a_table_name, table); // Add table
}

Vector<Real>
ItoPlasmaAir3LFA::computeItoDiffusion(const Real         a_time,
                                      const RealVect     a_pos,
                                      const RealVect     a_E,
                                      const Vector<Real> a_cdr_densities) const noexcept
{
  Vector<Real> D(m_numPlasmaSpecies, m_ion_D);
  D[m_ElectronIdx] = m_tables.at("diffco").getEntry<1>(a_E.vectorLength());

  return D;
}

void
ItoPlasmaAir3LFA::updateReactionRates(const RealVect a_E, const Real a_dx, const Real a_kappa) const noexcept
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

  m_kmcReactions[0]->rate() = alpha * velo;
  m_kmcReactions[1]->rate() = eta * velo;
  m_kmcReactions[2]->rate() = bpe;
  m_kmcReactions[3]->rate() = bpn;
  m_kmcReactions[4]->rate() = alpha * velo * xfactor;
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

  // Convert units
  m_pO2 = pressure * O2_frac * Units::atm2pascal;
  m_K1  = m_K1 * m_pO2;
  m_K2  = m_K2 * m_pO2;
}

ItoPlasmaAir3LFA::PhotonZ::~PhotonZ() {}

Real
ItoPlasmaAir3LFA::PhotonZ::getAbsorptionCoefficient(const RealVect a_pos) const
{
  const Real f = m_f1 + Random::getUniformReal01() * (m_f2 - m_f1);
  return m_K1 * pow(m_K2 / m_K1, (f - m_f1) / (m_f2 - m_f1));
}

#include <CD_NamespaceFooter.H>
