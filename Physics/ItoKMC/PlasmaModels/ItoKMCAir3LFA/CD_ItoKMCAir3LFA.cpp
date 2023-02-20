/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCAir3LFA.cpp
  @brief  Implementation of CD_ItoKMCAir3LFA.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoKMCAir3LFA.H>
#include <CD_Units.H>
#include <CD_DataParser.H>
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

ItoKMCAir3LFA::ItoKMCAir3LFA()
{
  CH_TIME("ItoKMCAir3LFA::ItoKMCAir3LFA");

  m_className = "ItoKMCAir3LFA";

  // Parse some transport settings
  this->parsePPC();
  this->parseAlgorithm();
  this->parseDebug();
  this->parseDx();
  this->parseTransport();

  // Set gas parameters.
  m_p         = 1.0 * Units::atm2pascal;
  m_T         = 300.0;
  m_N2frac    = 0.8;
  m_O2frac    = 0.2;
  m_pq        = m_pq * Units::atm2pascal;
  m_N         = m_p * Units::Na / (m_T * Units::R);
  m_ionDiffCo = m_ionMobility * Units::kb * m_T / Units::Qe;

  // Set up species
  m_itoSpecies.resize(1);
  m_cdrSpecies.resize(2);
  m_rtSpecies.resize(1);

  m_itoSpecies[0] = RefCountedPtr<ItoSpecies>(new Electron());
  m_cdrSpecies[0] = RefCountedPtr<CdrSpecies>(new Positive());
  m_cdrSpecies[1] = RefCountedPtr<CdrSpecies>(new Negative());
  m_rtSpecies[0]  = RefCountedPtr<RtSpecies>(new PhotonZ());

  List<ItoParticle>& electrons = m_itoSpecies[0]->getInitialParticles();

  electrons.clear();

  // Draw some initial electrons
  ParmParse pp("ItoKMCAir3LFA");

  int          initParticles;
  Real         initParticleWeight;
  Real         blobRadius;
  RealVect     blobCenter;
  Vector<Real> v;

  pp.get("init_particles", initParticles);
  pp.get("init_weights", initParticleWeight);
  pp.get("init_radius", blobRadius);
  pp.getarr("init_center", v, 0, SpaceDim);
  blobCenter = RealVect(D_DECL(v[0], v[1], v[2]));

  ParticleManagement::drawSphereParticles(electrons, initParticles, blobCenter, blobRadius);

  // Add reactions. These are
  //
  // e -> e + e + M+ (ionization)
  // e -> M- (attachment)
  // e + M+ -> null (electron-ion recombination)
  // M- + M+ -> null (ion-ion recombination)
  // e -> e + Y (photon generation)
  //
  auto r1 = std::make_shared<KMCReaction>(std::list<size_t>{0}, std::list<size_t>{0, 0, 1}, std::list<size_t>{});
  auto r2 = std::make_shared<KMCReaction>(std::list<size_t>{0}, std::list<size_t>{2}, std::list<size_t>{});
  auto r3 = std::make_shared<KMCReaction>(std::list<size_t>{0, 1}, std::list<size_t>{}, std::list<size_t>{});
  auto r4 = std::make_shared<KMCReaction>(std::list<size_t>{1, 2}, std::list<size_t>{}, std::list<size_t>{});
  auto r5 = std::make_shared<KMCReaction>(std::list<size_t>{0}, std::list<size_t>{0}, std::list<size_t>{0});

  m_kmcReactions.emplace_back(r1);
  m_kmcReactions.emplace_back(r2);
  m_kmcReactions.emplace_back(r3);
  m_kmcReactions.emplace_back(r4);
  m_kmcReactions.emplace_back(r5);

  // Photo-reactions
  auto y1 = std::make_shared<ItoKMCPhotoReaction>(0, std::list<size_t>{0, 1});
  m_photoReactions.emplace_back(y1);

  // Run the rest of the define method.
  this->define();

  // Read input data
  this->readTables();
}

ItoKMCAir3LFA::~ItoKMCAir3LFA() { CH_TIME("ItoKMCAir3LFA::~ItoKMCAir3LFA"); }

void
ItoKMCAir3LFA::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoKMCAir3LFA::parseRuntimeOptions");

  this->parsePPC();
  this->parseDebug();
  this->parseAlgorithm();
  this->parseDx();
}

void
ItoKMCAir3LFA::parseDx() noexcept
{
  CH_TIME("ItoKMCAir3LFA::parseDx");

  ParmParse pp(m_className.c_str());

  pp.get("dX", m_deltaX);
}

void
ItoKMCAir3LFA::parseTransport() noexcept
{
  CH_TIME("ItoKMCAir3LFA::parseTransport");

  ParmParse pp(m_className.c_str());

  // Tunable photoionization parameters.
  pp.get("quenching_pressure", m_pq);
  pp.get("photoi_factor", m_photoIonizationFactor);
  pp.get("ion_mobility", m_ionMobility);
  pp.get("impact_efficiency", m_ionImpactEfficiency);
  pp.get("quantum_efficiency", m_quantumEfficiency);
  pp.get("extrap_bc", m_extrapBC);

  m_ionImpactEfficiency = std::max(m_ionImpactEfficiency, 0.0);
  m_ionImpactEfficiency = std::min(m_ionImpactEfficiency, 1.0);

  m_quantumEfficiency = std::max(m_quantumEfficiency, 0.0);
  m_quantumEfficiency = std::min(m_quantumEfficiency, 1.0);
}

void
ItoKMCAir3LFA::readTables() noexcept
{
  CH_TIME("ItoKMCAir3LFA::readTables");

  // Read tables in
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

void
ItoKMCAir3LFA::addTable(const std::string a_table_name, const std::string a_file) noexcept
{
  CH_TIME("ItoKMCAir3LFA::addTable");

  LookupTable1D<2> table = DataParser::simpleFileReadASCII(a_file);

  table.sort();
  table.makeUniform(1000);
  m_tables.emplace(a_table_name, table);
}

Real
ItoKMCAir3LFA::computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_numParticles) const noexcept
{
  CH_TIME("ItoKMCAir3LFA::computeDt");
  const Real E     = a_E.vectorLength();
  const Real alpha = this->computeAlpha(a_E);
  const Real mu    = m_tables.at("mobility").getEntry<1>(E);

  return log(m_deltaX) / (mu * E * alpha);
}

Real
ItoKMCAir3LFA::computeAlpha(const RealVect a_E) const
{
  CH_TIME("ItoKMCAir3LFA::computeAlpha");

  return m_tables.at("alpha").getEntry<1>(a_E.vectorLength());
}

Vector<Real>
ItoKMCAir3LFA::computeMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept
{
  CH_TIME("ItoKMCAir3LFA::computeMobilities");

  Vector<Real> mobilities(3);

  mobilities[0] = m_tables.at("mobility").getEntry<1>(a_E.vectorLength());
  mobilities[1] = m_ionMobility;
  mobilities[2] = m_ionMobility;

  return mobilities;
}

Vector<Real>
ItoKMCAir3LFA::computeDiffusionCoefficients(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept
{
  Vector<Real> D(3);

  D[0] = m_tables.at("diffco").getEntry<1>(a_E.vectorLength());
  D[1] = m_ionDiffCo;
  D[2] = m_ionDiffCo;

  return D;
}

void
ItoKMCAir3LFA::secondaryEmissionEB(Vector<List<ItoParticle>>&       a_secondaryParticles,
                                   Vector<Real>&                    a_cdrFluxes,
                                   Vector<List<Photon>>&            a_secondaryPhotons,
                                   const Vector<List<ItoParticle>>& a_primaryParticles,
                                   const Vector<Real>&              a_extrapCDRFluxes,
                                   const Vector<List<Photon>>&      a_primaryPhotons,
                                   const RealVect&                  a_E,
                                   const RealVect&                  a_cellCenter,
                                   const RealVect&                  a_cellCentroid,
                                   const RealVect&                  a_bndryCentroid,
                                   const RealVect&                  a_bndryNormal,
                                   const Real                       a_bndryArea,
                                   const Real                       a_dx,
                                   const Real                       a_dt,
                                   const bool                       a_isDielectric,
                                   const int                        a_matIndex) const noexcept
{
  CH_TIME("ItoKMCAir3LFA::injectParticlesEB");
#if 0
  List<ItoParticle>&       ingoingElectrons = a_secondaryParticles[0];
  const List<ItoParticle>& outgoingIons     = a_primaryParticles[1];
  const List<Photon>&      outgoingPhotons  = a_primaryPhotons[0];

  if (m_extrapBC) {
    const Real diff   = a_oldNumParticles[0] - a_newNumParticles[0];
    const Real diffLL = (long long)diff;
    if (diffLL > 0LL) {
      ingoingElectrons.add(ItoParticle(diff, a_cellCenter + a_cellCentroid * a_dx));
    }
  }
  else {
    const bool isCathode = a_E.dotProduct(a_bndryNormal) < 0.0;
    const bool isAnode   = a_E.dotProduct(a_bndryNormal) > 0.0;

    const RealVect bndryPos = a_cellCenter + a_dx * a_bndryCentroid;

    if (isCathode) {
      // SEE from ion impact
      for (ListIterator<ItoParticle> lit(outgoingIons); lit.ok(); ++lit) {
        const long long success = Random::getBinomial<long long>(llround(lit().weight()), m_ionImpactEfficiency);
        if (success > 0LL) {
          ingoingElectrons.add(ItoParticle(1.0 * success, bndryPos));
        }
      }

      // SEE from photo-electric effect
      for (ListIterator<Photon> lit(outgoingPhotons); lit.ok(); ++lit) {
        const long long success = Random::getBinomial<long long>(llround(lit().weight()), m_quantumEfficiency);
        if (success > 0LL) {
          ingoingElectrons.add(ItoParticle(1.0 * success, bndryPos));
        }
      }
    }
    else if (isAnode) {
    }
  }
#endif
}

void
ItoKMCAir3LFA::updateReactionRates(const RealVect a_E, const Real a_dx, const Real a_kappa) const noexcept
{
  CH_TIME("ItoKMCAir3LFA::updateReactionRates");

  // Compute the reaction rates.
  const Real dV      = pow(a_dx, SpaceDim);
  const Real E       = a_E.vectorLength();
  const Real alpha   = m_tables.at("alpha").getEntry<1>(E);
  const Real eta     = m_tables.at("eta").getEntry<1>(E);
  const Real Te      = m_tables.at("Te").getEntry<1>(E);
  const Real mu      = m_tables.at("mobility").getEntry<1>(E);
  const Real velo    = mu * E;
  const Real xfactor = (m_pq / (m_p + m_pq)) * excitationRates(E) * sergeyFactor(m_O2frac) * m_photoIonizationFactor;
  const Real bpn     = 2E-13 * sqrt(300 / m_T) / dV;
  const Real bpe     = 1.138E-11 * pow(Te, -0.7) / dV;

  m_kmcReactions[0]->rate() = alpha * velo;
  m_kmcReactions[1]->rate() = eta * velo;
  m_kmcReactions[2]->rate() = bpe;
  m_kmcReactions[3]->rate() = bpn;
  m_kmcReactions[4]->rate() = alpha * velo * xfactor;
}

Real
ItoKMCAir3LFA::excitationRates(const Real a_E) const
{
  CH_TIME("ItoKMCAir3LFA::excitationRates");

  const Real Etd = a_E / (m_N * Units::Td);

  Real y = 1.0;
  if (Etd > 100.0) {
    y = 0.1 * exp(233 / Etd);
  }

  return y;
}

Real
ItoKMCAir3LFA::sergeyFactor(const Real a_O2frac) const
{
  CH_TIME("ItoKMCAir3LFA::sergeyFactor");

  return 3E-2 + 0.4 * pow(a_O2frac, 0.6);
}

ItoKMCAir3LFA::Electron::Electron()
{
  CH_TIME("ItoKMCAir3LFA::Electron::Electron");

  m_isMobile     = true;
  m_isDiffusive  = true;
  m_name         = "Electron";
  m_chargeNumber = -1;
}

ItoKMCAir3LFA::Electron::~Electron() { CH_TIME("ItoKMCAir3LFA::Electron::Electron"); }

ItoKMCAir3LFA::Positive::Positive()
{
  CH_TIME("ItoKMCAir3LFA::Positive::Positive");

  m_isMobile     = false;
  m_isDiffusive  = false;
  m_name         = "Positive";
  m_chargeNumber = 1;

  ParmParse pp("ItoKMCAir3LFA");
  pp.get("ion_transport", m_isMobile);
  pp.get("ion_transport", m_isDiffusive);
}

ItoKMCAir3LFA::Positive::~Positive() { CH_TIME("ItoKMCAir3LFA::Positive::~Positive"); }

Real
ItoKMCAir3LFA::Positive::initialData(const RealVect a_pos, const Real a_time) const
{
  return 0.0;
}

ItoKMCAir3LFA::Negative::Negative()
{
  CH_TIME("ItoKMCAir3LFA::Negative::Negative");

  m_isMobile     = false;
  m_isDiffusive  = false;
  m_name         = "Negative";
  m_chargeNumber = -1;

  ParmParse pp("ItoKMCAir3LFA");
  pp.get("ion_transport", m_isMobile);
  pp.get("ion_transport", m_isDiffusive);
}

ItoKMCAir3LFA::Negative::~Negative() { CH_TIME("ItoKMCAir3LFA::Negative::~Negative"); }

Real
ItoKMCAir3LFA::Negative::initialData(const RealVect a_pos, const Real a_time) const
{
  return 0.0;
}

ItoKMCAir3LFA::PhotonZ::PhotonZ()
{
  CH_TIME("ItoKMCAir3LFA::PhotonZ::PhotonZ");

  m_name = "PhotonZ";

  ParmParse pp("ItoKMCAir3LFA");

  pp.get("photoi_f1", m_f1);
  pp.get("photoi_f2", m_f2);
  pp.get("photoi_K1", m_K1);
  pp.get("photoi_K2", m_K2);

  constexpr Real fracO2   = 0.2;
  constexpr Real pressure = 1.0;
  constexpr Real pO2      = fracO2 * pressure * Units::atm2pascal;

  m_K1 = m_K1 * pO2;
  m_K2 = m_K2 * pO2;
}

ItoKMCAir3LFA::PhotonZ::~PhotonZ() { CH_TIME("ItoKMCAir3LFA::PhotonZ::~PhotonZ"); }

Real
ItoKMCAir3LFA::PhotonZ::getAbsorptionCoefficient(const RealVect a_pos) const
{
  // Draw random frequency in interval.
  const Real f = m_f1 + Random::getUniformReal01() * (m_f2 - m_f1);

  // Return random absorption coefficient.
  return m_K1 * pow(m_K2 / m_K1, (f - m_f1) / (m_f2 - m_f1));
}

#include <CD_NamespaceFooter.H>
