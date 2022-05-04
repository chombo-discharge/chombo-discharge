/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir9EedBourdon.cpp
  @brief  Implementation of CD_CdrPlasmaAir9EedBourdon.H
  @author Robert Marskar
*/

// Std includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Chombo includes
#include <PolyGeom.H>
#include <ParmParse.H>

// Our includes
#include <CD_CdrPlasmaAir9EedBourdon.H>
#include <CD_CdrPlasmaAir9EedBourdonSpecies.H>
#include <CD_DataParser.H>
#include <CD_Units.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

std::string CdrPlasmaAir9EedBourdon::s_bolsig_energy_E = "Energy (eV) 	Electric field / N (Td)";
std::string CdrPlasmaAir9EedBourdon::s_bolsig_mobility = "Energy (eV)	Mobility *N (1/m/V/s)";
std::string CdrPlasmaAir9EedBourdon::s_bolsig_N2_alpha = "C25   N2    Ionization    15.60 eV";
std::string CdrPlasmaAir9EedBourdon::s_bolsig_O2_alpha = "C42   O2    Ionization    12.06 eV";
std::string CdrPlasmaAir9EedBourdon::s_bolsig_townsend = "Energy (eV)	Townsend ioniz. coef. alpha/N (m2)";

CdrPlasmaAir9EedBourdon::CdrPlasmaAir9EedBourdon()
{
  m_numCdrSpecies = 9; // 8 reactive ones plus the eed
  m_numRtSpecies  = 3; // Bourdon model for Photons
  //  m_eed_solve   = true; // Yes, we're doing an EED solve so we must have a Poisson solution first
  //  m_eed_index   = 0;    // Index for the EED equation

  parseTransportFile();
  parseGas_parameters(m_Tg, m_p, m_N, m_O2frac, m_N2frac); // Get gas parameters
  parsePhotoi();
  parseSEE();
  parseTransport();

  initSpecies();

  // Read stuff from BOLSIG+ output files and make the tables uniform and with the correct units
  read_Electron_mobility();
  read_initEed();
  read_e_N2_alpha();
  read_e_O2_alpha();
  read_townsend();
}

CdrPlasmaAir9EedBourdon::~CdrPlasmaAir9EedBourdon() {}

void
CdrPlasmaAir9EedBourdon::parseTransportFile()
{
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  pp.get("transport_file", m_transport_file);
  pp.get("uniform_tables", m_uniform_tables);

  std::ifstream infile(m_transport_file);
  if (!infile.good()) {
    MayDay::Abort("CdrPlasmaAir9EedBourdon::parseTransportFile - could not find transport data");
  }
  else {
    infile.close();
  }
}

void
CdrPlasmaAir9EedBourdon::parseGas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac)
{
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  pp.get("gas_pressure", a_p);

  // This has a hard definition from BOLSIG+
  a_Tg     = 300.;
  a_O2frac = 0.20;
  a_N2frac = 0.80;

  const Real tot_frac = a_O2frac + a_N2frac;
  a_p                 = a_p * Units::atm2pascal;
  a_O2frac            = a_O2frac / tot_frac; // Normalize to one
  a_N2frac            = a_N2frac / tot_frac;
  a_N                 = a_p * Units::Na / (a_Tg * Units::R);
}

void
CdrPlasmaAir9EedBourdon::parsePhotoi()
{
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  pp.get("photoionization_efficiency", m_photoionization_efficiency);
  pp.get("excitation_efficiency", m_excitation_efficiency);
  pp.get("quenching_pressure", m_pq);

  m_pq *= Units::atm2pascal;
}

void
CdrPlasmaAir9EedBourdon::parseSEE()
{
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  pp.get("electrode_townsend2", m_townsend2_electrode);
  pp.get("dielectric_townsend2", m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency", m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
  pp.get("cathode_work", m_cathode_work);
  pp.get("dielectric_work", m_dielectric_work);
}

void
CdrPlasmaAir9EedBourdon::parseTransport()
{
  ParmParse pp("CdrPlasmaAir9EedBourdon");

  std::string str;

  pp.get("diffusive_electrons", str);
  m_isDiffusive_Electrons = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);
  m_isDiffusive_ions = (str == "true") ? true : false;
  pp.get("mobile_ions", str);
  m_isMobile_ions = (str == "true") ? true : false;

  pp.get("ion_mobility", m_ion_mobility);

  m_ion_diffusion = m_ion_mobility * (Units::kb * m_Tg) / Units::Qe;
}

void
CdrPlasmaAir9EedBourdon::initSpecies()
{
  m_cdrSpecies.resize(m_numCdrSpecies);
  m_eed_idx      = 0;
  m_Electron_idx = 1;
  m_N2plus_idx   = 2;
  m_N4plus_idx   = 3;
  m_O2plus_idx   = 4;
  m_O4plus_idx   = 5;
  m_O2plusN2_idx = 6;
  m_O2minus_idx  = 7;
  m_Ominus_idx   = 8;

  m_cdrSpecies[m_eed_idx]      = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::eed());
  m_cdrSpecies[m_Electron_idx] = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::Electron());
  m_cdrSpecies[m_N2plus_idx]   = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::N2plus());
  m_cdrSpecies[m_N4plus_idx]   = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::N4plus());
  m_cdrSpecies[m_O2plus_idx]   = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::O2plus());
  m_cdrSpecies[m_O4plus_idx]   = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::O4plus());
  m_cdrSpecies[m_O2plusN2_idx] = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::O2plusN2());
  m_cdrSpecies[m_O2minus_idx]  = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::O2minus());
  m_cdrSpecies[m_Ominus_idx]   = RefCountedPtr<CdrSpecies>(new CdrPlasmaAir9EedBourdon::Ominus());

  // Instantiate Photon solvers
  m_rtSpecies.resize(m_numRtSpecies);
  m_Photon1_idx = 0;
  m_Photon2_idx = 1;
  m_Photon3_idx = 2;

  m_rtSpecies[m_Photon1_idx] = RefCountedPtr<RtSpecies>(new CdrPlasmaAir9EedBourdon::PhotonOne());
  m_rtSpecies[m_Photon2_idx] = RefCountedPtr<RtSpecies>(new CdrPlasmaAir9EedBourdon::PhotonTwo());
  m_rtSpecies[m_Photon3_idx] = RefCountedPtr<RtSpecies>(new CdrPlasmaAir9EedBourdon::PhotonThree());
}

void
CdrPlasmaAir9EedBourdon::read_e_N2_alpha()
{

  // Read file entries
  m_e_N2_alpha = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir9EedBourdon::s_bolsig_N2_alpha, "");
  m_e_N2_alpha.sort();
  m_e_N2_alpha.makeUniform(m_uniform_tables);
}

void
CdrPlasmaAir9EedBourdon::read_e_O2_alpha()
{

  // Read file entries
  m_e_O2_alpha = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir9EedBourdon::s_bolsig_O2_alpha, "");
  m_e_O2_alpha.sort();
  m_e_O2_alpha.makeUniform(m_uniform_tables);
}

void
CdrPlasmaAir9EedBourdon::read_Electron_mobility()
{

  // Read file entries
  m_e_mobility = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir9EedBourdon::s_bolsig_mobility, "");
  m_e_mobility.sort();

  // Scale with density and make a uniform table (there's no guarantee that BOLSIG output is uniform!)
  m_e_mobility.scale<1>(1. / m_N);
  m_e_mobility.makeUniform(m_uniform_tables);
}

void
CdrPlasmaAir9EedBourdon::read_townsend()
{

  // Read file entries
  m_alpha_townsend =
    DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir9EedBourdon::s_bolsig_townsend, "");
  m_alpha_townsend.sort();

  // Scale with density and make a uniform table (there's no guarantee that BOLSIG output is uniform!)
  m_alpha_townsend.makeUniform(m_uniform_tables);
  m_alpha_townsend.scale<1>(m_N);

  //  m_alpha_townsend.dumpTable();
}

void
CdrPlasmaAir9EedBourdon::read_initEed()
{
  m_initEed = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir9EedBourdon::s_bolsig_energy_E, "");

  // Input table is in reverse order. Then make it uniform.
  m_initEed.swap(0, 1);
  m_initEed.sort();
  m_initEed.makeUniform(m_uniform_tables);
}

Real
CdrPlasmaAir9EedBourdon::computeAlpha(const Real a_E, const RealVect a_pos) const
{
  const Real EbyN   = a_E / (m_N * Units::Td);
  const Real energy = m_initEed.getEntry<1>(EbyN);
  const Real alpha  = m_alpha_townsend.getEntry<1>(energy);

  return alpha;
}

Real
CdrPlasmaAir9EedBourdon::compute_Electron_energy(const Real a_energy_density, const Real a_Electron_density) const
{
  const Real factor = 1.E8;
  return a_energy_density / (factor + a_Electron_density);
}

Vector<Real>
CdrPlasmaAir9EedBourdon::computeCdrDiffusionCoefficients(const Real         a_time,
                                                         const RealVect     a_pos,
                                                         const RealVect     a_E,
                                                         const Vector<Real> a_cdr_densities) const
{

  Vector<Real> diffco(m_numCdrSpecies, 0.0);

  const Real Electron_energy = compute_Electron_energy(a_cdr_densities[m_eed_idx], a_cdr_densities[m_Electron_idx]);
  const Real EbyN            = (a_E / (m_N * Units::Td)).vectorLength();

  diffco[m_eed_idx]      = this->compute_eed_diffco(Electron_energy);
  diffco[m_Electron_idx] = this->compute_e_diffco(Electron_energy);
  diffco[m_N2plus_idx]   = this->compute_N2plus_diffco();
  diffco[m_N4plus_idx]   = this->compute_N4plus_diffco();
  diffco[m_O2plus_idx]   = this->compute_O2plus_diffco();
  diffco[m_O4plus_idx]   = this->compute_O4plus_diffco();
  diffco[m_O2plusN2_idx] = this->compute_O2plusN2_diffco();
  diffco[m_O2minus_idx]  = this->compute_O2minus_diffco();

  return diffco;
}

Vector<RealVect>
CdrPlasmaAir9EedBourdon::computeCdrDriftVelocities(const Real         a_time,
                                                   const RealVect     a_pos,
                                                   const RealVect     a_E,
                                                   const Vector<Real> a_cdr_densities) const
{
  Vector<RealVect> velocities(m_numCdrSpecies, RealVect::Zero);

  const Real Electron_energy = compute_Electron_energy(a_cdr_densities[m_eed_idx], a_cdr_densities[m_Electron_idx]);
  const Real EbyN            = (a_E / (m_N * Units::Td)).vectorLength();

  velocities[m_eed_idx]      = this->compute_eed_mobility(Electron_energy) * (-a_E);
  velocities[m_Electron_idx] = this->compute_e_mobility(Electron_energy) * (-a_E);
  if (m_isMobile_ions) {
    velocities[m_N2plus_idx]   = this->compute_N2plus_mobility(EbyN) * a_E;
    velocities[m_N4plus_idx]   = this->compute_N4plus_mobility(EbyN) * a_E;
    velocities[m_O2plus_idx]   = this->compute_O2plus_mobility(EbyN) * a_E;
    velocities[m_O4plus_idx]   = this->compute_O4plus_mobility(EbyN) * a_E;
    velocities[m_O2plusN2_idx] = this->compute_O2plusN2_mobility(EbyN) * a_E;
    velocities[m_O2minus_idx]  = this->compute_O2minus_mobility(EbyN) * (-a_E);
    velocities[m_Ominus_idx]   = this->compute_Ominus_mobility(EbyN) * (-a_E);
  }

  return velocities;
}

void
CdrPlasmaAir9EedBourdon::advanceReactionNetwork(Vector<Real>&          a_cdr_sources,
                                                Vector<Real>&          a_rte_sources,
                                                const Vector<Real>     a_cdr_densities,
                                                const Vector<RealVect> a_cdr_gradients,
                                                const Vector<Real>     a_rte_densities,
                                                const RealVect         a_E,
                                                const RealVect         a_pos,
                                                const Real             a_dx,
                                                const Real             a_dt,
                                                const Real             a_time,
                                                const Real             a_kappa) const
{

  for (int i = 0; i < a_cdr_sources.size(); i++) {
    a_cdr_sources[i] = 0.0;
  }

  const Real Electron_energy = compute_Electron_energy(a_cdr_densities[m_eed_idx], a_cdr_densities[m_Electron_idx]);
  const Real Te              = Max(300., 2.0 * (Electron_energy * Units::Qe) / (3.0 * Units::kb)); // Kelvin
  const Real EbyN            = (a_E / (m_N * Units::Td)).vectorLength();

#if 0 // Debug
  if(Electron_energy > 1.E6) std::cout << Electron_energy << std::endl;
#endif

  // Room for improvement: The best thing would be to store the rate coefficients as matrices and then do S = K*n

  // Get all rate constant
  const Real k1  = this->compute_Electron_N2_alpha(Electron_energy);
  const Real k2  = this->compute_Electron_O2_alpha(Electron_energy);
  const Real k3  = this->compute_N2plus_N2_M_to_N4plus_M();
  const Real k4  = this->compute_N4plus_O2_to_O2_2N2();
  const Real k5  = this->compute_N2plus_O2_to_O2plus_N2(m_Tg);
  const Real k6  = this->compute_O2plus_2N2_to_O2plusN2_N2(m_Tg);
  const Real k7  = this->compute_O2plusN2_N2_to_O2plus_2N2(m_Tg);
  const Real k8  = this->compute_O2plusN2_O2_to_O4plus_N2();
  const Real k9  = this->compute_O2plus_O2_M_to_O4plus_M(m_Tg);
  const Real k10 = this->compute_e_O4plus_to_2O2(Te);
  const Real k11 = this->compute_e_O2plus_to_O2(Te);
  const Real k12 = this->compute_e_2O2_to_O2minus_O2(Te);
  const Real k13 = this->compute_O2minus_O4plus_to_3O2();
  const Real k14 = this->compute_O2minus_O4plus_M_to_3O2_M(m_Tg);
  const Real k15 = this->compute_O2minus_O2plus_M_to_2O2_M(m_Tg);
  const Real k16 = this->compute_e_O2_to_e_2O_c1(Electron_energy);
  const Real k17 = this->compute_e_O2_to_e_2O_c2(Electron_energy);
  const Real k18 = this->compute_e_O2_to_Ominus_O(Electron_energy);
  const Real k19 = this->compute_Oplus_O2_to_O_O2(m_Tg);
  const Real k20 = this->compute_e_N2_to_e_N2(Electron_energy);
  const Real k21 = this->compute_e_O2_to_e_O2(Electron_energy);

  // Electron energy losses for Electron collisions
  const Real dE_k1  = this->compute_e_N2_ionization_loss();
  const Real dE_k2  = this->compute_e_O2_ionization_loss();
  const Real dE_k16 = this->compute_e_O2_dissociation_loss_c1();
  const Real dE_k17 = this->compute_e_O2_dissociation_loss_c2();
  const Real dE_k18 = this->compute_e_O2_dissociative_attachment_loss();
  const Real dE_k20 = this->compute_e_N2_scattering_loss();
  const Real dE_k21 = this->compute_e_O2_scattering_loss();

  const Real n_e     = a_cdr_densities[m_Electron_idx];
  const Real n_N2    = m_N * m_N2frac;
  const Real n_O2    = m_N * m_O2frac;
  const Real n_N2p   = a_cdr_densities[m_N2plus_idx];
  const Real n_N4p   = a_cdr_densities[m_N4plus_idx];
  const Real n_O2p   = a_cdr_densities[m_O2plus_idx];
  const Real n_O4p   = a_cdr_densities[m_O4plus_idx];
  const Real n_O2pN2 = a_cdr_densities[m_O2plusN2_idx];
  const Real n_O2m   = a_cdr_densities[m_O2minus_idx];
  const Real n_Om    = a_cdr_densities[m_Ominus_idx];

  Real loss;
  Real products;

  // Compute Electron velocity, both drift and diffusion
  const RealVect ve = this->compute_e_mobility(Electron_energy) * (-a_E);
  const Real     De = this->compute_e_diffco(Electron_energy);
  const RealVect je = ve * n_e - De * a_cdr_gradients[m_Electron_idx];

  // Joule heating
  loss = PolyGeom::dot(-je, a_E);
  a_cdr_sources[m_eed_idx] += loss;

  // k1 reaction
  loss     = dE_k1;
  products = k1 * n_e * n_N2;
  a_cdr_sources[m_eed_idx] -= products * loss;
  a_cdr_sources[m_Electron_idx] += products;
  a_cdr_sources[m_N2plus_idx] += products;

  // k2 reaction
  loss     = dE_k2;
  products = k2 * n_e * n_O2;
  a_cdr_sources[m_eed_idx] -= products * loss;
  a_cdr_sources[m_Electron_idx] += products;
  a_cdr_sources[m_O2plus_idx] += products;

  // k3 reaction.
  products = k3 * n_N2p * n_N2 * m_N;
  a_cdr_sources[m_N2plus_idx] -= products;
  a_cdr_sources[m_N4plus_idx] += products;

  // k4 reaction
  products = k4 * n_N4p * n_O2;
  a_cdr_sources[m_N4plus_idx] -= products;
  a_cdr_sources[m_O2plus_idx] += products;

  // k5 reaction
  products = k5 * n_N2p * n_O2;
  a_cdr_sources[m_N2plus_idx] -= products;
  a_cdr_sources[m_O2plus_idx] += products;

  // k6 reaction
  products = k6 * n_O2p * n_N2 * n_N2;
  a_cdr_sources[m_O2plus_idx] -= products;
  a_cdr_sources[m_O2plusN2_idx] += products;

  // k7 reaction
  products = k7 * n_O2pN2 * n_N2;
  a_cdr_sources[m_O2plusN2_idx] -= products;
  a_cdr_sources[m_O2plus_idx] += products;

  // k8 reaction
  products = k8 * n_O2pN2 * n_O2;
  a_cdr_sources[m_O2plusN2_idx] -= products;
  a_cdr_sources[m_O4plus_idx] += products;

  // k9 reaction
  products = k9 * n_O2p * n_O2 * m_N;
  a_cdr_sources[m_O2plus_idx] -= products;
  a_cdr_sources[m_O4plus_idx] += products;

  // k10 reaction
  products = k10 * n_e * n_O4p;
  a_cdr_sources[m_Electron_idx] -= products;
  a_cdr_sources[m_O4plus_idx] -= products;

  // k11 reaction
  products = k11 * n_e * n_O2p;
  a_cdr_sources[m_Electron_idx] -= products;
  a_cdr_sources[m_O2plus_idx] -= products;

  // k12 reaction
  products = k12 * n_e * n_O2 * n_O2;
  a_cdr_sources[m_Electron_idx] -= products;
  a_cdr_sources[m_O2minus_idx] += products;

  // k13 reaction
  products = k13 * n_O2m * n_O4p;
  a_cdr_sources[m_O2minus_idx] -= products;
  a_cdr_sources[m_O4plus_idx] -= products;

  // k14 reaction
  products = k14 * n_O2m * n_O4p * m_N;
  a_cdr_sources[m_O2minus_idx] -= products;
  a_cdr_sources[m_O4plus_idx] -= products;

  // k15 reaction
  products = k15 * n_O2m * n_O2p * m_N;
  a_cdr_sources[m_O2minus_idx] -= products;
  a_cdr_sources[m_O2plus_idx] -= products;

  // This is here because the adapted functions
#if 1
  // k16 reaction
  loss     = dE_k16;
  products = k16 * n_e * n_O2;
  a_cdr_sources[m_eed_idx] -= products * loss;

  // k17 reaction
  loss     = dE_k17;
  products = k17 * n_e * n_O2;
  a_cdr_sources[m_eed_idx] -= products * loss;

  // k18 reaction
  loss     = dE_k18;
  products = k18 * n_e * n_O2;
  a_cdr_sources[m_eed_idx] -= products * loss;
  a_cdr_sources[m_Electron_idx] -= products;
  a_cdr_sources[m_Ominus_idx] += products;

  // k19 reaction
  products = k19 * n_Om * n_O2p;
  a_cdr_sources[m_Ominus_idx] -= products;
  a_cdr_sources[m_O2plus_idx] -= products;

  // k20 reaction
  loss     = dE_k20;
  products = k20 * n_e * n_N2;
  a_cdr_sources[m_eed_idx] -= products * loss;

  // k21 reaction
  loss     = dE_k21;
  products = k21 * n_e * n_O2;
  a_cdr_sources[m_eed_idx] -= products * loss;
#endif

  // Photoionization gamma + O2 -> e + O2+
  const CdrPlasmaAir9EedBourdon::PhotonOne* Photon1 =
    static_cast<CdrPlasmaAir9EedBourdon::PhotonOne*>(&(*m_rtSpecies[m_Photon1_idx]));
  const CdrPlasmaAir9EedBourdon::PhotonTwo* Photon2 =
    static_cast<CdrPlasmaAir9EedBourdon::PhotonTwo*>(&(*m_rtSpecies[m_Photon2_idx]));
  const CdrPlasmaAir9EedBourdon::PhotonThree* Photon3 =
    static_cast<CdrPlasmaAir9EedBourdon::PhotonThree*>(&(*m_rtSpecies[m_Photon3_idx]));
  products = m_photoionization_efficiency * Units::c * m_O2frac * m_p *
             (Photon1->getA() * a_rte_densities[m_Photon1_idx] + Photon2->getA() * a_rte_densities[m_Photon2_idx] +
              Photon3->getA() * a_rte_densities[m_Photon3_idx]);

  a_cdr_sources[m_Electron_idx] += products;
  a_cdr_sources[m_O2plus_idx] += products;

  // Photon source terms
  const Real kphot = k1 * a_cdr_densities[m_Electron_idx] * m_N * m_N2frac;

  a_rte_sources[m_Photon1_idx] = kphot * m_excitation_efficiency * (m_pq / (m_pq + m_p));
  a_rte_sources[m_Photon2_idx] = kphot * m_excitation_efficiency * (m_pq / (m_pq + m_p));
  a_rte_sources[m_Photon3_idx] = kphot * m_excitation_efficiency * (m_pq / (m_pq + m_p));

#if 0 // Debug
  Real sum = 0.0;
  sum += -a_cdr_sources[1] + a_cdr_sources[2] + a_cdr_sources[3] + a_cdr_sources[4] + a_cdr_sources[5] + a_cdr_sources[6] - a_cdr_sources[7] - a_cdr_sources[8];
  std::cout << sum << std::endl;
#endif
}

Vector<Real>
CdrPlasmaAir9EedBourdon::computeCdrDomainFluxes(const Real           a_time,
                                                const RealVect       a_pos,
                                                const int            a_dir,
                                                const Side::LoHiSide a_side,
                                                const RealVect       a_E,
                                                const Vector<Real>   a_cdr_densities,
                                                const Vector<Real>   a_cdr_velocities,
                                                const Vector<Real>   a_cdr_gradients,
                                                const Vector<Real>   a_rte_fluxes,
                                                const Vector<Real>   a_extrap_cdr_fluxes) const
{
  Vector<Real> fluxes(m_numCdrSpecies, 0.0);

  return a_extrap_cdr_fluxes;

  const int sgn = sign(a_side);
  for (int i = 0; i < fluxes.size(); i++) {
    fluxes[i] = sgn * Max(sgn * a_extrap_cdr_fluxes[i], 0.); // Outflow
  }
  return fluxes;
}

Vector<Real>
CdrPlasmaAir9EedBourdon::computeCdrDielectricFluxes(const Real         a_time,
                                                    const RealVect     a_pos,
                                                    const RealVect     a_normal,
                                                    const RealVect     a_E,
                                                    const Vector<Real> a_cdr_densities,
                                                    const Vector<Real> a_cdr_velocities,
                                                    const Vector<Real> a_cdr_gradients,
                                                    const Vector<Real> a_rte_fluxes,
                                                    const Vector<Real> a_extrap_cdr_fluxes) const
{

  return this->computeCdrFluxes(a_time,
                                a_pos,
                                a_normal,
                                a_E,
                                a_cdr_densities,
                                a_cdr_velocities,
                                a_cdr_gradients,
                                a_rte_fluxes,
                                a_extrap_cdr_fluxes,
                                m_townsend2_dielectric,
                                m_dielectric_quantum_efficiency);
}

Real
CdrPlasmaAir9EedBourdon::initialSigma(const Real a_time, const RealVect a_pos) const
{
  return 0.0;
}

Real
CdrPlasmaAir9EedBourdon::Electron_energy(const Real a_energy, const Real a_density)
{
  return a_energy / (1.0 + a_density);
}

Real
CdrPlasmaAir9EedBourdon::compute_eed_mobility(const Real a_energy) const
{
  return (5.0 / 3.0) * compute_e_mobility(a_energy);
}
Real
CdrPlasmaAir9EedBourdon::compute_e_mobility(const Real a_energy) const
{
  return m_e_mobility.getEntry<1>(a_energy);
}
Real
CdrPlasmaAir9EedBourdon::compute_N2plus_mobility(const Real a_EbyN) const
{
  return m_ion_mobility;
}
Real
CdrPlasmaAir9EedBourdon::compute_N4plus_mobility(const Real a_EbyN) const
{
  return m_ion_mobility;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2plus_mobility(const Real a_EbyN) const
{
  return m_ion_mobility;
}
Real
CdrPlasmaAir9EedBourdon::compute_O4plus_mobility(const Real a_EbyN) const
{
  return m_ion_mobility;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2plusN2_mobility(const Real a_EbyN) const
{
  return m_ion_mobility;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2minus_mobility(const Real a_EbyN) const
{
  return m_ion_mobility;
}
Real
CdrPlasmaAir9EedBourdon::compute_Ominus_mobility(const Real a_EbyN) const
{
  return m_ion_mobility;
}

Real
CdrPlasmaAir9EedBourdon::compute_eed_diffco(const Real a_energy) const
{
  return (5.0 / 3.0) * this->compute_e_diffco(a_energy);
}
Real
CdrPlasmaAir9EedBourdon::compute_e_diffco(const Real a_energy) const
{
  return (2.0 / 3.0) * a_energy * this->compute_e_mobility(a_energy);
}
Real
CdrPlasmaAir9EedBourdon::compute_N2plus_diffco() const
{
  return m_ion_diffusion;
}
Real
CdrPlasmaAir9EedBourdon::compute_N4plus_diffco() const
{
  return m_ion_diffusion;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2plus_diffco() const
{
  return m_ion_diffusion;
}
Real
CdrPlasmaAir9EedBourdon::compute_O4plus_diffco() const
{
  return m_ion_diffusion;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2plusN2_diffco() const
{
  return m_ion_diffusion;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2minus_diffco() const
{
  return m_ion_diffusion;
}
Real
CdrPlasmaAir9EedBourdon::compute_Ominus_diffco() const
{
  return m_ion_diffusion;
}

Real
CdrPlasmaAir9EedBourdon::compute_Electron_N2_alpha(const Real a_energy) const
{
  return m_e_N2_alpha.getEntry<1>(a_energy);
}
Real
CdrPlasmaAir9EedBourdon::compute_Electron_O2_alpha(const Real a_energy) const
{
  return m_e_O2_alpha.getEntry<1>(a_energy);
}
Real
CdrPlasmaAir9EedBourdon::compute_N2plus_N2_M_to_N4plus_M() const
{
  return 5.E-41;
}
Real
CdrPlasmaAir9EedBourdon::compute_N4plus_O2_to_O2_2N2() const
{
  return 2.5E-16;
}
Real
CdrPlasmaAir9EedBourdon::compute_N2plus_O2_to_O2plus_N2(const Real a_Tg) const
{
  return 1.05E-15 / sqrt(a_Tg);
}
Real
CdrPlasmaAir9EedBourdon::compute_O2plus_2N2_to_O2plusN2_N2(const Real a_Tg) const
{
  return 8.1E-38 / (a_Tg * a_Tg);
}
Real
CdrPlasmaAir9EedBourdon::compute_O2plusN2_N2_to_O2plus_2N2(const Real a_Tg) const
{
  return 14.8 * pow(a_Tg, -5.3) * exp(-2357.0 / a_Tg);
}
Real
CdrPlasmaAir9EedBourdon::compute_O2plusN2_O2_to_O4plus_N2() const
{
  return 1.E-15;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2plus_O2_M_to_O4plus_M(const Real a_Tg) const
{
  return 2.03E-34 * pow(a_Tg, -3.2);
}
Real
CdrPlasmaAir9EedBourdon::compute_e_O4plus_to_2O2(const Real a_Te) const
{
  return 2.42E-11 / (sqrt(a_Te));
}
Real
CdrPlasmaAir9EedBourdon::compute_e_O2plus_to_O2(const Real a_Te) const
{
  return 6.E-11 / a_Te;
}
Real
CdrPlasmaAir9EedBourdon::compute_e_2O2_to_O2minus_O2(const Real a_Te) const
{
  return 6.E-39 / a_Te;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2minus_O4plus_to_3O2() const
{
  return 1.E-13;
}
Real
CdrPlasmaAir9EedBourdon::compute_O2minus_O4plus_M_to_3O2_M(const Real a_Tg) const
{
  return 3.12E-31 * pow(a_Tg, -2.5);
}
Real
CdrPlasmaAir9EedBourdon::compute_O2minus_O2plus_M_to_2O2_M(const Real a_Tg) const
{
  return 3.12E-31 * pow(a_Tg, -2.5);
}
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_to_e_2O_c1(const Real a_energy) const
{
  Real k = 0.0;
  return 0.0;

  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.6937E-18;
  const Real max_energy_coeff = 0.8357E-14;

  if (a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if (a_energy > max_energy) {
    k = max_energy_coeff;
  }
  else {
    const Real A = -29.04;
    const Real B = 0.8249;
    const Real C = -17.40;
    const Real D = 7.278;
    const Real E = -2.656;

    const Real x = a_energy;
    k            = exp(A + B * log(x) + C / x + D / (x * x) + E / (x * x * x));
  }

  return k;
} // TABLE
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_to_e_2O_c2(const Real a_energy) const
{
  Real k = 0.0;
  return k;

  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.1149E-17;
  const Real max_energy_coeff = 0.2596E-12;

  if (a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if (a_energy > max_energy) {
    k = max_energy_coeff;
  }
  else {
    const Real A = -29.54;
    const Real B = 0.2660;
    const Real C = -10.13;
    const Real D = -2.733;
    const Real E = 1.101;

    const Real x = a_energy;
    k            = exp(A + B * log(x) + C / x + D / (x * x) + E / (x * x * x));
  }

  return k;
} // TABLE
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_to_Ominus_O(const Real a_energy) const
{
  Real k = 0.0;
  return k;

  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.1456E-17;
  const Real max_energy_coeff = 0.3192E-15;

  if (a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if (a_energy > max_energy) {
    k = max_energy_coeff;
  }
  else {
    const Real A = -37.26;
    const Real B = 0.4130;
    const Real C = 5.960;
    const Real D = -19.30;
    const Real E = 9.528;

    const Real x = a_energy;
    k            = exp(A + B * log(x) + C / x + D / (x * x) + E / (x * x * x));
  }

  return k;
} // TABLE
Real
CdrPlasmaAir9EedBourdon::compute_Oplus_O2_to_O_O2(const Real a_Tg) const
{
  return 3.46E-12 / sqrt(a_Tg);
}
Real
CdrPlasmaAir9EedBourdon::compute_e_N2_to_e_N2(const Real a_energy) const
{

  Real k = 0.0;
  return k;

  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.4298e-17;
  const Real max_energy_coeff = 0.3746e-15;

  if (a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if (a_energy > max_energy) {
    k = max_energy_coeff;
  }
  else {
    const Real A = -39.56;
    const Real B = 1.204;
    const Real C = -1.689;
    const Real D = 3.688;
    const Real E = -2.424;

    const Real x = a_energy;
    k            = exp(A + B * log(x) + C / x + D / (x * x) + E / (x * x * x));
  }

  return k;
} // TABLE
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_to_e_O2(const Real a_energy) const
{
  Real k = 0.0;
  return k;

  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.1561E-17;
  const Real max_energy_coeff = 0.3645E-15;

  if (a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if (a_energy > max_energy) {
    k = max_energy_coeff;
  }
  else {
    const Real A = -40.17;
    const Real B = 1.390;
    const Real C = -3.345;
    const Real D = 4.603;
    const Real E = -2.090;

    const Real x = a_energy;
    k            = exp(A + B * log(x) + C / x + D / (x * x) + E / (x * x * x));
  }

  return k;
} // TABLE

// Electron losses
Real
CdrPlasmaAir9EedBourdon::compute_e_N2_ionization_loss() const
{
  return 15.6;
}
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_ionization_loss() const
{
  return 12.07;
}
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_dissociation_loss_c1() const
{
  return 5.58;
}
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_dissociation_loss_c2() const
{
  return 8.4;
}
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_dissociative_attachment_loss() const
{
  return 3.6;
}
Real
CdrPlasmaAir9EedBourdon::compute_e_O2_scattering_loss() const
{
  return 1;
}
Real
CdrPlasmaAir9EedBourdon::compute_e_N2_scattering_loss() const
{
  return 1;
}

Real
CdrPlasmaAir9EedBourdon::initEed(const RealVect a_pos, const Real a_time, const RealVect a_E)
{
  const Real EbyN = (a_E / (m_N * Units::Td)).vectorLength();
  return m_initEed.getEntry<1>(EbyN) * m_cdrSpecies[m_Electron_idx]->initialData(a_pos, a_time);
}

Vector<Real>
CdrPlasmaAir9EedBourdon::computeCdrFluxes(const Real         a_time,
                                          const RealVect     a_pos,
                                          const RealVect     a_normal,
                                          const RealVect     a_E,
                                          const Vector<Real> a_cdr_densities,
                                          const Vector<Real> a_cdr_velocities,
                                          const Vector<Real> a_cdr_gradients,
                                          const Vector<Real> a_rte_fluxes,
                                          const Vector<Real> a_extrap_cdr_fluxes,
                                          const Real         a_townsend2,
                                          const Real         a_quantum_efficiency) const
{

  Vector<Real> fluxes(m_numCdrSpecies, 0.0);

#if 1 // debug
  return a_extrap_cdr_fluxes;
#endif

  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;

  const Real Electron_energy = a_cdr_densities[m_eed_idx] / (1.0 + a_cdr_densities[m_Electron_idx]);
  const Real Te              = 2.0 * Electron_energy * Units::Qe / (3.0 * Units::kb);
  const Real EbyN            = (a_E / (m_N * Units::Td)).vectorLength();
  const Real ion_mass        = 2.65E-26;                                        // kg
  const Real vth_g           = sqrt(Units::kb * m_Tg / (Units::pi * ion_mass)); // Ion thermal velocity
  const Real vth_e           = sqrt(Units::kb * Te / (Units::pi * Units::me));  // Electron thermal velocity

  // Switch for setting drift flux to zero for charge species
  Vector<Real> aj(m_numCdrSpecies, 0.0);
  for (int i = 0; i < m_numCdrSpecies; i++) {
    if (DataOps::sgn(m_cdrSpecies[i]->getChargeNumber()) * PolyGeom::dot(a_E, a_normal) < 0) {
      aj[i] = 1.0;
    }
    else {
      aj[i] = 0.0;
    }
  }

  // Drift outflow for now
  for (int i = 0; i < m_numCdrSpecies; i++) {
    //fluxes[i] = Max(0.0, aj[i]*a_extrap_cdr_fluxes[i]);
    //    fluxes[i] = Max(0.0, a_extrap_cdr_fluxes[i]);
  }

  return fluxes;
}

Vector<Real>
CdrPlasmaAir9EedBourdon::computeCdrElectrodeFluxes(const Real         a_time,
                                                   const RealVect     a_pos,
                                                   const RealVect     a_normal,
                                                   const RealVect     a_E,
                                                   const Vector<Real> a_cdr_densities,
                                                   const Vector<Real> a_cdr_velocities,
                                                   const Vector<Real> a_cdr_gradients,
                                                   const Vector<Real> a_rte_fluxes,
                                                   const Vector<Real> a_extrap_cdr_fluxes) const
{

  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;

  Vector<Real> ret(m_numCdrSpecies, 0.0);
  if (cathode) {
    ret = this->computeCathodeFluxes(a_time,
                                     a_pos,
                                     a_normal,
                                     a_E,
                                     a_cdr_densities,
                                     a_cdr_velocities,
                                     a_cdr_gradients,
                                     a_rte_fluxes,
                                     a_extrap_cdr_fluxes,
                                     m_townsend2_electrode,
                                     m_electrode_quantum_efficiency);
  }
  else if (anode) {
    ret = this->computeAnodeFluxes(a_time,
                                   a_pos,
                                   a_normal,
                                   a_E,
                                   a_cdr_densities,
                                   a_cdr_velocities,
                                   a_cdr_gradients,
                                   a_rte_fluxes,
                                   a_extrap_cdr_fluxes,
                                   m_townsend2_electrode,
                                   m_electrode_quantum_efficiency);
  }

  return ret;
}

Vector<Real>
CdrPlasmaAir9EedBourdon::computeAnodeFluxes(const Real         a_time,
                                            const RealVect     a_pos,
                                            const RealVect     a_normal,
                                            const RealVect     a_E,
                                            const Vector<Real> a_cdr_densities,
                                            const Vector<Real> a_cdr_velocities,
                                            const Vector<Real> a_cdr_gradients,
                                            const Vector<Real> a_rte_fluxes,
                                            const Vector<Real> a_extrap_cdr_fluxes,
                                            const Real         a_townsend2,
                                            const Real         a_quantum_efficiency) const
{
  return a_extrap_cdr_fluxes;
}

Vector<Real>
CdrPlasmaAir9EedBourdon::computeCathodeFluxes(const Real         a_time,
                                              const RealVect     a_pos,
                                              const RealVect     a_normal,
                                              const RealVect     a_E,
                                              const Vector<Real> a_cdr_densities,
                                              const Vector<Real> a_cdr_velocities,
                                              const Vector<Real> a_cdr_gradients,
                                              const Vector<Real> a_rte_fluxes,
                                              const Vector<Real> a_extrap_cdr_fluxes,
                                              const Real         a_townsend2,
                                              const Real         a_quantum_efficiency) const
{

  Vector<Real> ret(m_numCdrSpecies, 0.0);

  // Drift outflow for positive species. No inflow for the others just yet.
  for (int i = 0; i < m_numCdrSpecies; i++) {
    ret[i] = (m_cdrSpecies[i]->getChargeNumber() > 0) ? a_extrap_cdr_fluxes[i] : 0.0;
  }

  // Electron inflow due to ion impingement
  ret[m_Electron_idx] = 0.0;
  ret[m_eed_idx]      = 0.0;
  for (int i = 0; i < m_numCdrSpecies; i++) {
    const int q = m_cdrSpecies[i]->getChargeNumber();
    if (q > 0) {
      const Real flx = a_extrap_cdr_fluxes[i] * m_townsend2_electrode;
      ret[m_Electron_idx] -= flx;
      ret[m_eed_idx] -= (5. / 3.) * m_cathode_work * flx;
    }
  }

  // Secondary emission from Photons
  for (int j = 0; j < m_numRtSpecies; j++) {
  }

  return ret;
  //
  //  MayDay::Abort("CdrPlasmaAir9EedBourdon::computeCathodeFluxes - negative streamer BC not implemented");
}

#include <CD_NamespaceFooter.H>
