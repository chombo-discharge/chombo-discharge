/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir3Bourdon.cpp
  @brief  Implementation of CD_CdrPlasmaAir3Bourdon.H
  @author Robert Marskar
*/

// Std includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

// Chombo includes
#include <ParmParse.H>
#include <PolyGeom.H>

// Our includes
#include <CD_CdrPlasmaAir3Bourdon.H>
#include <CD_CdrPlasmaAir3BourdonSpecies.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_DataParser.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

std::string CdrPlasmaAir3Bourdon::s_bolsig_mobility = "# Electron mobility (E/N, mu*N)";
std::string CdrPlasmaAir3Bourdon::s_bolsig_diffco   = "# Electron diffusion coefficient (E/N, D*N)";
std::string CdrPlasmaAir3Bourdon::s_bolsig_alpha    = "# Townsend ionization coeff (E/N, alpha/N)";
std::string CdrPlasmaAir3Bourdon::s_bolsig_eta      = "# Townsend attachment coeff (E/N, eta/N)";
std::string CdrPlasmaAir3Bourdon::s_bolsig_energy   = "# Electron mean energy (E/N, eV)";

CdrPlasmaAir3Bourdon::CdrPlasmaAir3Bourdon() {

  ParmParse pp("CdrPlasmaAir3Bourdon");

  // Get input parameters
  pp.get("pressure",                      m_p);
  pp.get("quenching_pressure",            m_pq);
  pp.get("temperature",                   m_T);
  pp.get("transport_file",                m_transport_file);
  pp.get("uniform_tables",                m_uniform_entries);
  pp.get("use_alpha_corr",                m_alpha_corr);      
  pp.get("mobile_electrons",              m_isMobile_Electrons);    
  pp.get("diffusive_electrons",           m_isDiffusive_Electrons);
  pp.get("diffusive_ions",                m_isDiffusive_ions);
  pp.get("mobile_ions",                   m_isMobile_ions);
  pp.get("ion_mobility",                  m_ion_mobility);

  pp.get("excitation_efficiency",         m_photoexc_eff);
  pp.get("photoionization_efficiency",    m_photoi_eff);

  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);

  // Normalize units.
  m_N2frac = 0.8;
  m_O2frac = 0.2;
  m_p     *= Units::atm2pascal;
  m_pq    *= Units::atm2pascal;
  m_N      = m_p*Units::Na/(m_T*Units::R);
  
  m_ion_diffusion = m_ion_mobility*(Units::kb*m_T)/Units::Qe; // Einstein relation. 


  // Read from input file and put in lookup tables. 
  std::ifstream infile(m_transport_file);
  if(infile.good()){
    infile.close();

    constexpr Real eVtoK = 2.0*Units::Qe/(3.0*Units::kb);

    // Read data into LookupTables with 2 columns. After reading, we call LookupTable::sort which is a kind of requirement
    // that ensures that we can make the table uniform and used constant spacing. Note that the input electric field is
    // in Td so we also need to scale by N * 1.E-21. Likewise, transport coefficients are mu*N while ionization and
    // attachment coefficients are alpha/N. So, scale that out as well. 
    m_e_mobility = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Bourdon::s_bolsig_mobility, "");
    m_e_diffco   = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Bourdon::s_bolsig_diffco,   "");
    m_e_alpha    = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Bourdon::s_bolsig_alpha,    "");
    m_e_eta      = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Bourdon::s_bolsig_eta,      "");
    m_e_energy   = DataParser::fractionalFileReadASCII(m_transport_file, CdrPlasmaAir3Bourdon::s_bolsig_energy,   "");    

    m_e_mobility.sort();
    m_e_diffco.  sort();
    m_e_alpha.   sort();
    m_e_eta.     sort();
    m_e_energy.  sort();    
    
    m_e_mobility.scale<0>(m_N*Units::Td);
    m_e_diffco.  scale<0>(m_N*Units::Td);
    m_e_alpha.   scale<0>(m_N*Units::Td);
    m_e_eta.     scale<0>(m_N*Units::Td);
    m_e_energy.  scale<0>(m_N*Units::Td);    
    
    m_e_mobility.scale<1>(1./m_N); 
    m_e_diffco.  scale<1>(1./m_N); 
    m_e_alpha.   scale<1>(   m_N); 
    m_e_eta.     scale<1>(   m_N); 
    m_e_energy.  scale<1>( eVtoK);    

    m_e_diffco.  makeUniform(m_uniform_entries);
    m_e_mobility.makeUniform(m_uniform_entries);
    m_e_alpha.   makeUniform(m_uniform_entries);
    m_e_eta.     makeUniform(m_uniform_entries);
    m_e_energy.  makeUniform(m_uniform_entries);
  }
  else{
    MayDay::Abort("CdrPlasmaAir3Bourdon::parseTransportFile - could not find transport data");
  }


  // Parse domain boundary conditions. 
  parseDomainBc();

  // Instantiate species. 
  initSpecies();
}

CdrPlasmaAir3Bourdon::~CdrPlasmaAir3Bourdon() {

}

void CdrPlasmaAir3Bourdon::initSpecies(){
  m_numCdrSpecies = 3;
  m_numRtSpecies = 3;

  m_elec_idx = 0;
  m_plus_idx = 1;
  m_minu_idx = 2;
  m_pho1_idx = 0;
  m_pho2_idx = 1;
  m_pho3_idx = 2;

  m_CdrSpecies.resize(m_numCdrSpecies);
  m_CdrSpecies[m_elec_idx]  = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir3Bourdon::Electron());
  m_CdrSpecies[m_plus_idx]  = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir3Bourdon::MPlus());
  m_CdrSpecies[m_minu_idx]  = RefCountedPtr<CdrSpecies> (new CdrPlasmaAir3Bourdon::MMinus());

  m_RtSpecies.resize(m_numRtSpecies);
  m_RtSpecies[m_pho1_idx] = RefCountedPtr<RtSpecies> (new CdrPlasmaAir3Bourdon::PhotonOne());
  m_RtSpecies[m_pho2_idx] = RefCountedPtr<RtSpecies> (new CdrPlasmaAir3Bourdon::PhotonTwo());
  m_RtSpecies[m_pho3_idx] = RefCountedPtr<RtSpecies> (new CdrPlasmaAir3Bourdon::PhotonThree());
}

void CdrPlasmaAir3Bourdon::parseDomainBc(){

  ParmParse pp("CdrPlasmaAir3Bourdon");
  std::string str;

  m_wallBc.resize(2*SpaceDim, 0); 
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
	
      std::string str_dir;
      if(dir == 0){
	str_dir = "x";
      }
      else if(dir == 1){
	str_dir = "y";
      }
      else if(dir == 2){
	str_dir = "z";
      }

      // Check for wall BCs
      if(side == Side::Lo){
	std::string type;
	std::string bc_string = "domain_bc_" + str_dir + "_lo";
	if(pp.contains(bc_string.c_str())){
	  pp.get(bc_string.c_str(), type);
	  const int idx = 2*dir;
	  if(type == "wall"){
	    m_wallBc[idx] = 1;
	  }
	}
      }
      else if(side == Side::Hi){
	std::string type;
	std::string bc_string = "domain_bc_" + str_dir + "_hi";
	if(pp.contains(bc_string.c_str())){
	  pp.get(bc_string.c_str(), type);
	  const int idx = 2*dir + 1;
	  if(type == "wall"){
	    m_wallBc[idx] = 1;
	  }
	}
      }
    }
  }
}

void CdrPlasmaAir3Bourdon::advanceReactionNetwork(Vector<Real>&          a_particle_sources,
						  Vector<Real>&          a_Photon_sources,
						  const Vector<Real>     a_particle_densities,
						  const Vector<RealVect> a_particle_gradients,
						  const Vector<Real>     a_Photon_densities,
						  const RealVect         a_E,
						  const RealVect         a_pos,
						  const Real             a_dx,
						  const Real             a_dt,
						  const Real             a_time,
						  const Real             a_kappa) const {
  const Real E      = a_E.vectorLength();
  const Real ve     = E*m_e_mobility.getEntry<MU>(E);
  
  // Ionization and attachment coefficients
  Real alpha  = m_e_alpha.getEntry<ALPHA>(E);
  Real eta    = m_e_eta.getEntry<ETA>(E);

  // Modify alpha
  if(m_alpha_corr){
    const RealVect Eunit = a_E/a_E.vectorLength();
    const Real De        = m_e_diffco.getEntry<DIFFCO>(E);
    const RealVect gNe   = a_particle_gradients[m_elec_idx];

    Real fcorr = 1.0;
    fcorr = 1.0 + PolyGeom::dot(Eunit, De*gNe)/(1.0+a_particle_densities[m_elec_idx]*ve);
    fcorr = Min(fcorr, 1.0);
    fcorr = Max(0.0, fcorr);

    alpha = alpha*fcorr;
  }

  const Real Te = m_e_energy.getEntry<1>(E);

  const Real R1 = alpha*ve*a_particle_densities[m_elec_idx];
  const Real R2 = eta*ve*a_particle_densities[m_elec_idx];
  const Real R3 = 4.2E-12*sqrt(300/Te) * a_particle_densities[m_elec_idx] * a_particle_densities[m_plus_idx];

  Real& Se = a_particle_sources[m_elec_idx];
  Real& Sp = a_particle_sources[m_plus_idx];
  Real& Sm = a_particle_sources[m_minu_idx];

  Se = 0.0;
  Sp = 0.0;
  Sm = 0.0;

  // e + M => e + e + M+
  Se += R1;
  Sp += R1;

  // e + M => M-
  Se -= R2;
  Sm += R2;

  // e + M+ => M
  Se -= R3;
  Sp -= R3;

  // Photoionization, M + y => e + M+
  const CdrPlasmaAir3Bourdon::PhotonOne*   Photon1 = static_cast<CdrPlasmaAir3Bourdon::PhotonOne*>   (&(*m_RtSpecies[m_pho1_idx]));
  const CdrPlasmaAir3Bourdon::PhotonTwo*   Photon2 = static_cast<CdrPlasmaAir3Bourdon::PhotonTwo*>   (&(*m_RtSpecies[m_pho2_idx]));
  const CdrPlasmaAir3Bourdon::PhotonThree* Photon3 = static_cast<CdrPlasmaAir3Bourdon::PhotonThree*> (&(*m_RtSpecies[m_pho3_idx]));

  const Real Sph = m_photoi_eff*Units::c*m_O2frac*m_p*(Photon1->get_A()*a_Photon_densities[m_pho1_idx]
						       + Photon2->get_A()*a_Photon_densities[m_pho2_idx]
						       + Photon3->get_A()*a_Photon_densities[m_pho3_idx]);
  Se += Sph;
  Sp += Sph;

  // Photo-excitation and emission. 
  const Real quench  = m_pq/(m_p + m_pq);
  const Real Rgamma  = R1*m_photoexc_eff*quench;
  a_Photon_sources[m_pho1_idx] = Rgamma;
  a_Photon_sources[m_pho2_idx] = Rgamma;
  a_Photon_sources[m_pho3_idx] = Rgamma;

  CH_assert(std::abs(Se + Sp + Sm) <= 1.E-10);    

  return;
}

Vector<Real> CdrPlasmaAir3Bourdon::computeCdrDiffusionCoefficients(const Real         a_time,
								   const RealVect     a_pos,
								   const RealVect     a_E,
								   const Vector<Real> a_cdr_densities) const {

  Vector<Real> dco(m_numCdrSpecies, 0.0);
  dco[m_elec_idx] = m_e_diffco.getEntry<DIFFCO>(a_E.vectorLength());
  dco[m_plus_idx] = m_ion_diffusion;
  dco[m_minu_idx] = m_ion_diffusion;
  
  return dco;
}
  
Vector<RealVect> CdrPlasmaAir3Bourdon::computeCdrDriftVelocities(const Real         a_time,
								 const RealVect     a_pos,
								 const RealVect     a_E,
								 const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_numCdrSpecies, RealVect::Zero);

  vel[m_elec_idx] = -a_E*m_e_mobility.getEntry<MU>(a_E.vectorLength());
  vel[m_plus_idx] =  a_E*m_ion_mobility;
  vel[m_minu_idx] = -a_E*m_ion_mobility;
  
  return vel;
}
  
Vector<Real> CdrPlasmaAir3Bourdon::computeCdrDomainFluxes(const Real           a_time,
							  const RealVect       a_pos,
							  const int            a_dir,
							  const Side::LoHiSide a_side,
							  const RealVect       a_E,
							  const Vector<Real>   a_cdr_densities,
							  const Vector<Real>   a_cdr_velocities,
							  const Vector<Real>   a_cdr_gradients,
							  const Vector<Real>   a_rte_fluxes,
							  const Vector<Real>   a_extrap_cdr_fluxes) const{
  Vector<Real> fluxes(m_numCdrSpecies, 0.0);

  int idx, sgn;
  if(a_side == Side::Lo){
    sgn = -1;
    idx = 2*a_dir;
  }
  else{
    sgn = 1;
    idx = 2*a_dir + 1;
  }

  if(m_wallBc[idx] == 0){ // Inflow/outflow
    for (int i = 0; i < fluxes.size(); i++){
      fluxes[i] = a_extrap_cdr_fluxes[i];
    }
  }
  else if(m_wallBc[idx] == 1){ // wall
    for (int i = 0; i < fluxes.size(); i++){
      fluxes[i] = 0.0;
    }
  }
  else{
    MayDay::Abort("morrow_jiang::computeCdrDomainFluxes - uknown domain bc requested");
  }

  
  return fluxes;
}
  
Vector<Real> CdrPlasmaAir3Bourdon::computeCdrElectrodeFluxes(const Real         a_time,
							     const RealVect     a_pos,
							     const RealVect     a_normal,
							     const RealVect     a_E,
							     const Vector<Real> a_cdr_densities,
							     const Vector<Real> a_cdr_velocities,
							     const Vector<Real> a_cdr_gradients,
							     const Vector<Real> a_rte_fluxes,
							     const Vector<Real> a_extrap_cdr_fluxes) const{
  return computeCdrFluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			  a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
}

Vector<Real> CdrPlasmaAir3Bourdon::computeCdrDielectricFluxes(const Real         a_time,
							      const RealVect     a_pos,
							      const RealVect     a_normal,
							      const RealVect     a_E,
							      const Vector<Real> a_cdr_densities,
							      const Vector<Real> a_cdr_velocities,
							      const Vector<Real> a_cdr_gradients,
							      const Vector<Real> a_rte_fluxes,
							      const Vector<Real> a_extrap_cdr_fluxes) const{
  return computeCdrFluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			  a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

Vector<Real> CdrPlasmaAir3Bourdon::computeCdrFluxes(const Real         a_time,
						    const RealVect     a_pos,
						    const RealVect     a_normal,
						    const RealVect     a_E,
						    const Vector<Real> a_cdr_densities,
						    const Vector<Real> a_cdr_velocities,
						    const Vector<Real> a_cdr_gradients,
						    const Vector<Real> a_rte_fluxes,
						    const Vector<Real> a_extrap_cdr_fluxes,
						    const Real         a_townsend2,
						    const Real         a_quantum_efficiency) const{
  Vector<Real> fluxes(m_numCdrSpecies, 0.0);

  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;

  // Switch for setting drift flux to zero for charge species
  Vector<Real> aj(m_numCdrSpecies, 0.0);
  for (int i = 0; i < m_numCdrSpecies; i++){
    if(DataOps::sgn(m_CdrSpecies[i]->getChargeNumber())*PolyGeom::dot(a_E, a_normal) < 0){
      aj[i] = 1.0;
    }
    else {
      aj[i] = 0.0;
    }
  }

  // Drift outflow for now
  for (int i = 0; i < m_numCdrSpecies; i++){
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  return fluxes;
}

Real CdrPlasmaAir3Bourdon::initialSigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

Real CdrPlasmaAir3Bourdon::computeAlpha(const RealVect a_E) const{
  const Real E     = a_E.vectorLength();
  const Real alpha = m_e_alpha.getEntry<ALPHA>(E);

  return alpha;
}

#include <CD_NamespaceFooter.H>
