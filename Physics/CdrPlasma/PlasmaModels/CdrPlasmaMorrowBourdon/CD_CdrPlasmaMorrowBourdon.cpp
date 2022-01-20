/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaMorrowBourdon.cpp
  @brief  Implementation of CD_CdrPlasmaMorrowBourdon.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <PolyGeom.H>

// Our includes
#include <CD_CdrPlasmaMorrowBourdon.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaMorrowBourdon::CdrPlasmaMorrowBourdon(){
  CH_TIME("CdrPlasmaMorrowBourdon::CdrPlasmaMorrowBourdon");

  parseGas();
  parsePhotoi();
  parseSEE();
  parseBc();

  initSpecies();
}

CdrPlasmaMorrowBourdon::~CdrPlasmaMorrowBourdon(){


}

void CdrPlasmaMorrowBourdon::parseGas(){
  ParmParse pp("CdrPlasmaMorrowBourdon");
  pp.get("gas_temperature", m_temp);
  pp.get("gas_N2_frac",     m_fracN2);
  pp.get("gas_O2_frac",     m_fracO2);
  pp.get("gas_pressure",    m_p);

  // Convert to correct units and compute necessary things
  m_p  *= Units::atm2pascal;
  m_N   = m_p*Units::Na/(m_temp*Units::R);
}

void CdrPlasmaMorrowBourdon::parsePhotoi(){
  ParmParse pp("CdrPlasmaMorrowBourdon");
  pp.get("gas_quenching_pressure",     m_pq);
  pp.get("excitation_efficiency",      m_exc_eff);
  pp.get("photoionization_efficiency", m_photo_eff);

  m_pq *= Units::atm2pascal;
}

void CdrPlasmaMorrowBourdon::parseSEE(){
  ParmParse pp("CdrPlasmaMorrowBourdon");
  pp.get("electrode_townsend2",           m_townsend2_conductor);
  pp.get("electrode_quantum_efficiency",  m_electrode_yield);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("dielectric_quantum_efficiency", m_dielectric_yield);
}

void CdrPlasmaMorrowBourdon::parseBc(){

  m_wallBc.resize(2*SpaceDim, 0); 
  ParmParse pp("CdrPlasmaMorrowBourdon");
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();

      // Identifier
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


      if(side == Side::Lo){
	std::string type;
	std::string bc_string = "domain_bc_" + str_dir + "_lo";
	pp.get(bc_string.c_str(), type);
	const int idx = 2*dir;
	if(type == "wall"){
	  m_wallBc[idx] = 1;
	}
      }
      else if(side == Side::Hi){
	std::string type;
	std::string bc_string = "domain_bc_" + str_dir + "_hi";
	pp.get(bc_string.c_str(), type);
	const int idx = 2*dir + 1;
	if(type == "wall"){
	  m_wallBc[idx] = 1;
	}
      }
    }
  }
}

void CdrPlasmaMorrowBourdon::initSpecies(){
  m_numCdrSpecies = 3;
  m_numRtSpecies = 3;

  m_CdrSpecies.resize(m_numCdrSpecies);
  m_RtSpecies.resize(m_numRtSpecies);

  m_nelec_idx   = 0;
  m_nplus_idx   = 1;
  m_nminu_idx   = 2;
  
  m_Photon1_idx = 0;
  m_Photon2_idx = 1;
  m_Photon3_idx = 2;
  
  m_CdrSpecies[m_nelec_idx]    = RefCountedPtr<CdrSpecies> (new CdrPlasmaMorrowBourdon::Electron());
  m_CdrSpecies[m_nplus_idx]    = RefCountedPtr<CdrSpecies> (new CdrPlasmaMorrowBourdon::PositiveSpecies());
  m_CdrSpecies[m_nminu_idx]    = RefCountedPtr<CdrSpecies> (new CdrPlasmaMorrowBourdon::negative_species());
  
  m_RtSpecies[m_Photon1_idx]  = RefCountedPtr<RtSpecies> (new CdrPlasmaMorrowBourdon::PhotonOne());
  m_RtSpecies[m_Photon2_idx]  = RefCountedPtr<RtSpecies> (new CdrPlasmaMorrowBourdon::PhotonTwo());
  m_RtSpecies[m_Photon3_idx]  = RefCountedPtr<RtSpecies> (new CdrPlasmaMorrowBourdon::PhotonThree());
}

void CdrPlasmaMorrowBourdon::advanceReactionNetwork(Vector<Real>&          a_particle_sources,
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

  const Real E = a_E.vectorLength();
  const Real alpha  = computeAlpha(E, a_pos); // Ionization coefficient
  const Real eta    = computeEta(a_E);   // Attachment coefficient
  const Real beta   = computeBeta(a_E);  // Recombination coefficient

  // Cast so we can get A-coefficients
  const CdrPlasmaMorrowBourdon::PhotonOne*   Photon1 = static_cast<CdrPlasmaMorrowBourdon::PhotonOne*>   (&(*m_RtSpecies[m_Photon1_idx]));
  const CdrPlasmaMorrowBourdon::PhotonTwo*   Photon2 = static_cast<CdrPlasmaMorrowBourdon::PhotonTwo*>   (&(*m_RtSpecies[m_Photon2_idx]));
  const CdrPlasmaMorrowBourdon::PhotonThree* Photon3 = static_cast<CdrPlasmaMorrowBourdon::PhotonThree*> (&(*m_RtSpecies[m_Photon3_idx]));
  
  // Densities and velocities
  const Real Ne  = a_particle_densities[m_nelec_idx]; 
  const Real Np  = a_particle_densities[m_nplus_idx];
  const Real Nn  = a_particle_densities[m_nminu_idx];
  const Real Ve  = computeVe(a_E).vectorLength();
  const Real Sph = m_photo_eff*Units::c*m_fracO2*m_p*(Photon1->getA()*a_Photon_densities[m_Photon1_idx]
						      + Photon2->getA()*a_Photon_densities[m_Photon2_idx]
						      + Photon3->getA()*a_Photon_densities[m_Photon3_idx]);


  Real& Se = a_particle_sources[m_nelec_idx];
  Real& Sp = a_particle_sources[m_nplus_idx];
  Real& Sn = a_particle_sources[m_nminu_idx];

  Se = alpha*Ne*Ve - eta*Ne*Ve   - beta*Ne*Np + Sph;
  Sp = alpha*Ne*Ve - beta*Np*Nn  - beta*Ne*Np + Sph;
  Sn = eta*Ne*Ve   - beta*Np*Nn;

  const Real tmp = std::max(0.0, alpha*Ne*Ve*m_exc_eff*(m_pq/(m_pq + m_p)));

  a_Photon_sources[m_Photon1_idx] = tmp;
  a_Photon_sources[m_Photon2_idx] = tmp;
  a_Photon_sources[m_Photon3_idx] = tmp;
}

Vector<RealVect> CdrPlasmaMorrowBourdon::computeCdrDriftVelocities(const Real         a_time,
								   const RealVect     a_pos,
								   const RealVect     a_E,
								   const Vector<Real> a_cdr_densities) const {

  Vector<RealVect> velocities(m_numCdrSpecies);
  
  velocities[m_nelec_idx] = computeVe(a_E);
  velocities[m_nplus_idx] = computeVp(a_E);
  velocities[m_nminu_idx] = computeVn(a_E);

  return velocities;
}

RealVect CdrPlasmaMorrowBourdon::computeVe(const RealVect a_E) const{
  RealVect ve = RealVect::Zero;

  const RealVect E = a_E*1.E-2;          // Morrow-Lowke wants E in V/cm
  const Real Emag  = E.vectorLength();   //
  const Real N     = m_N*1.E-6;          // Morrow-Lowke wants N in cm^3
  const Real EbyN  = Emag/N;             //

  const Real lim0 = 2.6E-17;
  const Real lim1 = 1.E-16;
  const Real lim2 = 2.0E-15;
  const Real E_SI = a_E.vectorLength();
  if(EbyN <= lim0){
    ve = -E/Emag*(6.87E22*EbyN + 3.38E4);
  }
  else if(EbyN > lim0 && EbyN <= lim1){
    ve = -E/Emag*(7.293E21*EbyN + 1.63E6);
  }
  else if(EbyN > lim1 && EbyN <= lim2){
    ve = -E/Emag*(1.03E22*EbyN + 1.3E6);
  }
  else if(EbyN > lim2){
    ve = -E/Emag*(7.4E21*EbyN + 7.1E6);
  }

  ve *= 0.01; // Morrow-Lowke expressions are in cm/s
  return ve;
}

RealVect CdrPlasmaMorrowBourdon::computeVp(const RealVect a_E) const{
  const RealVect E = a_E*1.E-2;           // E in V/cm
  RealVect vp = 2.34*E*m_p/Units::atm2pascal;  // Morrow-Lowke wants V/cm
  vp *= 0.01;                             // Morrow-Lowke expression is in cm/s

  return vp;  
}

RealVect CdrPlasmaMorrowBourdon::computeVn(const RealVect a_E) const{
  RealVect vn = RealVect::Zero;

  const RealVect E = a_E*1.E-2;       // Morrow-Lowke wants E in V/cm
  const Real Emag = E.vectorLength(); // 
  const Real N    = m_N*1.E-6;        // Morrow-Lowke weants N in cm^3
  const Real EbyN = Emag/N;           //

  const Real lim0 = 5.0E-16;
  if(EbyN <= lim0){
    vn = -2.7*E*m_p/Units::atm2pascal;
  }
  else{
    vn = -1.86*E*m_p/Units::atm2pascal;
  }

  vn *= 0.01; // Morrow-Lowke expression is in cm/s
  return vn;
}

Real CdrPlasmaMorrowBourdon::computeAlpha(const Real a_E, const RealVect a_pos) const{
  Real alpha    = 0.;
  Real alphabyN = 0.;

  const Real Emag  = a_E*1.E-2;        // Morrow-Lowke wants E in V/cm
  const Real N     = m_N*1.E-6;        // Morrow-Lowke wants N in cm^3
  const Real EbyN  = Emag/N;           // This is now V/cm^2

  const Real lim = 1.5E-15;
  if(EbyN > lim){
    alphabyN = 2.0E-16*exp(-7.248E-15/EbyN);
  }
  else{
    alphabyN = 6.619E-17*exp(-5.593E-15/EbyN);
  }

  alphabyN *= 1.E-4; // Morrow-Lowke expression is in cm^2
  alpha = alphabyN*m_N;

  return alpha;
}

Real CdrPlasmaMorrowBourdon::computeEta(const RealVect a_E) const{

  const Real eta2 = this->computeEta2(a_E); 
  const Real eta3 = this->computeEta3(a_E);
  const Real eta  = eta2 + eta3;

  return eta;
}

Real CdrPlasmaMorrowBourdon::computeEta2(const RealVect a_E) const{
  Real eta2    = 0.;
  Real eta2byN = 0.;

  //
  const RealVect E = a_E*1.E-2;        // Morrow-Lowke wants E in V/cm
  const Real Emag  = E.vectorLength(); //
  const Real N     = m_N*1.E-6;        // Morrow-Lowke weants N in cm^3
  const Real EbyN  = Emag/N;           // This is now V/cm^2

  const Real lim = 1.05E-15;
  if(EbyN > lim){
    eta2byN = 8.889E-5*EbyN + 2.567E-19;
  }
  else{
    eta2byN = 6.089E-4*EbyN - 2.893E-19;
  }

  eta2byN *= 1.E-4;   // Morrow-Lowke expression is in cm^2, make it m^2
  eta2 = eta2byN*m_N; //

  return eta2;
}

Real CdrPlasmaMorrowBourdon::computeEta3(const RealVect a_E) const{
  const RealVect E = a_E*1.E-2;         // Morrow-Lowke wants E in V/cm
  const Real Emag  = E.vectorLength();  //
  const Real N     = m_N*1.E-6;         // Morrow-Lowke weants N in cm^3
  const Real EbyN  = Emag/N;            // This is now V/cm^2

  //
  Real eta3    = 0.;
  Real eta3byN = 4.7778E-59*pow(EbyN, -1.2749);

  //
  eta3byN *= 1.E-10; // Morrow-Lowke expression is in cm^5. Make it m^5
  eta3 = eta3byN*m_N*m_N;

  return eta3;
}

Real CdrPlasmaMorrowBourdon::computeBeta(const RealVect a_E) const{
  Real beta = 2.0E-7;
  beta *= 1.E-6; // Morrow-Lowke expression is in cm^3. Make it m^3
  return beta;
}

Real CdrPlasmaMorrowBourdon::computeDe(const RealVect a_E) const{
  const RealVect E  = a_E*1.E-2;                 // Morrow-Lowke wants E in V/cm
  const Real Emag   = E.vectorLength();          //
  const Real N      = m_N*1.E-6;                 // Morrow-Lowke weants N in cm^3
  const Real EbyN   = Emag/N;                    //

  //
  const RealVect Ve = this->computeVe(a_E);     // Does it's own conversion, comes out in m/s (aka. mentally sane units)
  const Real ve     = Ve.vectorLength()*1.E-2;   // Make it cm/s
  Real De = 0.3341E9*pow(EbyN, 0.54069)*ve/(1.E-4 + Emag);
  
  De *= 1.E-4; // Morrow-Lowke expression is in cm^2/s. Make it m^2/s
  
  return De;
}

Vector<Real> CdrPlasmaMorrowBourdon::computeCdrDiffusionCoefficients(const Real         a_time,
								     const RealVect     a_pos,
								     const RealVect     a_E,
								     const Vector<Real> a_cdr_densities) const {

  Vector<Real> diffCo(m_numCdrSpecies, 0.0);
  diffCo[m_nelec_idx] = computeDe(a_E);
  
  return diffCo;
}

Vector<Real> CdrPlasmaMorrowBourdon::computeCdrDielectricFluxes(const Real         a_time,
								const RealVect     a_pos,
								const RealVect     a_normal,
								const RealVect     a_E,
								const Vector<Real> a_cdr_densities,
								const Vector<Real> a_cdr_velocities,
								const Vector<Real> a_cdr_gradients,
								const Vector<Real> a_rte_fluxes,
								const Vector<Real> a_extrap_cdr_fluxes) const {
  // Outflux of species
  Vector<Real> fluxes(m_numCdrSpecies, 0.0);

  if(PolyGeom::dot(a_E, a_normal) > 0.0){ // Field points into gas phase
    fluxes[m_nelec_idx] = Max(0.0, a_extrap_cdr_fluxes[m_nelec_idx]); // Outflow for Electrons
    fluxes[m_nminu_idx] = Max(0.0, a_extrap_cdr_fluxes[m_nminu_idx]); // Outflow for negative species
  }
  else if(PolyGeom::dot(a_E, a_normal) < 0.0){ // Field points into dielectric
    fluxes[m_nplus_idx] = Max(0.0, a_extrap_cdr_fluxes[m_nplus_idx]); // Outflow for positive species
  }
  
  // Add in photoelectric effect and ion bombardment for Electrons by positive ions
  if(PolyGeom::dot(a_E, a_normal) < 0.){
    fluxes[m_nelec_idx] += -a_rte_fluxes[m_Photon1_idx]*m_dielectric_yield;
    fluxes[m_nelec_idx] += -a_rte_fluxes[m_Photon2_idx]*m_dielectric_yield;
    fluxes[m_nelec_idx] += -a_rte_fluxes[m_Photon3_idx]*m_dielectric_yield;
    fluxes[m_nelec_idx] += -Max(0.0, a_extrap_cdr_fluxes[m_nplus_idx])*m_townsend2_dielectric;
  }


  return fluxes;
}

Vector<Real> CdrPlasmaMorrowBourdon::computeCdrElectrodeFluxes(const Real         a_time,
							       const RealVect     a_pos,
							       const RealVect     a_normal,
							       const RealVect     a_E,
							       const Vector<Real> a_cdr_densities,
							       const Vector<Real> a_cdr_velocities,
							       const Vector<Real> a_cdr_gradients,
							       const Vector<Real> a_rte_fluxes,
							       const Vector<Real> a_extrap_cdr_fluxes) const {

  Vector<Real> fluxes(m_numCdrSpecies, 0.0);

  // Treat anode and cathode differently
  const bool is_cathode = PolyGeom::dot(a_E, a_normal) < 0.;
  const bool is_anode   = PolyGeom::dot(a_E, a_normal) > 0.;
  if(is_cathode){
    fluxes = this->computeCathodeFlux(a_extrap_cdr_fluxes,
				      a_cdr_densities,
				      a_cdr_velocities,
				      a_rte_fluxes,
				      a_E,
				      a_pos,
				      a_normal,
				      a_time);
  }
  else if(is_anode){
    fluxes = this->computeAnodeFlux(a_extrap_cdr_fluxes,
				    a_cdr_densities,
				    a_cdr_velocities,
				    a_rte_fluxes,
				    a_E,
				    a_pos,
				    a_normal,
				    a_time);
  }

  return fluxes;
}


Vector<Real> CdrPlasmaMorrowBourdon::computeCathodeFlux(const Vector<Real> a_extrapolated_fluxes,
							const Vector<Real> a_ion_densities,
							const Vector<Real> a_ion_velocities,
							const Vector<Real> a_Photon_fluxes,
							const RealVect     a_E,
							const RealVect     a_pos,
							const RealVect     a_normal,
							const Real         a_time) const{
  Vector<Real> fluxes(m_numCdrSpecies);

  // Set everything to outflow
  for (int i = 0; i < m_numCdrSpecies; i++){
    fluxes[i] = Max(0., a_extrapolated_fluxes[i]);
  }

  // For Electrons, we add ion bombardment of positive ions and the photoelectric effect
  fluxes[m_nelec_idx] = 0.;
  fluxes[m_nelec_idx] += -Max(0., a_extrapolated_fluxes[m_nplus_idx])*m_townsend2_conductor;

  // Photoelectric effect
  fluxes[m_nelec_idx] += -a_Photon_fluxes[m_Photon1_idx]*m_electrode_yield;
  fluxes[m_nelec_idx] += -a_Photon_fluxes[m_Photon2_idx]*m_electrode_yield;
  fluxes[m_nelec_idx] += -a_Photon_fluxes[m_Photon3_idx]*m_electrode_yield;

  return fluxes;
}

Vector<Real> CdrPlasmaMorrowBourdon::computeAnodeFlux(const Vector<Real> a_extrapolated_fluxes,
						      const Vector<Real> a_ion_densities,
						      const Vector<Real> a_ion_velocities,
						      const Vector<Real> a_Photon_fluxes,
						      const RealVect     a_E,
						      const RealVect     a_pos,
						      const RealVect     a_normal,
						      const Real         a_time) const{
  Vector<Real> fluxes(m_numCdrSpecies);

  // Set to outflux
  for (int i = 0; i < m_numCdrSpecies; i++){
    fluxes[i] = Max(0., a_extrapolated_fluxes[i]);
  }
  fluxes[m_nplus_idx] = a_extrapolated_fluxes[m_nplus_idx];

  return fluxes;
}

Vector<Real> CdrPlasmaMorrowBourdon::computeCdrDomainFluxes(const Real           a_time,
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

  int idx;
  int sgn;
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
      fluxes[i] = sgn*Max(0.0, sgn*a_extrap_cdr_fluxes[i]);
    }
  }
  else if(m_wallBc[idx] == 1){ // wall
    for (int i = 0; i < fluxes.size(); i++){
      fluxes[i] = 0.0;
    }
  }
  else{
    MayDay::Abort("CdrPlasmaMorrowBourdon::computeCdrDomainFluxes - uknown domain bc requested");
  }

  
  return fluxes;
}

Real CdrPlasmaMorrowBourdon::initialSigma(const Real a_time, const RealVect a_pos) const{
  return 0.;
}

CdrPlasmaMorrowBourdon::Electron::Electron(){
  m_name      = "Electron";
  m_chargeNumber    = -1;
  m_isDiffusive = true;
  m_isMobile    = true;


  ParmParse pp("CdrPlasmaMorrowBourdon");
  std::string str = "true";
  pp.get("uniform_density",    m_uniform_density);
  pp.get("seed_density",       m_seed_density);
  pp.get("seed_radius",        m_seed_radius);
  pp.get("Electron_diffusion", str);
  if(str == "true"){
    m_isDiffusive = true;
  }
  else if(str == "false"){
    m_isDiffusive = false;
  }
  Vector<Real> pos(SpaceDim);
  pp.getarr("seed_position", pos, 0, SpaceDim);
  m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
}

CdrPlasmaMorrowBourdon::Electron::~Electron(){
}

Real CdrPlasmaMorrowBourdon::Electron::initialData(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);

  return seed + m_uniform_density;
}

CdrPlasmaMorrowBourdon::PositiveSpecies::PositiveSpecies(){
  m_name      = "PositiveSpecies";
  m_chargeNumber    = 1;
  m_isDiffusive = false;
  m_isMobile    = true;

  Vector<Real> pos(SpaceDim);
  std::string str;
  
  ParmParse pp("CdrPlasmaMorrowBourdon");
  pp.get("uniform_density",  m_uniform_density);
  pp.get("seed_density",     m_seed_density);
  pp.get("seed_radius",      m_seed_radius);
  pp.get("mobile_ions", str); m_isMobile = (str == "true") ? true : false;
  pp.getarr("seed_position", pos, 0, SpaceDim);  m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
}

CdrPlasmaMorrowBourdon::PositiveSpecies::~PositiveSpecies(){
}

Real CdrPlasmaMorrowBourdon::PositiveSpecies::initialData(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);
  
  return seed + m_uniform_density;
}


CdrPlasmaMorrowBourdon::negative_species::negative_species(){
  m_name      = "negative_species";
  m_chargeNumber    = -1;
  m_isDiffusive = false;
  m_isMobile    = true;

  ParmParse pp("CdrPlasmaMorrowBourdon");
    
  // Turn off ion mobility
  std::string str;
  pp.get("mobile_ions", str); m_isMobile = (str == "true") ? true : false;
}

CdrPlasmaMorrowBourdon::negative_species::~negative_species(){
}

Real CdrPlasmaMorrowBourdon::negative_species::initialData(const RealVect a_pos, const Real a_time) const {
  return 0.;
}


CdrPlasmaMorrowBourdon::PhotonOne::PhotonOne(){
  m_name     = "PhotonOne";

  Real O2_frac, pressure;
  ParmParse pp("CdrPlasmaMorrowBourdon");
  pp.get("Photon1_A_coeff",      m_A);
  pp.get("Photon1_lambda_coeff", m_lambda);
  pp.get("gas_O2_frac",  O2_frac);
  pp.get("gas_pressure", pressure);
  
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
}

CdrPlasmaMorrowBourdon::PhotonOne::~PhotonOne(){
  
}

Real CdrPlasmaMorrowBourdon::PhotonOne::getAbsorptionCoefficient(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.
}

CdrPlasmaMorrowBourdon::PhotonTwo::PhotonTwo(){
  m_name     = "PhotonTwo";

  Real O2_frac, pressure;
  ParmParse pp("CdrPlasmaMorrowBourdon");
  pp.get("Photon2_A_coeff",      m_A);
  pp.get("Photon2_lambda_coeff", m_lambda);
  pp.get("gas_O2_frac",  O2_frac);
  pp.get("gas_pressure", pressure);
  
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
}

CdrPlasmaMorrowBourdon::PhotonTwo::~PhotonTwo(){
}

Real CdrPlasmaMorrowBourdon::PhotonTwo::getAbsorptionCoefficient(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.
}

CdrPlasmaMorrowBourdon::PhotonThree::PhotonThree(){
  m_name     = "PhotonThree";

  Real O2_frac, pressure;
  ParmParse pp("CdrPlasmaMorrowBourdon");
  pp.get("Photon3_A_coeff",      m_A);
  pp.get("Photon3_lambda_coeff", m_lambda);
  pp.get("gas_O2_frac",  O2_frac);
  pp.get("gas_pressure", pressure);

  m_pO2 = pressure*O2_frac*Units::atm2pascal;  
}

CdrPlasmaMorrowBourdon::PhotonThree::~PhotonThree(){
}

Real CdrPlasmaMorrowBourdon::PhotonThree::getAbsorptionCoefficient(const RealVect a_pos) const {
  return m_lambda*m_pO2/sqrt(3.0); // I think this is correct.

}

#include <CD_NamespaceFooter.H>
