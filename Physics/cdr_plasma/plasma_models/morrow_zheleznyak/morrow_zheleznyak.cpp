/*!
  @file   morrow_zheleznyak.cpp
  @brief  Implementation of morrow_zheleznyak.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Really, really need to revise these functions since we've changed the scaling for the rte equations
*/

#include "morrow_zheleznyak.H"
#include <CD_Units.H>
#include <CD_DataOps.H>

#include <ParmParse.H>
#include <PolyGeom.H>

#include <chrono>

#include "CD_NamespaceHeader.H"
using namespace physics::cdr_plasma;

morrow_zheleznyak::morrow_zheleznyak(){
  m_num_CdrSpecies = 3;
  m_num_RtSpecies = 1;

  m_CdrSpecies.resize(m_num_CdrSpecies);
  m_RtSpecies.resize(m_num_RtSpecies);

  m_nelec_idx   = 0;
  m_nplus_idx   = 1;
  m_nminu_idx   = 2;
  m_Photon1_idx = 0;

  // Instantiate species
  m_CdrSpecies[m_nelec_idx]    = RefCountedPtr<CdrSpecies> (new morrow_zheleznyak::electron());
  m_CdrSpecies[m_nplus_idx]    = RefCountedPtr<CdrSpecies> (new morrow_zheleznyak::positive_species());
  m_CdrSpecies[m_nminu_idx]    = RefCountedPtr<CdrSpecies> (new morrow_zheleznyak::negative_species());
  m_RtSpecies[m_Photon1_idx]  = RefCountedPtr<RtSpecies> (new morrow_zheleznyak::uv_Photon());

  // Parse some basic settings
  parse_gas_params();
  parse_see();
  parseDomainBc();
  parse_ebbc();
  parse_reaction_settings();
  parse_alpha_corr();

  // Init RNG
  if(m_seed < 0) {
    m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  m_udist11 = new std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_rng     = new std::mt19937_64(m_seed);

  // Parse initial data
  parse_initial_particles();

  // Convert to correct units and compute necessary things
  m_p  *= Units::atm2pascal;
  m_pq *= Units::atm2pascal;
  m_N   = m_p*Units::Na/(m_T*Units::R);
}

morrow_zheleznyak::~morrow_zheleznyak(){


}

void morrow_zheleznyak::advance_reaction_network(Vector<Real>&          a_particle_sources,
						 Vector<Real>&          a_Photon_sources,
						 const Vector<Real>     a_particle_densities,
						 const Vector<RealVect> a_particle_gradients,
						 const Vector<Real>     a_Photon_densities,
						 const RealVect         a_E,
						 const RealVect         a_pos,
						 const Real             a_dx,
						 const Real             a_dt,
						 const Real             a_time,
						 const Real             a_kappa) const{
  for (int i = 0; i < a_particle_densities.size(); i++){
    if(a_particle_densities[i] < 0.0) MayDay::Abort("morrow_zheleznyak::advance_network - shouldn't happen");
  }
  // Six reactions for this plasma model:
  // ===================================
  // 1. e + M   => 2e + M+ 
  // 2. e + M   => M-
  // 3. e + M+  => 0
  // 4. M- + M+ => 0
  // 5  y + M   => e + M+
  // 6. e + M   => e + M + y
  //
  // The various forms of the chemistry step is given in the routines below
  if(m_scomp == source_comp::ssa){ // SSA algorithm
    network_ssa(a_particle_sources, a_Photon_sources, a_particle_densities, a_particle_gradients, a_Photon_densities,
		a_E, a_pos, a_dx, a_dt, a_time, a_kappa);
  }
  else if(m_scomp == source_comp::tau){ // Tau leaping
    network_tau(a_particle_sources, a_Photon_sources, a_particle_densities, a_particle_gradients, a_Photon_densities,
		a_E, a_pos, a_dx, a_dt, a_time, a_kappa);
  }
  else if(m_scomp == source_comp::rre){ // Reaction rate equation
    network_rre(a_particle_sources, a_Photon_sources, a_particle_densities, a_particle_gradients, a_Photon_densities,
		a_E, a_pos, a_dx, a_dt, a_time, a_kappa);
  }

  return;
}

Real morrow_zheleznyak::excitation_rates(const Real a_E) const{
  const Real Etd = a_E/(m_N*Units::Td);

  Real y = 1.0;
  if(Etd > 100){
    y = 0.1*exp(233/Etd);
  }

  return y;
}

Real morrow_zheleznyak::sergey_factor(const Real a_O2frac) const{
  return 3E-2 + 0.4*pow(a_O2frac, 0.6);
}



// Deterministic reaction-rate equation
void morrow_zheleznyak::network_rre(Vector<Real>&          a_particle_sources,
				    Vector<Real>&          a_Photon_sources,
				    const Vector<Real>     a_particle_densities,
				    const Vector<RealVect> a_particle_gradients,
				    const Vector<Real>     a_Photon_densities,
				    const RealVect         a_E,
				    const RealVect         a_pos,
				    const Real             a_dx,
				    const Real             a_dt,
				    const Real             a_time,
				    const Real             a_kappa) const{
  const Real volume = pow(a_dx, SpaceDim);

  Vector<Real> rest(m_num_CdrSpecies, 0.0);
  Vector<int> particle_numbers(m_num_CdrSpecies, 0);
  Vector<int> Photon_numbers(m_num_RtSpecies, 0);

  // Get some aux stuff for propensity functions
  const RealVect Ve = compute_ve(a_E);
  const Real ve     = Ve.vectorLength();
  const Real E      = a_E.vectorLength();
  const Real alpha  = compute_alpha(a_E);
  const Real eta    = compute_eta(a_E);
  const Real beta   = compute_beta(a_E);

  const Real Xe = a_particle_densities[m_nelec_idx];
  const Real Xp = a_particle_densities[m_nplus_idx];
  const Real Xm = a_particle_densities[m_nminu_idx];

  Real& Se = a_particle_sources[m_nelec_idx];
  Real& Sp = a_particle_sources[m_nplus_idx];
  Real& Sm = a_particle_sources[m_nminu_idx];

  Se = 0.0;
  Sp = 0.0;
  Sm = 0.0;

  Real fcorr = 1.0;
  if(m_alpha_corr){
    const Real De        = compute_De(a_E);
    const Real ve        = compute_ve(a_E).vectorLength();;
    const RealVect gNe  = a_particle_gradients[m_nelec_idx];
    const RealVect Eunit = a_E/a_E.vectorLength();

    fcorr = 1.0 + PolyGeom::dot(Eunit, De*gNe)/((1.0 + Xe)*ve);
    fcorr = Min(fcorr, 1.0);
    fcorr = Max(0.0, fcorr);
  }
  
  // Reaction 1: e + M => 2e + M+
  const Real a1 = Xe*alpha*fcorr*ve;
  Se += a1;
  Sp += a1;
  
  // Reaction 2: e + M   => M-
  const Real a2 = Xe*eta*fcorr*ve;
  Se -= a2;
  Sm += a2;

  // Reaction 3: e + M+  => 0
  const Real a3 = Xe*Xp*beta;
  Se -= a3;
  Sp -= a3;

  // Reaction 4: M+ + M-  => 0
  const Real a4 = Xp*Xm*beta;
  Sp -= a4;
  Sm -= a4;

  // Reaction 5: Y + M => e + M+
  const Real a5 = a_Photon_densities[0]/a_dt;
  Se += a5;
  Sp += a5;

  // Reaction 6: e + M => e + M + y
  //  const Real a6 = a1*m_exc_eff*m_pq/(m_p+m_pq);
  const Real a6 = a1*volume*(m_pq/(m_p+m_pq))*sergey_factor(m_fracO2)*excitation_rates(E);
  const int S6  = poisson_reaction(a6, a_dt);
  a_Photon_sources[0] = 1.0*S6;


  return;
}

void morrow_zheleznyak::network_tau(Vector<Real>&          a_particle_sources,
				    Vector<Real>&          a_Photon_sources,
				    const Vector<Real>     a_particle_densities,
				    const Vector<RealVect> a_particle_gradients,
				    const Vector<Real>     a_Photon_densities,
				    const RealVect         a_E,
				    const RealVect         a_pos,
				    const Real             a_dx,
				    const Real             a_dt,
				    const Real             a_time,
				    const Real             a_kappa) const{
  const Real volume = pow(a_dx, SpaceDim);

  Vector<int> particle_numbers(m_num_CdrSpecies, 0);
  Vector<int> Photon_numbers(m_num_RtSpecies, 0);

  Vector<Real> x(m_num_CdrSpecies, 0.0);
  Vector<int> X(m_num_CdrSpecies, 0);
  Vector<int> Y(m_num_RtSpecies, 0);

  const Real thresh = 1.E-2;
  
  for (int i = 0; i < m_num_CdrSpecies; i++){
    X[i] = floor(thresh + a_particle_densities[i]*volume); // Integer particles
    x[i] = a_particle_densities[i] - X[i]*volume;          // "Partial particles"
  }

  for (int i = 0; i < m_num_RtSpecies; i++){
    Y[i] = floor(thresh + a_Photon_densities[i]);
  }

  int se = 0;
  int sp = 0;
  int sm = 0;

  // Get some aux stuff for propensity functions
  const RealVect Ve = compute_ve(a_E);
  const Real ve     = Ve.vectorLength();
  const Real E      = a_E.vectorLength();
  const Real alpha  = compute_alpha(a_E);
  const Real eta    = compute_eta(a_E);
  const Real beta   = compute_beta(a_E);

  const Real xe = x[m_nelec_idx];
  const Real xp = x[m_nplus_idx];
  const Real xm = x[m_nminu_idx];
  
  const Real Xe = X[m_nelec_idx];
  const Real Xp = X[m_nplus_idx];
  const Real Xm = X[m_nminu_idx];

  // Reaction 1: e + M => 2e + M+
  const Real A1 = Xe*alpha*ve;
  const Real a1 = xe*alpha*ve;
  const int S1  = poisson_reaction(A1, a_dt);
  se += S1;
  sp += S1;
  
  // Reaction 2: e + M   => M-
  const Real A2 = Xe*eta*ve;
  const Real a2 = xe*eta*ve;
  const int S2  = poisson_reaction(A2, a_dt);
  se -= S2;
  sm += S2;

  // Reaction 3: e + M+  => 0
  const Real A3 = Xe*Xp*beta/volume;
  const Real a3 = xe*xp*beta;
  const int S3  = poisson_reaction(A3, a_dt);
  // se -= S3;
  // sp -= S3;

  // Reaction 4: M+ + M-  => 0
  const Real A4 = Xp*Xm*beta/volume;
  const Real a4 = xp*xm*beta;
  const int  S4 = poisson_reaction(A4, a_dt);
  // sp -= S4;
  // sm -= S4;

  // Reaction 5: Y + M => e + M+
  se += Y[0];
  sp += Y[0];

  // Reaction 6: e + M => e + M + y
  const Real quench = m_pq/(m_p+m_pq);
  const Real a6     = Xe*alpha*ve*sergey_factor(m_fracO2)*excitation_rates(E);
  const int S6      = poisson_reaction(a6, a_dt);
  a_Photon_sources[0] = 1.0*S6;

  // Do some scaling
  const Real factor = 1./(volume*a_dt);
  Real& Se = a_particle_sources[m_nelec_idx];
  Real& Sp = a_particle_sources[m_nplus_idx];
  Real& Sm = a_particle_sources[m_nminu_idx];

  Se = se*factor + (a1 - a2 - a3);
  Sp = sp*factor + (a1 - a3 - a4);
  Sm = sm*factor + (a2 - a4);

#if 1 // Debug
  const int res = -se + sp - sm;
  if(res != 0) MayDay::Abort("nope");
#endif
  

  return;
}

// Tau leaping method
void morrow_zheleznyak::network_ssa(Vector<Real>&          a_particle_sources,
				    Vector<Real>&          a_Photon_sources,
				    const Vector<Real>     a_particle_densities,
				    const Vector<RealVect> a_particle_gradients,
				    const Vector<Real>     a_Photon_densities,
				    const RealVect         a_E,
				    const RealVect         a_pos,
				    const Real             a_dx,
				    const Real             a_dt,
				    const Real             a_time,
				    const Real             a_kappa) const{

  const Real volume = pow(a_dx, 3);

  Vector<Real> x(m_num_CdrSpecies, 0.0);
  Vector<int>  X(m_num_CdrSpecies, 0);
  Vector<int>  Y(m_num_RtSpecies, 0);

  //
  const RealVect Ve = compute_ve(a_E);
  const Real ve     = Ve.vectorLength();
  const Real E      = a_E.vectorLength();
  const Real alpha  = compute_alpha(a_E);
  const Real eta    = compute_eta(a_E);
  const Real beta   = compute_beta(a_E);
  const Real eff    = sergey_factor(m_fracO2)*excitation_rates(E)*m_pq/(m_p+m_pq);


  // Six reactions for this plasma model:
  // ===================================
  // 1. e + M   => 2e + M+ 
  // 2. e + M   => M-
  // 3. e + M+  => 0
  // 4. M- + M+ => 0
  // 5. e + M   => e + M + y
  // 6  y + M   => e + M+
  const Real thresh = 1.E-2;
  
  // Initial particle densities
  for (int i = 0; i < m_num_CdrSpecies; i++){
    x[i] = a_particle_densities[i] - X[i]*volume;
    X[i] = floor(thresh + a_particle_densities[i]*volume);
  }
  const Vector<int> X0 = X;


  // Initial Photon densities
  for (int i = 0; i < m_num_RtSpecies; i++){
    Y[i] = floor(thresh + a_Photon_densities[i]);
  }

  // Deposit Photons. This is reaction #6 in the list above. We will only do the other 5 reactions.
  // X[m_nelec_idx] += Y[0];
  // X[m_nplus_idx] += Y[0];
  Y[0] = 0;

  const Real xe = x[m_nelec_idx];
  const Real xp = x[m_nplus_idx];
  const Real xm = x[m_nminu_idx];
  
  int& Xe = X[m_nelec_idx];
  int& Xp = X[m_nplus_idx];
  int& Xm = X[m_nminu_idx];

  // Do SSA algorithm for INTEGER particles
  Real tau  = 0.0;
  while (tau <= a_dt){

    Vector<Real> ar(5);
    
    ar[0] = Xe*alpha*ve;        // Impact ionization
    ar[1] = Xe*eta*ve;          // Electron attachment
#if 0
    ar[2] = Xe*Xp*beta/volume;  // Electron-ion recombination
    ar[3] = Xp*Xm*beta/volume;  // Ion-ion recombination
#else
    ar[2] = 0.0;
    ar[3] = 0.0;
#endif
    ar[4] = ar[0]*eff;          // Photon generation

    // Total reaction rate
    Real A = 0;
    for (int i=0; i<ar.size(); i++){
      A += ar[i];
    }

    if(A > 0){

      // Normalize ar
      //      ar *= 1./A;
      
      const Real u1 = (*m_udist01)(*m_rng);
      const Real u2 = (*m_udist01)(*m_rng);
      const Real dtau = (1./A)*log(1/u1);

      // Type of reaction
      int r = 0;
      Real sum = 0.0;
      for (int i=0; i < ar.size(); i++){
	sum += ar[i];
	if(u2*A < sum) {
	  r = i;
	  break;
	}
      }

      // Advance reactions
      if(r == 0){ // Impact ionization
	Xe += 1;
	Xp += 1;
      }
      else if(r == 1){ // Attachment
	Xe -= 1;
	Xm += 1;
      }
      else if(r == 2){ // Electron-ion recombination
	Xe -= 1;
	Xp -= 1;
      }
      else if(r == 3){ // Ion-ion recombination
	Xp -=1;
	Xe -=1;
      }
      else { // Send out a Photon
	Y[0] += 1;
      }

      // Increment with time 
      tau += dtau;
    }
    else{
      break;
    }
  }
  
  // Now we need to normalize some stuff
  a_Photon_sources[0] = 1.0*Y[0];

  // Non-integer reactions
  const Real a1 = xe*alpha*ve;
  const Real a2 = xe*eta*ve;
  const Real a3 = xe*xp*beta;
  const Real a4 = xp*xm*beta;
  
  const Real factor = 1./(volume*a_dt);
  for (int i = 0; i < m_num_CdrSpecies; i++){
    a_particle_sources[i] = factor*(X[i] - X0[i]);
  }

  Real& Se = a_particle_sources[m_nelec_idx];
  Real& Sp = a_particle_sources[m_nplus_idx];
  Real& Sm = a_particle_sources[m_nminu_idx];

  Se += (a1 - a2 - a3);
  Sp += (a1 - a3 - a4);
  Sm += (a2 - a4);
  

  return;
}

Vector<RealVect> morrow_zheleznyak::compute_cdr_velocities(const Real         a_time,
							   const RealVect     a_pos,
							   const RealVect     a_E,
							   const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> velocities(m_num_CdrSpecies, RealVect::Zero);
  
  velocities[m_nelec_idx] = this->compute_ve(a_E);
  velocities[m_nplus_idx] = this->compute_vp(a_E);
  velocities[m_nminu_idx] = this->compute_vn(a_E);

  return velocities;
}

RealVect morrow_zheleznyak::compute_ve(const RealVect a_E) const{
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

RealVect morrow_zheleznyak::compute_vp(const RealVect a_E) const{
  const RealVect E = a_E*1.E-2;           // E in V/cm
  RealVect vp = 2.34*E*m_p/Units::atm2pascal;  // Morrow-Lowke wants V/cm
  vp *= 0.01;                             // Morrow-Lowke expression is in cm/s
  return vp;  
}

RealVect morrow_zheleznyak::compute_vn(const RealVect a_E) const{
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

Real morrow_zheleznyak::compute_alpha(const RealVect a_E) const{
  Real alpha    = 0.;
  Real alphabyN = 0.;

  const RealVect E = a_E*1.E-2;        // Morrow-Lowke wants E in V/cm
  const Real Emag  = E.vectorLength(); // 
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

Real morrow_zheleznyak::compute_eta(const RealVect a_E) const{

  const Real eta2 = this->compute_eta2(a_E); 
  const Real eta3 = this->compute_eta3(a_E);
  const Real eta  = eta2 + eta3;

  return eta;
}

Real morrow_zheleznyak::compute_eta2(const RealVect a_E) const{
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
    eta2byN = 6.089E-4*EbyN + 2.893E-19;
  }

  eta2byN *= 1.E-4;   // Morrow-Lowke expression is in cm^2, make it m^2
  eta2 = eta2byN*m_N; //

  return eta2;
}

Real morrow_zheleznyak::compute_eta3(const RealVect a_E) const{
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

Real morrow_zheleznyak::compute_beta(const RealVect a_E) const{
  Real beta = 2.0E-7;
  beta *= 1.E-6; // Morrow-Lowke expression is in cm^3. Make it m^3
  return beta;
}

Real morrow_zheleznyak::compute_De(const RealVect a_E) const{
  const RealVect E  = a_E*1.E-2;                 // Morrow-Lowke wants E in V/cm
  const Real Emag   = E.vectorLength();          //
  const Real N      = m_N*1.E-6;                 // Morrow-Lowke weants N in cm^3
  const Real EbyN   = Emag/N;                    //

  //
  const RealVect Ve = this->compute_ve(a_E);     // Does it's own conversion, comes out in m/s (aka. mentally sane units)
  const Real ve     = Ve.vectorLength()*1.E-2;   // Make it cm/s
  Real De = 0.3341E9*pow(EbyN, 0.54069)*ve/Emag;
  
  De *= 1.E-4; // Morrow-Lowke expression is in cm^2/s. Make it m^2/s

  return De;
}


Vector<Real> morrow_zheleznyak::compute_cdr_diffusion_coefficients(const Real         a_time,
								   const RealVect     a_pos,
								   const RealVect     a_E,
								   const Vector<Real> a_cdr_densities) const{
  Vector<Real> diffCo(m_num_CdrSpecies, 0.0);
  diffCo[m_nelec_idx] = this->compute_De(a_E);
  diffCo[m_nplus_idx] = 0.;
  diffCo[m_nminu_idx] = 0.;
  
  return diffCo;
}

Vector<Real> morrow_zheleznyak::compute_cdr_fluxes(const Real         a_time,
						   const RealVect     a_pos,
						   const RealVect     a_normal,
						   const RealVect     a_E,
						   const Vector<Real> a_cdr_densities,
						   const Vector<Real> a_cdr_velocities,
						   const Vector<Real> a_cdr_gradients,
						   const Vector<Real> a_rte_fluxes,
						   const Vector<Real> a_extrap_cdr_fluxes,
						   const Real         a_townsend2,
						   const Real         a_quantum_efficiency) const {
  Vector<Real> fluxes(m_num_CdrSpecies, 0.0);
  
  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;

  // Switch for setting drift flux to zero for charge species
  Vector<Real> aj(m_num_CdrSpecies, 0.0);
  for (int i = 0; i < m_num_CdrSpecies; i++){
    if(DataOps::sgn(m_CdrSpecies[i]->getChargeNumber())*PolyGeom::dot(a_E, a_normal) < 0){
      aj[i] = 1.0;
    }
    else {
      aj[i] = 0.0;
    }
  }

  // Drift outflow for now
  for (int i = 0; i < m_num_CdrSpecies; i++){
    fluxes[i] = aj[i]*(a_extrap_cdr_fluxes[i]);
  }
  
  return fluxes;
}

Vector<Real> morrow_zheleznyak::compute_cdr_domain_fluxes(const Real           a_time,
							  const RealVect       a_pos,
							  const int            a_dir,
							  const Side::LoHiSide a_side,
							  const RealVect       a_E,
							  const Vector<Real>   a_cdr_densities,
							  const Vector<Real>   a_cdr_velocities,
							  const Vector<Real>   a_cdr_gradients,
							  const Vector<Real>   a_rte_fluxes,
							  const Vector<Real>   a_extrap_cdr_fluxes) const{
  Vector<Real> fluxes(m_num_CdrSpecies, 0.0); 

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
      fluxes[i] = a_extrap_cdr_fluxes[i];
    }
  }
  else if(m_wallBc[idx] == 1){ // wall
    for (int i = 0; i < fluxes.size(); i++){
      fluxes[i] = 0.0;
    }
  }
  else{
    MayDay::Abort("morrow_zheleznyak::compute_cdr_domain_fluxes - uknown domain bc requested");
  }

  
  return fluxes;
}

Vector<Real> morrow_zheleznyak::compute_cdr_electrode_fluxes(const Real         a_time,
							     const RealVect     a_pos,
							     const RealVect     a_normal,
							     const RealVect     a_E,
							     const Vector<Real> a_cdr_densities,
							     const Vector<Real> a_cdr_velocities,
							     const Vector<Real> a_cdr_gradients,
							     const Vector<Real> a_rte_fluxes,
							     const Vector<Real> a_extrap_cdr_fluxes) const {
  if(m_extrap_electrode_ebbc){
    return a_extrap_cdr_fluxes;
  }
  else{
    return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients,
				    a_rte_fluxes, a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
  }
}

Vector<Real> morrow_zheleznyak::compute_cdr_dielectric_fluxes(const Real         a_time,
							      const RealVect     a_pos,
							      const RealVect     a_normal,
							      const RealVect     a_E,
							      const Vector<Real> a_cdr_densities,
							      const Vector<Real> a_cdr_velocities,
							      const Vector<Real> a_cdr_gradients,
							      const Vector<Real> a_rte_fluxes,
							      const Vector<Real> a_extrap_cdr_fluxes) const {
  if(m_extrap_dielectric_ebbc){
    return a_extrap_cdr_fluxes;
  }
  else{
    return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients,
				    a_rte_fluxes, a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
  }
}

int morrow_zheleznyak::poisson_reaction(const Real a_propensity, const Real a_dt) const{
  int value = 0;
  const Real mean = a_propensity*a_dt;

  if(mean < m_poiss_exp_swap){
    std::poisson_distribution<int> dist(mean);
    value = dist(*m_rng);
  }
  else{
    std::normal_distribution<double> dist(mean, sqrt(mean));
    value = dist(*m_rng);
  }

  return value;
}

Real morrow_zheleznyak::initial_sigma(const Real a_time, const RealVect a_pos) const{
  return 0.;
}

morrow_zheleznyak::electron::electron(){
  m_name      = "electron";
  m_chargeNumber    = -1;
  m_isDiffusive = true;
  m_isMobile    = true;
  m_unit      = "m-3";

  Vector<Real> pos(SpaceDim);
  std::string str;

  ParmParse pp("morrow_zheleznyak");
  pp.get("uniform_density",     m_uniform_density);
  pp.get("seed_density",        m_seed_density);
  pp.get("seed_radius",         m_seed_radius);
  pp.get("diffusive_electrons", str); m_isDiffusive = (str == "true") ? true : false;

  pp.getarr("seed_position", pos, 0, SpaceDim); m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
}

morrow_zheleznyak::positive_species::positive_species(){
  m_name      = "positive_species";
  m_chargeNumber    = 1;
  m_isDiffusive = false;
  m_isMobile    = false;
  m_unit      = "m-3";

  Vector<Real> pos(SpaceDim);
  std::string str;

  ParmParse pp("morrow_zheleznyak");
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density",    m_seed_density);
  pp.get("seed_radius",     m_seed_radius);
  pp.get("mobile_ions",     str); m_isMobile    = (str == "true") ? true : false;

  pp.getarr("seed_position", pos, 0, SpaceDim); m_seed_pos = RealVect(D_DECL(pos[0], pos[1], pos[2]));
}

morrow_zheleznyak::negative_species::negative_species(){
  m_name      = "negative_species";
  m_chargeNumber    = -1;
  m_isDiffusive = false;
  m_isMobile    = false;
  m_unit      = "m-3";

  std::string str;

  ParmParse pp("morrow_zheleznyak");
  pp.get("mobile_ions",    str); m_isMobile    = (str == "true") ? true : false;
}

morrow_zheleznyak::electron::~electron(){
  
}

morrow_zheleznyak::positive_species::~positive_species(){
  
}

morrow_zheleznyak::negative_species::~negative_species(){
  
}

Real morrow_zheleznyak::electron::initialData(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);

  return m_uniform_density + seed;
}

Real morrow_zheleznyak::positive_species::initialData(const RealVect a_pos, const Real a_time) const {
  const Real factor = (a_pos - m_seed_pos).vectorLength()/m_seed_radius;
  const Real seed   = m_seed_density*exp(-factor*factor);
  
  return m_uniform_density + seed;
}

Real morrow_zheleznyak::negative_species::initialData(const RealVect a_pos, const Real a_time) const {
  return 0.;
}

morrow_zheleznyak::uv_Photon::uv_Photon(){
  m_name   = "uv_Photon";

  Real pressure, O2_frac;
  
  ParmParse pp("morrow_zheleznyak");
  pp.get("photoi_f1",      m_f1);
  pp.get("photoi_f2",      m_f2);
  pp.get("photoi_K1",      m_K1);
  pp.get("photoi_K2",      m_K2);

  pp.get("gas_O2_frac",  O2_frac);
  pp.get("gas_pressure", pressure);
  pp.get("seed",         m_seed);

  // Convert units
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
  m_K1  = m_K1*m_pO2;
  m_K2  = m_K2*m_pO2;

  // Seed the RNG
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
}

morrow_zheleznyak::uv_Photon::~uv_Photon(){
  
}

Real morrow_zheleznyak::uv_Photon::getKappa(const RealVect a_pos) const {
  MayDay::Abort("morrow_zheleznyak::uv_Photon::getKappa - should not be called. morrow_zheleznyak is used with the McPhoto module");
}

Real morrow_zheleznyak::uv_Photon::get_random_kappa() const {
  const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
  return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
}

void morrow_zheleznyak::parse_gas_params(){
  ParmParse pp("morrow_zheleznyak");
  std::string str;
  pp.get("gas_temperature",            m_T);
  pp.get("gas_N2_frac",                m_fracN2);
  pp.get("gas_O2_frac",                m_fracO2);
  pp.get("gas_pressure",               m_p);
  pp.get("gas_quenching_pressure",     m_pq);
}

void morrow_zheleznyak::parse_alpha_corr(){
  ParmParse pp("morrow_zheleznyak");

  std::string str;
  pp.get("use_alpha_correction", str);

  m_alpha_corr = (str == "true") ? true : false;
}

void morrow_zheleznyak::parse_see(){

  ParmParse pp("morrow_zheleznyak");

  pp.get("electrode_townsend2",           m_townsend2_electrode);
  pp.get("dielectric_townsend2",          m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
}

void morrow_zheleznyak::parseDomainBc(){

  ParmParse pp("morrow_zheleznyak");
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

void morrow_zheleznyak::parse_ebbc(){

  ParmParse pp("morrow_zheleznyak");

  std::string str;

  pp.get("extrap_electrode", str);  m_extrap_electrode_ebbc  = (str == "true") ? true : false;
  pp.get("extrap_dielectric", str); m_extrap_dielectric_ebbc = (str == "true") ? true : false;
}

  

void morrow_zheleznyak::parse_reaction_settings(){
  ParmParse pp("morrow_zheleznyak");
  std::string str;

  pp.get("chemistry", str);
  if(str == "ssa"){
    m_scomp = source_comp::ssa;
  }
  else if(str == "tau"){
    m_scomp = source_comp::tau;
  }
  else if(str == "rre"){
    m_scomp = source_comp::rre;
  }
  else{
    MayDay::Abort("morrow_zheleznyak::parse_reaction_settings - stop!");
  }

  pp.get("seed",           m_seed);
  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  pp.get("cutoff_poisson", m_cutoff_poisson);
}

void morrow_zheleznyak::parse_initial_particles(){

  List<Particle> p;

  // Get types of particles
  add_uniform_particles(p);
  add_gaussian_particles(p);
  
  // Copy initial particles to various species
  m_CdrSpecies[m_nelec_idx]->getInitialParticles() = p;
  m_CdrSpecies[m_nplus_idx]->getInitialParticles() = p;

  // Get the initial deposition scheme
  ParmParse pp("morrow_zheleznyak");
  std::string str;
  InterpType deposition;

  pp.get("particle_deposition", str);
  if(str == "cic"){
    deposition = InterpType::CIC;
  }
  else if(str == "ngp"){
    deposition = InterpType::NGP;
  }
  else if(str == "tsc"){
    deposition = InterpType::TSC;
  }
  else{
    MayDay::Abort("morrow_zheleznyak::parse_initial_particles - unknown deposition type requested");
  }
  
  m_CdrSpecies[m_nelec_idx]->getDeposition() = deposition;
  m_CdrSpecies[m_nplus_idx]->getDeposition() = deposition;
}

void morrow_zheleznyak::add_uniform_particles(List<Particle>& a_particles){
  
  // Get lo/hi sides
  Real weight;
  RealVect lo, hi;
  Vector<Real> vec(SpaceDim);
  {
    ParmParse pp1("physical_domain");

    pp1.getarr("lo_corner", vec, 0, SpaceDim); lo = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    pp1.getarr("hi_corner", vec, 0, SpaceDim); hi = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  }


  // Make RNGs
  auto rngX = new std::uniform_real_distribution<Real>(lo[0], hi[0]);
  auto rngY = new std::uniform_real_distribution<Real>(lo[1], hi[1]);
#if CH_SPACEDIM==3
  auto rngZ = new std::uniform_real_distribution<Real>(lo[2], hi[2]);
#endif

  int num_uniform_particles;
  Real num_particles;
  
  // Create uniform particles
  ParmParse pp("morrow_zheleznyak");
  pp.get("uniform_particles",  num_particles);
  pp.get("particle_weight", weight);
  num_uniform_particles = round(num_particles);

  for (int i = 0; i < num_uniform_particles; i++){
    const Real x = (*rngX)(*m_rng);
    const Real y = (*rngY)(*m_rng);
#if CH_SPACEDIM==3
    const Real z = (*rngZ)(*m_rng);
#endif
    RealVect pos = RealVect(D_DECL(x, y, z));
    a_particles.add(Particle(weight, pos));
  }

  delete rngX;
  delete rngY;
#if CH_SPACEDIM==3
  delete rngZ;
#endif
}

void morrow_zheleznyak::add_gaussian_particles(List<Particle>& a_particles){

  // Create uniform particles
  ParmParse pp("morrow_zheleznyak");

  int num_gaussian_particles;
  Real num_particles;
  Real weight;
  Real gaussian_radius;
  RealVect gaussian_center;
  Vector<Real> vec(SpaceDim);
  // Create Gaussian seed particles
  pp.get("particle_weight", weight);
  pp.get("gaussian_particles", num_particles);
  pp.get("gaussian_radius",    gaussian_radius);
  pp.getarr("gaussian_center", vec, 0, SpaceDim); gaussian_center = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  num_gaussian_particles = round(num_particles);

  m_gauss = std::normal_distribution<Real>(0., gaussian_radius);
  for (int i = 0; i < num_gaussian_particles; i++){
    RealVect pos = gaussian_center + randomGaussian(gaussian_radius);
    a_particles.add(Particle(weight, pos));
  }

}
 
RealVect morrow_zheleznyak::randomGaussian(const Real a_rad){

  const Real rad = m_gauss(*m_rng);
  return rad*randomDirection();
}

RealVect morrow_zheleznyak::randomDirection(){
#if CH_SPACEDIM == 2
  return randomDirection2D();
#else
  return randomDirection3D();
#endif
}

#if CH_SPACEDIM == 2
RealVect morrow_zheleznyak::randomDirection2D(){
  const Real EPS = 1.E-8;
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = (*m_udist11)(*m_rng);
    x2 = (*m_udist11)(*m_rng);
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
}
#endif

#if CH_SPACEDIM==3
RealVect morrow_zheleznyak::randomDirection3D(){
  const Real EPS = 1.E-8;
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = (*m_udist11)(*m_rng);
    x2 = (*m_udist11)(*m_rng);
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
}
#endif

#include "CD_NamespaceFooter.H"
