#include "air3_bourdon.H"
#include "air3_bourdon_species.H"
#include <CD_DataOps.H>
#include "units.H" 

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include <ParmParse.H>
#include <PolyGeom.H>

#include "CD_NamespaceHeader.H"
using namespace physics::cdr_plasma;

std::string air3_bourdon::s_bolsig_mobility = "# Electron mobility (E/N, mu*N)";
std::string air3_bourdon::s_bolsig_diffco   = "# Electron diffusion coefficient (E/N, D*N)";
std::string air3_bourdon::s_bolsig_alpha    = "# Townsend ionization coeff (E/N, alpha/N)";
std::string air3_bourdon::s_bolsig_eta      = "# Townsend attachment coeff (E/N, eta/N)";

air3_bourdon::air3_bourdon() {

  ParmParse pp("air3_bourdon");

  // Get input parameters
  pp.get("pressure",                      m_p);
  pp.get("quenching_pressure",            m_pq);
  pp.get("temperature",                   m_T);
  pp.get("transport_file",                m_transport_file);
  pp.get("uniform_tables",                m_uniform_entries);
  pp.get("use_alpha_corr",                m_alpha_corr);      
  pp.get("mobile_electrons",              m_isMobile_electrons);    
  pp.get("diffusive_electrons",           m_isDiffusive_electrons);
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
  m_p     *= units::s_atm2pascal;
  m_pq    *= units::s_atm2pascal;
  m_N      = m_p*units::s_Na/(m_T*units::s_R);
  
  m_ion_diffusion = m_ion_mobility*(units::s_kb*m_T)/units::s_Qe; // Einstein relation. 


  // Read from input file and put in lookup tables. 
  std::ifstream infile(m_transport_file);
  if(infile.good()){
    infile.close();

    read_file_entries(m_e_mobility, air3_bourdon::s_bolsig_mobility);
    read_file_entries(m_e_diffco,   air3_bourdon::s_bolsig_diffco);
    read_file_entries(m_e_alpha,    air3_bourdon::s_bolsig_alpha);
    read_file_entries(m_e_eta,      air3_bourdon::s_bolsig_eta);
    
    m_e_mobility.scaleX(m_N*units::s_Td);
    m_e_diffco.scaleX(m_N*units::s_Td);
    m_e_alpha.scaleX(m_N*units::s_Td);
    m_e_eta.scaleX(m_N*units::s_Td);
    
    m_e_mobility.scaleY(1./m_N); 
    m_e_diffco.scaleY(1./m_N); 
    m_e_alpha.scaleY(m_N); 
    m_e_eta.scaleY(m_N);

    m_e_diffco.makeUniform(m_uniform_entries);
    m_e_mobility.makeUniform(m_uniform_entries);
    m_e_alpha.makeUniform(m_uniform_entries);
    m_e_eta.makeUniform(m_uniform_entries);

  }
  else{
    MayDay::Abort("air3_bourdon::parse_transport_file - could not find transport data");
  }


  // Parse domain boundary conditions. 
  parseDomainBc();

  // Instantiate species. 
  instantiate_species();
}

air3_bourdon::~air3_bourdon() {

}

void air3_bourdon::read_file_entries(LookupTable& a_table, const std::string a_string){
  Real x, y;
  bool read_line = false;
  std::ifstream infile(m_transport_file);
  std::string line;

  while (std::getline(infile, line)){

    // Right trim string
    line.erase(line.find_last_not_of(" \n\r\t")+1);

    if(line == a_string){ // Begin reading
      read_line = true;
    }
    else if(line == "" & read_line){ // Stop reading
      read_line = false;
    }

    if(read_line){
      std::istringstream iss(line);
      if (!(iss >> x >> y)) {
	continue;
      }
      a_table.addEntry(x, y);
    }
  }
  infile.close();
}

void air3_bourdon::instantiate_species(){
  m_num_CdrSpecies = 3;
  m_num_RtSpecies = 3;

  m_elec_idx = 0;
  m_plus_idx = 1;
  m_minu_idx = 2;
  m_pho1_idx = 0;
  m_pho2_idx = 1;
  m_pho3_idx = 2;

  m_CdrSpecies.resize(m_num_CdrSpecies);
  m_CdrSpecies[m_elec_idx]  = RefCountedPtr<CdrSpecies> (new air3_bourdon::electron());
  m_CdrSpecies[m_plus_idx]  = RefCountedPtr<CdrSpecies> (new air3_bourdon::M_plus());
  m_CdrSpecies[m_minu_idx]  = RefCountedPtr<CdrSpecies> (new air3_bourdon::M_minus());

  m_RtSpecies.resize(m_num_RtSpecies);
  m_RtSpecies[m_pho1_idx] = RefCountedPtr<RtSpecies> (new air3_bourdon::Photon_one());
  m_RtSpecies[m_pho2_idx] = RefCountedPtr<RtSpecies> (new air3_bourdon::Photon_two());
  m_RtSpecies[m_pho3_idx] = RefCountedPtr<RtSpecies> (new air3_bourdon::Photon_three());
}

void air3_bourdon::parseDomainBc(){

  ParmParse pp("air3_bourdon");
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

void air3_bourdon::advance_reaction_network(Vector<Real>&          a_particle_sources,
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
  const Real ve     = E*m_e_mobility.getEntry(E);
  
  // Ionization and attachment coefficients
  Real alpha  = m_e_alpha.getEntry(E);
  Real eta    = m_e_eta.getEntry(E);

  // Modify alpha
  if(m_alpha_corr){
    const RealVect Eunit = a_E/a_E.vectorLength();
    const Real De        = m_e_diffco.getEntry(E);
    const RealVect gNe   = a_particle_gradients[m_elec_idx];

    Real fcorr = 1.0;
    fcorr = 1.0 + PolyGeom::dot(Eunit, De*gNe)/(1.0+a_particle_densities[m_elec_idx]*ve);
    fcorr = Min(fcorr, 1.0);
    fcorr = Max(0.0, fcorr);

    alpha = alpha*fcorr;
  }

  const Real R1 = alpha*ve*a_particle_densities[m_elec_idx];
  const Real R2 = eta*ve*a_particle_densities[m_elec_idx];

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

  // Photoionization, M + y => e + M+
  const air3_bourdon::Photon_one*   Photon1 = static_cast<air3_bourdon::Photon_one*>   (&(*m_RtSpecies[m_pho1_idx]));
  const air3_bourdon::Photon_two*   Photon2 = static_cast<air3_bourdon::Photon_two*>   (&(*m_RtSpecies[m_pho2_idx]));
  const air3_bourdon::Photon_three* Photon3 = static_cast<air3_bourdon::Photon_three*> (&(*m_RtSpecies[m_pho3_idx]));

  const Real Sph = m_photoi_eff*units::s_c0*m_O2frac*m_p*(Photon1->get_A()*a_Photon_densities[m_pho1_idx]
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

  return;
}

Vector<Real> air3_bourdon::compute_cdr_diffusion_coefficients(const Real         a_time,
							      const RealVect     a_pos,
							      const RealVect     a_E,
							      const Vector<Real> a_cdr_densities) const {

  Vector<Real> dco(m_num_CdrSpecies, 0.0);
  dco[m_elec_idx] = m_e_diffco.getEntry(a_E.vectorLength());
  dco[m_plus_idx] = m_ion_diffusion;
  dco[m_minu_idx] = m_ion_diffusion;
  
  return dco;
}
  
Vector<RealVect> air3_bourdon::compute_cdr_velocities(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities) const{
  Vector<RealVect> vel(m_num_CdrSpecies, RealVect::Zero);

  vel[m_elec_idx] = -a_E*m_e_mobility.getEntry(a_E.vectorLength());
  vel[m_plus_idx] =  a_E*m_ion_mobility;
  vel[m_minu_idx] = -a_E*m_ion_mobility;
  
  return vel;
}
  
Vector<Real> air3_bourdon::compute_cdr_domain_fluxes(const Real           a_time,
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
    MayDay::Abort("morrow_jiang::compute_cdr_domain_fluxes - uknown domain bc requested");
  }

  
  return fluxes;
}
  
Vector<Real> air3_bourdon::compute_cdr_electrode_fluxes(const Real         a_time,
							const RealVect     a_pos,
							const RealVect     a_normal,
							const RealVect     a_E,
							const Vector<Real> a_cdr_densities,
							const Vector<Real> a_cdr_velocities,
							const Vector<Real> a_cdr_gradients,
							const Vector<Real> a_rte_fluxes,
							const Vector<Real> a_extrap_cdr_fluxes) const{
  return compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			    a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
}

Vector<Real> air3_bourdon::compute_cdr_dielectric_fluxes(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_normal,
							 const RealVect     a_E,
							 const Vector<Real> a_cdr_densities,
							 const Vector<Real> a_cdr_velocities,
							 const Vector<Real> a_cdr_gradients,
							 const Vector<Real> a_rte_fluxes,
							 const Vector<Real> a_extrap_cdr_fluxes) const{
  return compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
			    a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

Vector<Real> air3_bourdon::compute_cdr_fluxes(const Real         a_time,
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
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  return fluxes;
}

Real air3_bourdon::initial_sigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

Real air3_bourdon::compute_alpha(const RealVect a_E) const{
  const Real E     = a_E.vectorLength();
  const Real alpha = m_e_alpha.getEntry(E);

  return alpha;
}
#include "CD_NamespaceFooter.H"
