/*!
  @file   air_bolsig.cpp
  @brief  Implementation of air_bolsig.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "air_bolsig.H"
#include "units.H"
#include "perlin_if.H"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ParmParse.H>

std::string air_bolsig::s_script_file  = "bolsig_inputs.";
std::string air_bolsig::s_data_file    = "bolsig_outputs.";
std::string air_bolsig::s_bolsig_alpha = "E/N (Td)	Townsend ioniz. coef. alpha/N (m2)";
std::string air_bolsig::s_bolsig_mob   = "E/N (Td)	Mobility *N (1/m/V/s)";
std::string air_bolsig::s_bolsig_diff  = "E/N (Td)	Diffusion coefficient *N (1/m/s)";
std::string air_bolsig::s_bolsig_eta   = "E/N (Td)	Townsend attach. coef. eta/N (m2)";
std::string air_bolsig::s_bolsig_kex23 = "C23   N2    Excitation    12.25 eV";                              
std::string air_bolsig::s_bolsig_kex24 = "C24   N2    Excitation    13.00 eV";

air_bolsig::air_bolsig(){
  m_bolsig_path  = "./";
  m_lxcat_path   = "./";
  m_lxcat_file   = "LXCat-June2013.txt";
  m_min_townsend = 0.0;
  m_max_townsend = 200.0;
  m_grid_points  = 100;
  
  m_gas_temp = 300.;
  m_frac_O2  = 0.2;
  m_frac_N2  = 0.8;
  m_p        = 1.0;
  m_pq       = 0.03947;

  m_background_rate               = 1.E9;
  m_electron_recombination        = 5.E-14;
  m_electron_detachment           = 1.E-18;
  m_ion_recombination             = 2.07E-12;
  m_positive_species_mobility     = 2.E-4;
  m_negative_species_mobility     = 2.E-4;
  m_excitation_efficiency         = 0.6;
  m_photoionization_efficiency    = 0.1;
  m_townsend2_electrode           = 1.E-4;
  m_townsend2_dielectric          = 1.E-6;
  m_electrode_quantum_efficiency  = 1.E-4;
  m_dielectric_quantum_efficiency = 1.E-6;

  m_noise_amplitude   = 0.0;
  m_noise_octaves     = 1;
  m_noise_persistence = 0.5;
  m_noise_frequency   = RealVect::Unit;
      
      
  { // Get path to BOLSIG and database file
    ParmParse pp("air_bolsig");
    pp.query("bolsig_path",  m_bolsig_path);
    pp.query("lxcat_path",   m_lxcat_path);
    pp.query("lxcat_file",   m_lxcat_file);
    pp.query("grid_points",  m_grid_points);
    pp.query("min_townsend", m_min_townsend);
    pp.query("max_townsend", m_max_townsend);
  }

  { // Get gas parameters
    ParmParse pp("air_bolsig");
    pp.query("gas_temperature",        m_gas_temp);
    pp.query("gas_O2_frac",            m_frac_O2);
    pp.query("gas_N2_frac",            m_frac_N2);
    pp.query("gas_pressure",           m_p);
    pp.query("gas_quenching_pressure", m_pq);
  }

  { // Kinetic things not covered by BOLSIG
    ParmParse pp("air_bolsig");
    pp.query("ionization_rate",               m_background_rate);
    pp.query("electron_recombination",        m_electron_recombination);
    pp.query("electron_detachment",           m_electron_detachment);
    pp.query("ion_recombination",             m_ion_recombination);
    pp.query("positive_species_mobility",     m_positive_species_mobility);
    pp.query("negative_species_mobility",     m_negative_species_mobility);
    pp.query("excitation_efficiency",         m_excitation_efficiency);
    pp.query("photoionization_efficiency",    m_photoionization_efficiency);
    pp.query("electrode_townsend2"       ,    m_townsend2_electrode);
    pp.query("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
    pp.query("dielectric_townsend2"       ,   m_townsend2_dielectric);
    pp.query("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
  }

  { // Noise things. Noise function is passed by pointer to species since it must be unique
    ParmParse pp("air_bolsig");
    pp.query("noise_amplitude",   m_noise_amplitude);
    pp.query("noise_octaves",     m_noise_octaves);
    pp.query("noise_persistence", m_noise_persistence);
    if(pp.contains("noise_frequency")){
      Vector<Real> freq(SpaceDim);
      pp.getarr("noise_frequency", freq, 0, SpaceDim);
      m_noise_frequency = RealVect(D_DECL(freq[0], freq[1], freq[2]));
    }
  }

  { // Transform to SI units
    m_p  *= units::s_atm2pascal;
    m_pq *= units::s_atm2pascal;
    m_N   = m_p*units::s_Na/(m_gas_temp*units::s_R);
  }

  // Compute transport data for electrons by calling BOLSIG
  this->compute_transport_coefficients();

}

air_bolsig::~air_bolsig(){

}

Vector<RealVect> air_bolsig::compute_velocities(const RealVect& a_E) const{

}


Vector<Real> air_bolsig::compute_source_terms(const Vector<Real>& a_cdr_densities, 
					      const Vector<Real>& a_rte_densities,
					      const RealVect&     a_E) const{

}


Vector<Real> air_bolsig::compute_rte_source_terms(const Vector<Real>& a_cdr_densities, const RealVect& a_E) const{

}

Vector<Real> air_bolsig::compute_diffusion_coefficients(const RealVect& a_E) const{

}


Vector<Real> air_bolsig::compute_conductor_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						  const Vector<Real>& a_cdr_densities,
						  const Vector<Real>& a_cdr_velocities,
						  const Vector<Real>& a_rte_fluxes,
						  const RealVect&     a_E,
						  const RealVect&     a_pos,
						  const RealVect&     a_normal,
						  const Real&         a_time) const{

}


Vector<Real> air_bolsig::compute_dielectric_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						   const Vector<Real>& a_cdr_densities,
						   const Vector<Real>& a_cdr_velocities,
						   const Vector<Real>& a_rte_fluxes,
						   const RealVect&     a_E,
						   const RealVect&     a_pos,
						   const RealVect&     a_normal,
						   const Real&         a_time) const{

}

Real air_bolsig::initial_sigma(const RealVect& a_pos) const {

}

void air_bolsig::compute_transport_coefficients(){
  pout() << "Computing transport data using BOLSIG- ..." << endl;

  std::stringstream ss;

  ss << procID();

  m_local_lxcat_file  = m_lxcat_file  + ss.str();
  m_local_script_file = s_script_file + ss.str();
  m_local_data_file   = s_data_file   + ss.str();

  this->build_bolsig_script();   // Build input script for BOLSIG
  this->call_bolsig();     // Call bolsigminus and output data
  this->extract_bolsig_data();   // Extract data
  this->delete_bolsig_data();    // Cleanup, delete file created by bolsigminus
  this->delete_bolsig_script();  // Cleanup, delete script

  
  pout() << "Done computing transport data" << endl;
}

void air_bolsig::build_bolsig_script(){
  // Get LxCat file
  int success;
  std::string cmd;
  cmd = "cp " + m_lxcat_path + "/" + m_lxcat_file + " ./" + m_local_lxcat_file;
  success = system(cmd.c_str());
  if(success != 0){
    MayDay::Abort("air_bolsig::build_bolsig_script - Could not get LXCat file");
  }
  
  // Open file
  ofstream file;
  file.open(m_local_script_file.c_str());

  // Comment for debugging, and turn off output
  file << "/ Comment\n";
  file << std::endl;
  file << "NOSCREEN\n";
  file << std::endl;

  // Read collisions
  file << "READCOLLISIONS\n"; 
  file << m_local_lxcat_file.c_str() << "\t / File\n";
  file << "N2 O2"                    << "\t / Species\n";
  file << 1                          << "\t / Extrapolate: 0= No 1= Yes\n";
  file << std::endl;
  
  // Set conditions
  file << "CONDITIONS\n";
  file << 0.1        << "\t / Electric field / N (Td)\n";
  file << 0.0        << "\t / Angular field frequency / N (m3/s)\n";
  file << 0.0        << "\t / Cosine of E-B field angle\n";
  file << m_gas_temp << "\t / Gas temperature (K)\n";
  file << m_gas_temp << "\t / Excitation temperature (K)\n";
  file << 0.         << "\t / Transition energy (eV)\n";
  file << 0.         << "\t / Ionization degree\n";
  file << 1e18       << "\t / Plasma density (1/m3)\n";
  file << 1          << "\t / Ion charge parameter\n";
  file << 1          << "\t / Ion/neutral mass ratio\n";
  file << 1          << "\t / e-e momentum effects: 0=No; 1=Yes*\n";
  file << 1          << "\t / Energy sharing: 1=Equal*; 2=One takes all\n";
  file << 1          << "\t / Growth: 1=Temporal*; 2=Spatial; 3=Not included; 4=Grad-n expansion\n";
  file << 0.         << "\t / Maxwellian mean energy (eV) \n";
  file << 200        << "\t / # of grid points\n";
  file << 0          << "\t / Manual grid: 0=No; 1=Linear; 2=Parabolic \n";
  file << 200.       << "\t / Manual maximum energy (eV)\n";
  file << 1e-10      << "\t / Precision\n";
  file << 1e-4       << "\t / Convergence\n";
  file << 1000       << "\t / Maximum # of iterations\n";
  file << m_frac_N2  << " "
       << m_frac_O2  << "\t / Species fractions\n";
  file << 1          << "\t / Normalize composition to unity: 0=No; 1=Yes\n";
  file << std::endl;


  // Run series1           / Variable: 1=E/N; 2=Mean energy; 3=Maxwellian energy
  const int grid_points = m_grid_points;
  const Real min        = m_min_townsend;
  const Real max        = m_max_townsend;
  
  file << "RUNSERIES\n";
  file << 1           << "\t / Variable: 1=E/N; 2=Mean energy; 3=Maxwellian energy\n";
  file << min         << " "
       << max         << "\t / Min/max\n";
  file << grid_points << "\t / Number\n";
  file << 1           << "\t / Type: 1=Linear; 2=Quadratic; 3=Exponential\n";
  file << std::endl;

  // Save results to file
  const std::string output = m_local_data_file;
  file << "SAVERESULTS\n";
  file << output.c_str() << "\t / File\n";
  file << 3              << "\t / Format: 1=Run by run; 2=Combined; 3=E/N; 4=Energy; 5=SIGLO; 6=PLASIMO\n";
  file << 1              << "\t / Conditions: 0=No; 1=Yes\n";
  file << 1              << "\t / Transport coefficients: 0=No; 1=Yes\n";
  file << 1              << "\t / Rate coefficients: 0=No; 1=Yes\n";
  file << 0              << "\t / Reverse rate coefficients: 0=No; 1=Yes\n";
  file << 0              << "\t / Energy loss coefficients: 0=No; 1=Yes\n";
  file << 0              << "\t / Distribution function: 0=No; 1=Yes \n";
  file << 0              << "\t / Skip failed runs: 0=No; 1=Yes\n";

  file << "END\n";
  file.close();
}

void air_bolsig::call_bolsig(){

  int success;
  std::string cmd;

  cmd = "./" + m_bolsig_path + "/bolsigminus " + m_local_script_file;
  success = system(cmd.c_str());
  if(success != 0){
    MayDay::Abort("air_bolsig::call_bolsig - Could not call BOLSIG");
  }

  // Delete the transient copy of the data base when we are done. 
  cmd = "rm " + m_local_lxcat_file;
  success = system(cmd.c_str());
  if(success != 0){
    MayDay::Abort("AirBolsigPlasKin::call_bolsig - Could not delete temporary file");
  }

  // Delete the bolsig log file. 
  if(procID() == 0){
    const std::string cmd = "rm bolsiglog.txt";
    const int success = system(cmd.c_str());
    if(success != 0){
      MayDay::Abort("air_bolsig::call_bolsig - could not delete log file");
    }
  }
}

void air_bolsig::extract_bolsig_data(){

}

void air_bolsig::delete_bolsig_data(){
  const std::string cmd = "rm " + m_local_data_file;
  const int success = system(cmd.c_str());
  if(success != 0){
    MayDay::Abort("air_bolsig::delete_bolsig_script - Could not delete script file");
  }
}

void air_bolsig::delete_bolsig_script(){
  const std::string cmd = "rm " + m_local_script_file;
  const int success = system(cmd.c_str());
  if(success != 0){
    MayDay::Abort("air_bolsig::delete_bolsig_script - Could not delete script file");
  }
}
