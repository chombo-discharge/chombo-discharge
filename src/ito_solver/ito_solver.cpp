/*!
  @file   ito_solver.H
  @brief  Declaration of an abstract class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
*/

#include "ito_solver.H"

#include <ParmParse.H>

ito_solver::ito_solver(){
  m_name       = "ito_solver";
  m_class_name = "ito_solver";
}

ito_solver::~ito_solver(){

}

std::string ito_solver::get_name(){
  return m_name;
}

void ito_solvre::parse_options(){
  CH_TIME("ito_solver::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }

  this->parse_rng();
  this->parse_plot_vars();
  this->parse_deposition();
  this->parse_bisect_step();
  this->parse_domain_bc();
  this->parse_pvr_buffer();
}

void ito_solver::parse_rng(){
  CH_TIME("ito_solver::parse_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_rng" << endl;
  }
}

void ito_solver::parse_plot_vars(){
  CH_TIME("ito_solver::parse_plot_vars");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_plot_vars" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_class_name.c_str());
  pp.get("seed", m_seed);
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  
  m_rng = std::mt19937_64(m_seed);

  m_udist01 = std::uniform_real_distribution<Real>( 0.0, 1.0);
  m_udist11 = std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_gauss01 = std::normal_distribution<Real>(0.0, 1.0);
}

void ito_solver::parse_deposition(){
  CH_TIME("ito_solver::parse_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_rng" << endl;
  }
}

void ito_solver::parse_bisect_step(){
  CH_TIME("ito_solver::parse_bisect_step");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_bisect_step" << endl;
  }
}

void ito_solver::parse_domain_bc(){
  CH_TIME("ito_solver::parse_domain_bc");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_domain_bc" << endl;
  }
}

void ito_solver::parse_pvr_buffer(){
  CH_TIME("ito_solver::parse_pvr_buffer");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_pvr_buffer" << endl;
  }
}

Vector<std::string> ito_solver::get_plotvar_names() const {
  CH_TIME("ito_solver::get_plotvar_names");
  if(m_verbosity > 5){
    pout() << m_name + "::get_plotvar_names" << endl;
  }

  Vector<std::string> m_names(0);
  if(m_plot_phi) names.push_back(m_name + " phi");

  return names;
}

void ito_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("ito_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;
}

void ito_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("ito_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }

  m_amr = a_amr;
}

void set_phase(phase::which_phase a_phase){
  CH_TIME("ito_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void ito_solver::set_verbosity(const int a_verbosity){
  CH_TIME("ito_solver::set_verbosity");
  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
}
  
void set_plot_variables(){

}

void ito_solver::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("ito_solver::set_time");
  if(m_verbosity > 5){
    pout() << m_name + "::set_time" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}
