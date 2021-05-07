/*!
  @file air_zheleznyak_species.H
*/

#include "air7_zheleznyak.H"
#include "air7_zheleznyak_species.H"
#include "units.H"

#include <chrono>

namespace ChomboDischarge {
  using namespace physics::cdr_plasma;

  air7_zheleznyak::electron::electron(){
    m_name = "electron";
    m_charge = -1;
    ParmParse pp("air7_zheleznyak");
    std::string str;
    Vector<Real> vec(SpaceDim);
  
    pp.get("mobile_electrons", str);    m_mobile    = (str == "true") ? true : false;
    pp.get("diffusive_electrons", str); m_diffusive = (str == "true") ? true : false;

    pp.get("uniform_density", m_uniform_density);
    pp.get("seed_density", m_seed_density);
    pp.get("seed_radius",   m_seed_rad);
    pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
  }

  air7_zheleznyak::N2plus::N2plus(){
    m_name = "N2plus";
    m_charge = 1;
    ParmParse pp("air7_zheleznyak");
    std::string str;
    Vector<Real> vec(SpaceDim);
  
    pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
    pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;

    pp.get("frac_N2", m_frac);
    pp.get("uniform_density", m_uniform_density);
    pp.get("seed_density", m_seed_density);
    pp.get("seed_radius",   m_seed_rad);
    pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
  }

  air7_zheleznyak::O2plus::O2plus(){
    m_name = "O2plus";
    m_charge = 1;
    ParmParse pp("air7_zheleznyak");
    std::string str;
    Vector<Real> vec(SpaceDim);
  
    pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
    pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;

    pp.get("frac_O2", m_frac);
    pp.get("uniform_density", m_uniform_density);
    pp.get("seed_density", m_seed_density);
    pp.get("seed_radius",   m_seed_rad);
    pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
  }

  air7_zheleznyak::N4plus::N4plus(){
    m_name = "N4plus";
    m_charge = 1;
    ParmParse pp("air7_zheleznyak");
    std::string str;
  
    pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
    pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
  }

  air7_zheleznyak::O4plus::O4plus(){
    m_name = "O4plus";
    m_charge = 1;
    ParmParse pp("air7_zheleznyak");
    std::string str;
  
    pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
    pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
  }

  air7_zheleznyak::O2plusN2::O2plusN2(){
    m_name = "O2plusN2";
    m_charge = 1;
    ParmParse pp("air7_zheleznyak");
    std::string str;
  
    pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
    pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
  }

  air7_zheleznyak::O2minus::O2minus(){
    m_name = "O2minus";
    m_charge = -1;
    ParmParse pp("air7_zheleznyak");
    std::string str;
  
    pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
    pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
  }

  Real air7_zheleznyak::electron::initial_data(const RealVect a_pos, const Real a_time) const{
    const Real factor = (a_pos - m_seed_pos).vectorLength();
    return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
  }

  Real air7_zheleznyak::N2plus::initial_data(const RealVect a_pos, const Real a_time) const{
    const Real factor = (a_pos - m_seed_pos).vectorLength();
    return m_frac*(m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad)));
  }

  Real air7_zheleznyak::O2plus::initial_data(const RealVect a_pos, const Real a_time) const{
    const Real factor = (a_pos - m_seed_pos).vectorLength();
    return m_frac*(m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad)));
  }

  air7_zheleznyak::uv_photon::uv_photon(){
    m_name   = "uv_photon";

    Real pressure, O2_frac;
  
    ParmParse pp("air7_zheleznyak");
    pp.get("photoi_f1",      m_f1);
    pp.get("photoi_f2",      m_f2);
    pp.get("photoi_K1",      m_K1);
    pp.get("photoi_K2",      m_K2);

    pp.get("frac_O2",  O2_frac);
    pp.get("pressure", pressure);
    pp.get("photoi_seed",     m_seed);

    // Convert units
    m_pO2 = pressure*O2_frac*units::s_atm2pascal;
    m_K1  = m_K1*m_pO2;
    m_K2  = m_K2*m_pO2;

    // Seed the RNG
    if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
    m_rng = new std::mt19937_64(m_seed);
    m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
  }

  air7_zheleznyak::uv_photon::~uv_photon(){
  
  }

  Real air7_zheleznyak::uv_photon::get_kappa(const RealVect a_pos) const {
    return get_random_kappa();
  }

  Real air7_zheleznyak::uv_photon::get_random_kappa() const {
    const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
    return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
  }
}
