/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir7StephensSpecies.H
  @brief  Implementation of CD_CdrPlasmaAir7StephensSpecies.cpp 
  @author Robert Marskar
*/

// Std includes
#include <chrono>

// Chombo includes
#include <ParmParse.H>

#include <CD_CdrPlasmaAir7Stephens.H>
#include <CD_CdrPlasmaAir7StephensSpecies.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaAir7Stephens::Electron::Electron()
{
  m_name         = "Electron";
  m_chargeNumber = -1;
  ParmParse    pp("CdrPlasmaAir7Stephens");
  std::string  str;
  Vector<Real> vec(SpaceDim);

  pp.get("mobile_electrons", str);
  m_isMobile = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str);
  m_isDiffusive = (str == "true") ? true : false;

  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius", m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim);
  m_seed_pos = RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

CdrPlasmaAir7Stephens::N2plus::N2plus()
{
  m_name         = "N2plus";
  m_chargeNumber = 1;
  ParmParse    pp("CdrPlasmaAir7Stephens");
  std::string  str;
  Vector<Real> vec(SpaceDim);

  pp.get("mobile_ions", str);
  m_isMobile = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);
  m_isDiffusive = (str == "true") ? true : false;

  pp.get("frac_N2", m_frac);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius", m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim);
  m_seed_pos = RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

CdrPlasmaAir7Stephens::O2plus::O2plus()
{
  m_name         = "O2plus";
  m_chargeNumber = 1;
  ParmParse    pp("CdrPlasmaAir7Stephens");
  std::string  str;
  Vector<Real> vec(SpaceDim);

  pp.get("mobile_ions", str);
  m_isMobile = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);
  m_isDiffusive = (str == "true") ? true : false;

  pp.get("frac_O2", m_frac);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius", m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim);
  m_seed_pos = RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

CdrPlasmaAir7Stephens::N4plus::N4plus()
{
  m_name         = "N4plus";
  m_chargeNumber = 1;
  ParmParse   pp("CdrPlasmaAir7Stephens");
  std::string str;

  pp.get("mobile_ions", str);
  m_isMobile = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);
  m_isDiffusive = (str == "true") ? true : false;
}

CdrPlasmaAir7Stephens::O4plus::O4plus()
{
  m_name         = "O4plus";
  m_chargeNumber = 1;
  ParmParse   pp("CdrPlasmaAir7Stephens");
  std::string str;

  pp.get("mobile_ions", str);
  m_isMobile = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);
  m_isDiffusive = (str == "true") ? true : false;
}

CdrPlasmaAir7Stephens::O2plusN2::O2plusN2()
{
  m_name         = "O2plusN2";
  m_chargeNumber = 1;
  ParmParse   pp("CdrPlasmaAir7Stephens");
  std::string str;

  pp.get("mobile_ions", str);
  m_isMobile = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);
  m_isDiffusive = (str == "true") ? true : false;
}

CdrPlasmaAir7Stephens::O2minus::O2minus()
{
  m_name         = "O2minus";
  m_chargeNumber = -1;
  ParmParse   pp("CdrPlasmaAir7Stephens");
  std::string str;

  pp.get("mobile_ions", str);
  m_isMobile = (str == "true") ? true : false;
  pp.get("diffusive_ions", str);
  m_isDiffusive = (str == "true") ? true : false;
}

Real
CdrPlasmaAir7Stephens::Electron::initialData(const RealVect a_pos, const Real a_time) const
{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_uniform_density + m_seed_density * exp(-factor * factor / (m_seed_rad * m_seed_rad));
}

Real
CdrPlasmaAir7Stephens::N2plus::initialData(const RealVect a_pos, const Real a_time) const
{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_frac * (m_uniform_density + m_seed_density * exp(-factor * factor / (m_seed_rad * m_seed_rad)));
}

Real
CdrPlasmaAir7Stephens::O2plus::initialData(const RealVect a_pos, const Real a_time) const
{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_frac * (m_uniform_density + m_seed_density * exp(-factor * factor / (m_seed_rad * m_seed_rad)));
}

CdrPlasmaAir7Stephens::phot_c4v0_X1v0::phot_c4v0_X1v0()
{
  m_name = "phot_c4v0_X1v0";

  Real      p;
  ParmParse pp("CdrPlasmaAir7Stephens");

  pp.get("pressure", p);
  pp.get("c4v0_X1v0_beer", m_kappa);
  m_kappa = 1. / (m_kappa * p);
}

CdrPlasmaAir7Stephens::phot_c4v0_X1v1::phot_c4v0_X1v1()
{
  m_name = "phot_c4v0_X1v1";

  Real      p;
  ParmParse pp("CdrPlasmaAir7Stephens");

  pp.get("pressure", p);
  pp.get("c4v0_X1v1_beer", m_kappa);
  m_kappa = 1. / (m_kappa * p);
}

CdrPlasmaAir7Stephens::phot_c4v1_X1v0::phot_c4v1_X1v0()
{
  m_name = "phot_c4v1_X1v0";

  Real      p;
  ParmParse pp("CdrPlasmaAir7Stephens");

  pp.get("pressure", p);
  pp.get("c4v1_X1v0_beer", m_kappa);

  m_kappa = 1. / (m_kappa * p);
}

CdrPlasmaAir7Stephens::phot_c4v1_X1v1::phot_c4v1_X1v1()
{
  m_name = "phot_c4v1_X1v1";

  Real      p;
  ParmParse pp("CdrPlasmaAir7Stephens");

  pp.get("pressure", p);
  pp.get("c4v1_X1v1_beer", m_kappa);

  m_kappa = 1. / (m_kappa * p);
}

CdrPlasmaAir7Stephens::phot_c4v1_X1v2::phot_c4v1_X1v2()
{
  m_name = "phot_c4v1_X1v2";

  Real      p;
  ParmParse pp("CdrPlasmaAir7Stephens");

  pp.get("pressure", p);
  pp.get("c4v1_X1v2_beer", m_kappa);

  m_kappa = 1. / (m_kappa * p);
}

CdrPlasmaAir7Stephens::phot_c4v1_X1v3::phot_c4v1_X1v3()
{
  m_name = "phot_c4v1_X1v3";

  Real      p;
  ParmParse pp("CdrPlasmaAir7Stephens");

  pp.get("pressure", p);
  pp.get("c4v1_X1v3_beer", m_kappa);

  m_kappa = 1. / (m_kappa * p);
}

CdrPlasmaAir7Stephens::phot_b1v1_X1v0::phot_b1v1_X1v0()
{
  m_name = "phot_b1v1_X1v0";

  Real      p;
  ParmParse pp("CdrPlasmaAir7Stephens");

  pp.get("pressure", p);
  pp.get("b1v1_X1v0_beer", m_kappa);

  m_kappa = 1. / (m_kappa * p);
}

CdrPlasmaAir7Stephens::phot_b1v1_X1v1::phot_b1v1_X1v1()
{
  m_name = "phot_b1v1_X1v1";

  Real      p;
  ParmParse pp("CdrPlasmaAir7Stephens");

  pp.get("pressure", p);
  pp.get("b1v1_X1v1_beer", m_kappa);

  m_kappa = 1. / (m_kappa * p);
}

#include <CD_NamespaceFooter.H>
