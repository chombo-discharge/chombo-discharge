#include "Cosmology.H"
#include "CH_Timer.H"
#include <cmath>
#include "NamespaceHeader.H"

Cosmology::Cosmology(const Real& a_H0)
{
  CH_TIME("Cosmology::Cosmology");
  m_H0 = a_H0;
  m_Omega_m = 1.0;
}

void Cosmology::define(const Real& a_H0)
{
  CH_TIME("Cosmology::define");
  m_H0 = a_H0;
  m_Omega_m = 1.0;
}

Real Cosmology::a(const Real& a_time)
{
  CH_TIME("Cosmology::a");
  return pow(1.5*m_H0*a_time, 2.0/3.0);
}

Real Cosmology::adot(const Real& a_time)
{
  CH_TIME("Cosmology::adot");
  return (2.0/3.0)*pow(1.5*m_H0, 2.0/3.0) * pow(a_time, -1.0/3.0);
}

Real Cosmology::t(const Real& a_expansion_factor)
{
  CH_TIME("Cosmology::t");
  return 2.0*pow(a_expansion_factor, 1.5) / 3.0 / m_H0;
}

#include "NamespaceFooter.H"
