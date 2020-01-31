#ifndef _COSMOLOGY_H_
#define _COSMOLOGY_H_

#include "REAL.H"
#include "NamespaceHeader.H"

class Cosmology
{
public:

  Cosmology(){}

  /// This class stores all the cosmology parameters for
  /// the Zeldovich Pancake test problem. Basically, it 
  /// contains routines for converting between proper time
  /// and the expansion factor of the universe. Note that 
  /// currently it is hard-coded for a critical density, 
  /// matter-dominated universe (Omega_m = 1, others are 0)
  Cosmology(const Real& a_H0);

  /// same as constructor
  void define(const Real& a_H0);

  /// returns the expansion factor associated with a given
  /// proper time
  Real a(const Real& a_time);

  /// returns the time derivative of the expansion factor
  /// at a given proper time
  Real adot(const Real& a_time);

  /// returns the proper time associated with a given 
  /// expansion factor
  Real t(const Real& a_expansion_factor);

  /// the present-day Hubble parameter and the matter density
  /// of the universe in critical units
  Real m_H0, m_Omega_m;
};

#include "NamespaceFooter.H"
#endif // include guard
