#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ZELDOVICHIBC_H_
#define _ZELDOVICHIBC_H_

#include  <iostream>

#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"
#include "Particle.H"
#include "ListBox.H"
#include "ParticleData.H"
#include "Cosmology.H"

///
/**
   Initial conditions for Zel'dovich pancake.
   Standard test for cosmological code.

   Parameters:
     a_domain     - Represents the problem domain
     a_dx         - Poisson solve cell spacing
     a_offset     - Domain offset; the physical coordinate of the left edge of the box
     a_initialA   - Initial expansion factor
     a_collapseA  - Expansion factor of collapse; sets the perturb. amplitude
     a_partPerCell - number of particles per Poisson cell
     a_wavenum    - Wavenumber of initial sinusoidal plane wave
     a_cosmology - contains information on the cosmology model

 */
class ZeldovichIBC
{
public:

  /// Null constructor
  ZeldovichIBC();

  /// Full constructor
  ZeldovichIBC(const ProblemDomain& a_domain,
               const Real& a_dx,
               const RealVect& a_offset,
               const Real& a_initialA,
               const Real& a_collapseA,
               const Real& a_partPerCell,
               const Real& a_wavenum,
               const Cosmology& a_cosmology);
    
  /// Destructor
  ~ZeldovichIBC() {}

  /// Full constructor
  void define(const ProblemDomain& a_domain,
              const Real& a_dx,
              const RealVect& a_offset,
              const Real& a_initialA,
              const Real& a_collapseA,
              const Real& a_partPerCell,
              const Real& a_wavenum,
              const Cosmology& a_cosmology);
  
  /// Set up initial conditions for particles
  void initialize(ParticleData<Particle>& a_DM);

  /// factory method
  ZeldovichIBC* new_zeldovichIBC() const;

protected:
  ProblemDomain m_domain;
  Real m_dx;
  RealVect m_offset;
  Real m_initialA;
  Real m_collapseA;
  Real m_partPerCell;
  Real m_wavenum;
  Cosmology m_cosmology;
};

#endif // include guard
