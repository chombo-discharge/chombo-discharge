/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RandomImplem.H
  @brief  Implementation of CD_Random.H
  @author Robert Marskar
*/

#ifndef CD_RandomImplem_H
#define CD_RandomImplem_H

// Std includes
#include <chrono>

// Chombo includes
#include <SPMD.H>

// Our includes
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

std::mt19937_64                      Random::s_rng       = std::mt19937_64(0);
std::uniform_real_distribution<Real> Random::s_uniform01 = std::uniform_real_distribution<Real>( 0.0,  1.0);
std::uniform_real_distribution<Real> Random::s_uniform11 = std::uniform_real_distribution<Real>(-1.0,  1.0);
std::normal_distribution      <Real> Random::s_normal01  = std::normal_distribution      <Real>( 0.0,  1.0);

void Random::setSeed(const int a_seed) {
#ifdef CH_MPI
  const int seed = a_seed + procID();
#else
  const int seed = a_seed;
#endif

  s_rng = std::mt19937_64(seed);
}

void Random::init() {
  int seed = std::chrono::system_clock::now().time_since_epoch().count();

#ifdef CH_MPI
  seed += procID();
#endif
  
  s_rng = std::mt19937_64(seed);
}

Real Random::getUniformReal01() {
  return s_uniform01(s_rng);
}

Real Random::getUniformReal11() {
  return s_uniform11(s_rng);
}

Real Random::getNormal01() {
  return s_normal01(s_rng);
}

RealVect Random::getDirection() {
  constexpr Real safety = 1.E-12;

#if CH_SPACEDIM==2
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < safety){
    x1 = Random::getUniformReal11();
    x2 = Random::getUniformReal11();
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
#elif CH_SPACEDIM==3
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < safety){
    x1 = Random::getUniformReal11();
    x2 = Random::getUniformReal11();    
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
#endif
}

template<typename T>
Real Random::get(T& a_distribution){
  return a_distribution(s_rng);
}

#include <CD_NamespaceFooter.H>  

#endif
