/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaReactionImplem.H
  @brief  Implementation of CD_ItoPlasmaReaction.H
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaReactionImplem_H
#define CD_ItoPlasmaReactionImplem_H

// Our includes
#include <CD_ItoPlasmaReaction.H>
#include <CD_NamespaceHeader.H>

#define ITO_REACTION_DEBUG 0

using namespace Physics::ItoPlasma;

inline
void ItoPlasmaReaction::jumpState(Vector<long long>& a_particles,
				  const long long    a_num_reactions) const {
#if ITO_REACTION_DEBUG
  if(a_num_reactions < 0) {
    MayDay::Abort("ItoPlasmaReaction::jumpState - can't have negative number of reactions!");
  }
  for (int i = 0; i < a_particles.size(); i++){
    if(a_particles[i] < 0){
      MayDay::Abort("ItoPlasmaReaction::jumpState can't have negative particles on the way in");
    }
  }
#endif
  
  for (const auto& r : m_reactants){
    a_particles[r] -= a_num_reactions;
  }

  for (const auto& r : m_particle_products){
    a_particles[r] += a_num_reactions;
  }

#if ITO_REACTION_DEBUG
  for (int i = 0; i < a_particles.size(); i++){
    if(a_particles[i] < 0){
      MayDay::Abort("ItoPlasmaReaction::jumpState can't have negative particles on the way out");
    }
  }
#endif
}

inline
void ItoPlasmaReaction::jumpState(Vector<long long>& a_particles,
				  Vector<long long>& a_Photons,
				  const long long    a_num_reactions) const {

  this->jumpState(a_particles, a_num_reactions);

  for (const auto& r : m_photon_products){
    a_Photons[r] += a_num_reactions;;
  }
}

inline
void ItoPlasmaReaction::jumpEnergy(Vector<Real>& a_energies, const Vector<Real>& a_mean_energies, const long long a_num_reactions) const {
  for (const auto& r : m_energy_jumps){
    if(this->getStateChange(r.first) < 0) { // Consuming reaction, this eats a_mean_energies.  
      a_energies[r.first] -= a_mean_energies[r.first]*a_num_reactions;
    }
    else{
      a_energies[r.first] += r.second*a_num_reactions;
    }
  }
}

inline
void ItoPlasmaReaction::jumpEnergy(Vector<Real>&       a_energies,
				   const Vector<Real>& a_mean_energies,
				   const Vector<Real>& a_sources,
				   const long long     a_num_reactions,
				   const Real          a_dt) const {

  for (const auto& r : m_energy_jumps){
    if(this->getStateChange(r.first) < 0) { // Consuming reaction, this eats a_mean_energies.  
      a_energies[r.first] -= a_mean_energies[r.first]*a_num_reactions;
    }
    else{
      a_energies[r.first] += r.second*a_num_reactions;
    }

    a_energies[r.first] += a_sources[r.first]*a_dt;
  }
}

inline
Real& ItoPlasmaReaction::rate() const {
  return m_rate;
}

inline
Real ItoPlasmaReaction::propensity(const Vector<long long>& a_particles) const {
  Real a = m_rate;

  // Make a local copy of the particles
  Vector<long long> particles = a_particles;

  // Iterate over reactants, if a reactant appears twice the propensity function is X*(X-1)
  for (const auto& r : m_reactants){
    a *= particles[r];
    particles[r]--;
  }

  return Max(0.0, a);
}

inline
int ItoPlasmaReaction::getGi(const Vector<long long>& a_particles) const{
  return 1; // For simplicity
}

inline
const std::list<int>& ItoPlasmaReaction::getReactants() const{
  return m_reactants;
}

inline
const std::list<int>& ItoPlasmaReaction::getParticleProdcuts() const{
  return m_particle_products;
}

inline
const std::list<int>& ItoPlasmaReaction::getPhotonProducts() const{
  return m_photon_products;
}

inline
void ItoPlasmaReaction::computeStateChange() {

  // Consumed species
  for (const auto& r : m_reactants){
    if(m_phiChange.find(r) == m_phiChange.end()){
      m_phiChange.emplace(r, -1);
    }
    else{
      m_phiChange[r]--;
    }
  }

  // Produced species
  for (const auto& r : m_particle_products){
    if(m_phiChange.find(r) == m_phiChange.end()){
      m_phiChange.emplace(r, +1);
    }
    else{
      m_phiChange[r]++;
    }
  }
}

inline
void ItoPlasmaReaction::computeEnergyChange() {
  for (const auto& r : m_energy_jumps){
    if(m_energyChange.find(r.first) == m_energyChange.end()){
      m_energyChange.emplace(r.first, r.second);
    }
  }
}

inline
const std::map<int, int>& ItoPlasmaReaction::getStateChange() const {
  return m_phiChange;
}

inline
const int ItoPlasmaReaction::getStateChange(const int a_idx) const {
  int ret = 0;
  
  if(m_phiChange.find(a_idx) != m_phiChange.end()){
    ret = m_phiChange.at(a_idx);
  }

  return ret;
}

inline
const std::map<int, Real>& ItoPlasmaReaction::getEnergychange() const {
  return m_energyChange;
}

inline
const Real ItoPlasmaReaction::getEnergychange(const int a_idx) const {
  Real ret = 0;
  
  if(m_energyChange.find(a_idx) != m_energyChange.end()){
    ret = m_energyChange.at(a_idx);
  }

  return ret;
}

#include <CD_NamespaceFooter.H>

#endif
