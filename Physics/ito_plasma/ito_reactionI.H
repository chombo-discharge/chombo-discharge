/*!
  @file   ito_reactionI.H
  @brief  Implementation of ito_reaction.H
  @author Robert Marskar
  @date   Aug. 2020
*/

#ifndef _ITO_REACTIONI_H_
#define _ITO_REACTIONI_H_

#include "ito_reaction.H"

#define ITO_REACTION_DEBUG 0

#include "CD_NamespaceHeader.H"
using namespace Physics::ito_plasma;

inline
void ito_reaction::jump_state(Vector<long long>& a_particles,
			      const long long    a_num_reactions) const {
#if ITO_REACTION_DEBUG
  if(a_num_reactions < 0) {
    MayDay::Abort("ito_reaction::jump_state - can't have negative number of reactions!");
  }
  for (int i = 0; i < a_particles.size(); i++){
    if(a_particles[i] < 0){
      MayDay::Abort("ito_reaction::jump_state can't have negative particles on the way in");
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
      MayDay::Abort("ito_reaction::jump_state can't have negative particles on the way out");
    }
  }
#endif
}

inline
void ito_reaction::jump_state(Vector<long long>& a_particles,
			      Vector<long long>& a_Photons,
			      const long long    a_num_reactions) const {

  this->jump_state(a_particles, a_num_reactions);

  for (const auto& r : m_Photon_products){
    a_Photons[r] += a_num_reactions;;
  }
}

inline
void ito_reaction::jump_energy(Vector<Real>& a_energies, const Vector<Real>& a_mean_energies, const long long a_num_reactions) const {
  for (const auto& r : m_energy_jumps){
    if(this->getPhi_change(r.first) < 0) { // Consuming reaction, this eats a_mean_energies.  
      a_energies[r.first] -= a_mean_energies[r.first]*a_num_reactions;
    }
    else{
      a_energies[r.first] += r.second*a_num_reactions;
    }
  }
}

inline
void ito_reaction::jump_energy(Vector<Real>&       a_energies,
			       const Vector<Real>& a_mean_energies,
			       const Vector<Real>& a_sources,
			       const long long     a_num_reactions,
			       const Real          a_dt) const {

  for (const auto& r : m_energy_jumps){
    if(this->getPhi_change(r.first) < 0) { // Consuming reaction, this eats a_mean_energies.  
      a_energies[r.first] -= a_mean_energies[r.first]*a_num_reactions;
    }
    else{
      a_energies[r.first] += r.second*a_num_reactions;
    }

    a_energies[r.first] += a_sources[r.first]*a_dt;
  }
}

inline
Real& ito_reaction::rate() const {
  return m_rate;
}

inline
Real ito_reaction::propensity(const Vector<long long>& a_particles) const {
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
int ito_reaction::get_gi(const Vector<long long>& a_particles) const{
  return 1; // For simplicity
}

inline
const std::list<int>& ito_reaction::get_reactants() const{
  return m_reactants;
}

inline
const std::list<int>& ito_reaction::get_particle_products() const{
  return m_particle_products;
}

inline
const std::list<int>& ito_reaction::get_Photon_products() const{
  return m_Photon_products;
}

inline
void ito_reaction::compute_state_change() {

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
void ito_reaction::compute_energy_change() {
  for (const auto& r : m_energy_jumps){
    if(m_energyChange.find(r.first) == m_energyChange.end()){
      m_energyChange.emplace(r.first, r.second);
    }
  }
}

inline
const std::map<int, int>& ito_reaction::getPhi_change() const {
  return m_phiChange;
}

inline
const int ito_reaction::getPhi_change(const int a_idx) const {
  int ret = 0;
  
  if(m_phiChange.find(a_idx) != m_phiChange.end()){
    ret = m_phiChange.at(a_idx);
  }

  return ret;
}

inline
const std::map<int, Real>& ito_reaction::get_energy_change() const {
  return m_energyChange;
}

inline
const Real ito_reaction::get_energy_change(const int a_idx) const {
  Real ret = 0;
  
  if(m_energyChange.find(a_idx) != m_energyChange.end()){
    ret = m_energyChange.at(a_idx);
  }

  return ret;
}
#include "CD_NamespaceFooter.H"

#endif
