/*!
  @file   ito_reactionI.H
  @brief  Implementation of ito_reaction.H
  @author Robert Marskar
  @date   Aug. 2020
*/

#ifndef _ITO_REACTIONI_H_
#define _ITO_REACTIONI_H_

#include "ito_reaction.H"

using namespace physics::ito_plasma;

inline void ito_reaction::jump_state(Vector<int>& a_particles,
				     const int    a_num_reactions) const {
  for (const auto& r : m_reactants){
    a_particles[r] -= a_num_reactions;
  }

  for (const auto& r : m_particle_products){
    a_particles[r] += a_num_reactions;;
  }
}

inline void ito_reaction::jump_state(Vector<int>& a_particles,
				     Vector<int>& a_photons,
				     const int    a_num_reactions) const {

  jump_state(a_particles, a_num_reactions);

  for (const auto& r : m_photon_products){
    a_photons[r] += a_num_reactions;;
  }
}

inline Real& ito_reaction::rate() const {
  return m_rate;
}

inline Real ito_reaction::propensity(const Vector<int>& a_particles) const {
  Real a = m_rate;

  // Make a local copy of the particles
  Vector<int> particles = a_particles;

  // Iterate over reactants, if a reactant appears twice the propensity function is X*(X-1)
  for (int i = 0; i < m_reactants.size(); i++){
    const int reactant = m_reactants[i];
	  
    a *= particles[reactant];
    particles[reactant]--;
  }

  return Max(0.0, a);
}

#endif
