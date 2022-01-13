/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaReaction.cpp
  @brief  Implementation of CD_CdrPlasmaReaction.H
  @author Robert Marskar
*/

// Our includes
#include <CD_CdrPlasmaReaction.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaReaction::CdrPlasmaReaction() {
  m_isDefined = false;
}

CdrPlasmaReaction::CdrPlasmaReaction(const std::initializer_list<int> a_particleReactants,
				     const std::initializer_list<int> a_particleProducts,
				     const std::initializer_list<int> a_photonProducts)  :  m_particleReactants(a_particleReactants),
											    m_particleProducts (a_particleProducts ),
											    m_photonProducts   (a_photonProducts   )
{
  this->computeStateChanges();
}

CdrPlasmaReaction::CdrPlasmaReaction(const std::initializer_list<int> a_particleReactants,
				     const std::initializer_list<int> a_particleProducts) : m_particleReactants(a_particleReactants),
											    m_particleProducts (a_particleProducts ),
											    m_photonProducts   (0)
{
  this->computeStateChanges();
}

CdrPlasmaReaction::~CdrPlasmaReaction() {

}

void CdrPlasmaReaction::computeStateChanges() {

  // Consumed particle species. Decrease the particle number by one for every time a particle reactant
  // appears on the left hand side. 
  for (const auto& r : m_particleReactants){
    if(m_particleJump.find(r) == m_particleJump.end()){
      m_particleJump.emplace(r,-1);
    }
    else{
      m_particleJump[r]--;
    }
  }

  // Produced particle species
  for (const auto& r : m_particleProducts){
    if(m_particleJump.find(r) == m_particleJump.end()){
      m_particleJump.emplace(r,+1);
    }
    else{
      m_particleJump[r]++;
    }
  }
}

Real& CdrPlasmaReaction::rate() const {
  return m_reactionRate;
}

void CdrPlasmaReaction::fire(std::vector<Real>&       a_particleSources,
			     std::vector<Real>&       a_photonSources,
			     const std::vector<Real>& a_particleDensities) const {
  this->fire(a_particleSources, a_photonSources, a_particleDensities, m_reactionRate);
}

void CdrPlasmaReaction::fire(std::vector<Real>&       a_particleSources,
			     std::vector<Real>&       a_photonSources,
			     const std::vector<Real>& a_particleDensities,
			     const Real               a_reactionRate) const {

  // Compute the basic rate k * n[X1] * n[X2] etc.
  Real volumetricRate = a_reactionRate;
  for (const auto& r : m_particleReactants){
    volumetricRate *= a_particleDensities[r];
  }

  // Consumed mass on left-hand side. 
  for (const auto& r : m_particleReactants){
    a_particleSources[r] -= volumetricRate;
  }

  // Produced mass on right-hand side.
  for (const auto& p : m_particleProducts){
    a_particleSources[p] += volumetricRate;
  }

  for (const auto& p : m_photonProducts){
    a_photonSources[p] += volumetricRate;
  }  
}

#include <CD_NamespaceFooter.H>
