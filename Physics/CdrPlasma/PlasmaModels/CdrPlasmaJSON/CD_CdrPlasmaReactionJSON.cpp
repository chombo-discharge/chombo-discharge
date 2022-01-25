/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaReactionJSON.cpp
  @brief  Implementation of CD_CdrPlasmaReactionJSON.H
  @author Robert Marskar
*/

// Our includes
#include <CD_CdrPlasmaReactionJSON.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaReactionJSON::CdrPlasmaReactionJSON() {
  m_isDefined = false;
}

CdrPlasmaReactionJSON::CdrPlasmaReactionJSON(const std::list<int> a_plasmaReactants,
					     const std::list<int> a_neutralReactants,			    
					     const std::list<int> a_plasmaProducts,
					     const std::list<int> a_photonProducts) {
  m_plasmaReactants  = a_plasmaReactants ;
  m_neutralReactants = a_neutralReactants;
  m_plasmaProducts   = a_plasmaProducts  ;
  m_photonProducts   = a_photonProducts  ;

  m_isDefined        = true;
}


CdrPlasmaReactionJSON::~CdrPlasmaReactionJSON() {

}

const std::list<int>& CdrPlasmaReactionJSON::getPlasmaReactants() const {
  return m_plasmaReactants;
}

const std::list<int>& CdrPlasmaReactionJSON::getNeutralReactants() const {
  return m_neutralReactants;
}

const std::list<int>& CdrPlasmaReactionJSON::getPlasmaProducts() const {
  return m_plasmaProducts;
}

const std::list<int>& CdrPlasmaReactionJSON::getPhotonProducts() const {
  return m_photonProducts;
}

void CdrPlasmaReactionJSON::addProducts(std::vector<Real>&       a_plasmaSources,
					std::vector<Real>&       a_photonSources,
					const std::vector<Real>& a_plasmaDensities,
					const Real               a_reactionRate) const {

  // Compute the basic rate k * n[X1] * n[X2] etc.
  Real volumetricRate = a_reactionRate;
  for (const auto& r : m_plasmaReactants){
    volumetricRate *= a_plasmaDensities[r];
  }

  // Remove mass from the left hand side. 
  for (const auto& r : m_plasmaReactants){
    a_plasmaSources[r] -= volumetricRate;
  }

  // Add mass on the right-hand side. 
  for (const auto& p : m_plasmaProducts){
    a_plasmaSources[p] += volumetricRate;
  }

  for (const auto& p : m_photonProducts){
    a_photonSources[p] += volumetricRate;
  }  
}

#include <CD_NamespaceFooter.H>
