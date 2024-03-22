/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaSurfaceReactionJSON.cpp
  @brief  Implementation of CD_CdrPlasmaSurfaceReactionJSON.H
  @author Robert Marskar
*/

// Our includes
#include <CD_CdrPlasmaSurfaceReactionJSON.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaSurfaceReactionJSON::CdrPlasmaSurfaceReactionJSON(const std::list<int> a_plasmaReactants,
                                                           const std::list<int> a_photonReactants,
                                                           const std::list<int> a_plasmaProducts)
{
  m_plasmaReactants = a_plasmaReactants;
  m_photonReactants = a_photonReactants;
  m_plasmaProducts  = a_plasmaProducts;
}

CdrPlasmaSurfaceReactionJSON::~CdrPlasmaSurfaceReactionJSON()
{}

const std::list<int>&
CdrPlasmaSurfaceReactionJSON::getPlasmaReactants() const
{
  return m_plasmaReactants;
}

const std::list<int>&
CdrPlasmaSurfaceReactionJSON::getPhotonReactants() const
{
  return m_photonReactants;
}

const std::list<int>&
CdrPlasmaSurfaceReactionJSON::getPlasmaProducts() const
{
  return m_plasmaProducts;
}

#include <CD_NamespaceFooter.H>
