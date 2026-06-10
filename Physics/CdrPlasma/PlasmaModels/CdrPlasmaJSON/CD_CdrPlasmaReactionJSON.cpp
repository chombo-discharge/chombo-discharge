/*
 * SPDX-FileCopyrightText: 2022-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_CdrPlasmaReactionJSON.cpp
   @brief  Implementation of CD_CdrPlasmaReactionJSON.H
   @author Robert Marskar
*/

// Our includes
#include <CD_CdrPlasmaReactionJSON.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaReactionJSON::CdrPlasmaReactionJSON(const std::list<int>& a_plasmaReactants,
                                             const std::list<int>& a_neutralReactants,
                                             const std::list<int>& a_plasmaProducts,
                                             const std::list<int>& a_photonProducts)
  : m_plasmaReactants(a_plasmaReactants),
    m_neutralReactants(a_neutralReactants),
    m_photonProducts(a_photonProducts),
    m_plasmaProducts(a_plasmaProducts)
{}

CdrPlasmaReactionJSON::~CdrPlasmaReactionJSON() = default;

const std::list<int>&
CdrPlasmaReactionJSON::getPlasmaReactants() const
{
  return m_plasmaReactants;
}

const std::list<int>&
CdrPlasmaReactionJSON::getNeutralReactants() const
{
  return m_neutralReactants;
}

const std::list<int>&
CdrPlasmaReactionJSON::getPlasmaProducts() const
{
  return m_plasmaProducts;
}

const std::list<int>&
CdrPlasmaReactionJSON::getPhotonProducts() const
{
  return m_photonProducts;
}

#include <CD_NamespaceFooter.H>
