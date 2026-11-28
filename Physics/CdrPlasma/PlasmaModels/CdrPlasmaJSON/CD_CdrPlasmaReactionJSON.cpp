/*
 * SPDX-FileCopyrightText: 2022-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CdrPlasmaReactionJSON.cpp
 * @brief  Implementation of CD_CdrPlasmaReactionJSON.H
 * @author Robert Marskar
 */

// Our includes
#include <CD_CdrPlasmaReactionJSON.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaReactionJSON::CdrPlasmaReactionJSON(const std::list<int>& a_plasmaReactants,
                                             const std::list<int>& a_neutralReactants,
                                             const std::list<int>& a_plasmaProducts,
                                             const std::list<int>& a_photonProducts)
  : m_plasmaReactants(a_plasmaReactants.begin(), a_plasmaReactants.end()),
    m_neutralReactants(a_neutralReactants.begin(), a_neutralReactants.end()),
    m_plasmaProducts(a_plasmaProducts.begin(), a_plasmaProducts.end()),
    m_photonProducts(a_photonProducts.begin(), a_photonProducts.end())
{}

CdrPlasmaReactionJSON::~CdrPlasmaReactionJSON() = default;

const std::vector<int>&
CdrPlasmaReactionJSON::getPlasmaReactants() const
{
  return m_plasmaReactants;
}

const std::vector<int>&
CdrPlasmaReactionJSON::getNeutralReactants() const
{
  return m_neutralReactants;
}

const std::vector<int>&
CdrPlasmaReactionJSON::getPlasmaProducts() const
{
  return m_plasmaProducts;
}

const std::vector<int>&
CdrPlasmaReactionJSON::getPhotonProducts() const
{
  return m_photonProducts;
}

#include <CD_NamespaceFooter.H>
