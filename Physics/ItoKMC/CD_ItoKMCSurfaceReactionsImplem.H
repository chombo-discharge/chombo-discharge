/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCSurfaceReactionsImplem.H
  @brief  Implementation of CD_ItoKMCSurfaceReactions.H
  @author Robert Marskar
*/

#ifndef CD_ItoKMCSurfaceReactionsImplem_H
#define CD_ItoKMCSurfaceReactionsIMplem_H

// Our includes
#include <CD_ItoKMCSurfaceReactions.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

inline ItoKMCSurfaceReactions::ItoKMCSurfaceReactions() noexcept
{
  this->reset();
}

inline ItoKMCSurfaceReactions::~ItoKMCSurfaceReactions() noexcept
{
  this->reset();
}

inline void
ItoKMCSurfaceReactions::add(const ItoKMCSurfaceReaction& a_surfaceReaction) noexcept
{
  if (m_isDefined) {

    if (m_reactant != a_surfaceReaction.getReactant()) {
      MayDay::Warning("ItoKMCSurfaceReaction::addReaction - reactant indicies do not match!");
    }
  }
  else {
    this->reset();

    m_reactant = a_surfaceReaction.getReactant();
  }

  m_products.emplace_back(a_surfaceReaction.getProducts());
  m_efficiencies.emplace_back(a_surfaceReaction.getEfficiency());

  m_distribution = std::discrete_distribution<size_t>(std::begin(m_efficiencies), std::end(m_efficiencies));

  m_isDefined = true;
}

inline void
ItoKMCSurfaceReactions::reset() noexcept
{
  m_reactant = 0;
  m_products.clear();
  m_efficiencies.clear();

  m_isDefined = false;
}

inline const size_t&
ItoKMCSurfaceReactions::getReactant() const noexcept
{
  return m_reactant;
}

inline const std::vector<ItoKMCSurfaceReactions::Products>&
ItoKMCSurfaceReactions::getProducts() const noexcept
{
  return m_products;
}

inline const std::vector<Real>&
ItoKMCSurfaceReactions::getEfficiencies() const noexcept
{
  return m_efficiencies;
}

inline std::discrete_distribution<size_t>&
ItoKMCSurfaceReactions::getDistribution() const noexcept
{
  return m_distribution;
}

#include <CD_NamespaceFooter.H>

#endif
