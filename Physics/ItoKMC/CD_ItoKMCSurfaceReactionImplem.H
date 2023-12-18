/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCSurfaceReactionImplem.H
  @brief  Declaration of a class that encapsulates secondary emission at EBs and domain edges/faces
  @author Robert Marskar
*/

#ifndef CD_ItoKMCSurfaceReactionImplem_H
#define CD_ItoKMCSurfaceReactionImplem_H

// Our includes
#include <CD_ItoKMCSurfaceReaction.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

inline ItoKMCSurfaceReaction::ItoKMCSurfaceReaction(const size_t            a_reactant,
                                                    const std::list<size_t> a_products,
                                                    const std::list<Real>   a_efficiencies) noexcept
{
  this->define(a_reactant, a_products, a_efficiencies);
}

inline ItoKMCSurfaceReaction::~ItoKMCSurfaceReaction() noexcept {}

inline void
ItoKMCSurfaceReaction::define(const size_t            a_reactant,
                              const std::list<size_t> a_products,
                              const std::list<Real>   a_efficiencies) noexcept
{
  m_reactant     = a_reactant;
  m_products     = a_products;
  m_efficiencies = a_efficiencies;
}

inline const size_t&
ItoKMCSurfaceReaction::getReactant() const noexcept
{
  return m_reactant;
}

inline const std::list<size_t>&
ItoKMCSurfaceReaction::getProducts() const noexcept
{
  return m_products;
}

inline const std::list<Real>&
ItoKMCSurfaceReaction::getEfficiencies() const noexcept
{
  return m_efficiencies;
}

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCSurfaceReactionImplem.H>

#endif
