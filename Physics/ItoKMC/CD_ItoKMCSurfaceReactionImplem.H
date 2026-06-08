/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCSurfaceReactionImplem.H
  @brief  Implementation of CD_ItoKMCSurfaceReaction.H
  @author Robert Marskar
*/

#ifndef CD_ITOKMCSURFACEREACTIONIMPLEM_H
#define CD_ITOKMCSURFACEREACTIONIMPLEM_H

// Our includes
#include <CD_ItoKMCSurfaceReaction.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

inline ItoKMCSurfaceReaction::ItoKMCSurfaceReaction(const size_t            a_reactant,
                                                    const std::list<size_t> a_products,
                                                    const Real              a_efficiency) noexcept
{
  this->define(a_reactant, a_products, a_efficiency);
}

inline ItoKMCSurfaceReaction::~ItoKMCSurfaceReaction() noexcept
{}

inline void
ItoKMCSurfaceReaction::define(const size_t            a_reactant,
                              const std::list<size_t> a_products,
                              const Real              a_efficiency) noexcept
{
  m_reactant   = a_reactant;
  m_products   = a_products;
  m_efficiency = a_efficiency;
}

inline const size_t&
ItoKMCSurfaceReaction::getReactant() const noexcept
{
  return m_reactant;
}

inline size_t&
ItoKMCSurfaceReaction::getReactant() noexcept
{
  return m_reactant;
}

inline const std::list<size_t>&
ItoKMCSurfaceReaction::getProducts() const noexcept
{
  return m_products;
}

inline std::list<size_t>&
ItoKMCSurfaceReaction::getProducts() noexcept
{
  return m_products;
}

inline const Real&
ItoKMCSurfaceReaction::getEfficiency() const noexcept
{
  return m_efficiency;
}

inline Real&
ItoKMCSurfaceReaction::getEfficiency() noexcept
{
  return m_efficiency;
}

#include <CD_NamespaceFooter.H>

#endif
