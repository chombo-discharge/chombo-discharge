/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBHelmholtzLarsenDomainBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzLarsenDomainBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzLarsenDomainBCFactory.H>
#include <CD_EBHelmholtzLarsenDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzLarsenDomainBCFactory::EBHelmholtzLarsenDomainBCFactory(const RefCountedPtr<RtSpecies>& a_species,
                                                                   const Real                      a_r1,
                                                                   const Real                      a_r2,
                                                                   const SourceFunction            a_source)
  : m_species(a_species), m_r1(a_r1), m_r2(a_r2), m_source(a_source)
{}

EBHelmholtzLarsenDomainBCFactory::EBHelmholtzLarsenDomainBCFactory(const RefCountedPtr<RtSpecies>& a_species,
                                                                   const Real                      a_r1,
                                                                   const Real                      a_r2)
  : m_species(a_species), m_r1(a_r1), m_r2(a_r2)
{

  m_source = [](const RealVect& a_position) {
    return 0.0;
  };
}

EBHelmholtzLarsenDomainBCFactory::~EBHelmholtzLarsenDomainBCFactory() = default;

RefCountedPtr<EBHelmholtzDomainBC>
EBHelmholtzLarsenDomainBCFactory::create() const
{
  return RefCountedPtr<EBHelmholtzLarsenDomainBC>(new EBHelmholtzLarsenDomainBC(m_species, m_r1, m_r2, m_source));
}

#include <CD_NamespaceFooter.H>
