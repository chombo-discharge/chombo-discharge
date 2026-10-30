/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMassRedistribution.H
 * @brief  Declaration of the redistribution of mass deposited in cut cells.
 * @author Robert Marskar
 */

#ifndef CD_PARTICLEMASSREDISTRIBUTION_H
#define CD_PARTICLEMASSREDISTRIBUTION_H

// Chombo includes
#include <RefCountedPtr.H>
#include <Vector.H>

// Our includes
#include <CD_EBFluxRedistribution.H>
#include <CD_EBAMRData.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Redistribution of mass deposited in cut cells to the neighbouring cells.
 * @details When particles are deposited onto the mesh, part of the mass that lands in a cut cell can be moved to its
 * neighbours, e.g. the mass that a hybrid deposition does not keep in the cut cell. This class adds a mass per cut cell
 * to the neighbouring cells, on the cut cell's level and across refinement boundaries, with the stencils of the
 * EBFluxRedistribution operators of each level. It reuses those operators and builds nothing of its own.
 */
class ParticleMassRedistribution
{
public:
  /**
   * @brief Disallowed weak construction.
   */
  ParticleMassRedistribution() = delete;

  /**
   * @brief Constructor.
   * @param[in] a_fluxRedistribution Flux redistribution operators on each level, e.g. from
   * AmrMesh::getRedistributionOp.
   */
  explicit ParticleMassRedistribution(const Vector<RefCountedPtr<EBFluxRedistribution>>& a_fluxRedistribution) noexcept;

  /**
   * @brief Destructor.
   */
  ~ParticleMassRedistribution() noexcept;

  /**
   * @brief Add the mass of each cut cell to its neighbouring cells.
   * @details The mass of a cut cell on level l is added to cells on levels l-1, l and l+1, with the weights of the
   * EBFluxRedistribution stencils. a_phi is a density, so the mass is added as mass per unit volume.
   * @param[in,out] a_phi    Data that the mass is added to.
   * @param[in]     a_deltaM Mass per cut cell, divided by the volume of a regular cell on its level, i.e. in the units
   * of a_phi.
   */
  void
  redistribute(EBAMRCellData& a_phi, const EBAMRIVData& a_deltaM) const noexcept;

protected:
  /**
   * @brief Flux redistribution operators on each level.
   */
  Vector<RefCountedPtr<EBFluxRedistribution>> m_fluxRedistribution;
};

#include <CD_NamespaceFooter.H>

#endif
