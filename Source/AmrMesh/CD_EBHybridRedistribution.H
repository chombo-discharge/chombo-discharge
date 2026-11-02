/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBHybridRedistribution.H
 * @brief  Declaration of the hybrid divergence with flux redistribution as an EBRedistribution.
 * @author Robert Marskar
 */

#ifndef CD_EBHYBRIDREDISTRIBUTION_H
#define CD_EBHYBRIDREDISTRIBUTION_H

// Std includes
#include <string>

// Our includes
#include <CD_EBRedistribution.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief The hybrid divergence with flux redistribution, as an EBRedistribution.
 * @details In a cut cell with volume fraction kappa, the right-hand side is the hybrid divergence
 *
 *    R = kappa*D + (1-kappa)*D_nc,
 *
 * where D is the conservative divergence and D_nc the non-conservative divergence, a kappa-weighted average of kappa*D
 * over the neighbouring cells. The mass that this misses, (1-kappa)*(kappa*D - kappa*D_nc), is redistributed to the
 * neighbouring cells through the EBFluxRedistribution operators of the AmrMesh, on this level and across refinement
 * boundaries. The result depends on neither the state nor the time step.
 */
class EBHybridRedistribution : public EBRedistribution
{
public:
  /**
   * @brief Disallowed weak construction.
   */
  EBHybridRedistribution() = delete;

  /**
   * @brief Constructor. Allocates the scratch data on the current grids, so the object must be rebuilt after a regrid.
   * @details The realm must have the operators for the non-conservative divergence and, with a_redistribute, for the
   * redistribution registered.
   * @param[in] a_amr                AMR mesh
   * @param[in] a_realm              Realm
   * @param[in] a_phase              Phase
   * @param[in] a_blendConservation  If false, the non-conservative divergence is zero, and R = kappa*D in cut cells.
   * @param[in] a_redistribute       If false, the missed mass is not redistributed, and the update is not
   * conservative.
   */
  EBHybridRedistribution(const RefCountedPtr<AmrMesh>& a_amr,
                         const std::string&            a_realm,
                         const phase::which_phase      a_phase,
                         const bool                    a_blendConservation,
                         const bool                    a_redistribute) noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~EBHybridRedistribution() noexcept;

  /**
   * @brief Turn kappa*D into the hybrid divergence and redistribute the missed mass.
   * @details See EBRedistribution::redistribute. The state and the time step are not used.
   * @param[in,out] a_divG On input kappa*D, on output R.
   * @param[in]     a_phi  State that the update starts from (not used).
   * @param[in]     a_dt   Time step of the update (not used).
   */
  void
  redistribute(EBAMRCellData& a_divG, const EBAMRCellData& a_phi, const Real a_dt) noexcept override;

protected:
  /**
   * @brief AMR mesh
   */
  RefCountedPtr<AmrMesh> m_amr;

  /**
   * @brief Realm
   */
  std::string m_realm;

  /**
   * @brief Phase
   */
  phase::which_phase m_phase;

  /**
   * @brief If false, the non-conservative divergence is zero.
   */
  bool m_blendConservation;

  /**
   * @brief If false, the missed mass is not redistributed.
   */
  bool m_redistribute;

  /**
   * @brief Non-conservative divergence in the cut cells.
   */
  EBAMRIVData m_nonConservativeDivergence;

  /**
   * @brief Mass that the hybrid divergence misses in the cut cells.
   */
  EBAMRIVData m_massDifference;

  /**
   * @brief Make the hybrid divergence on one level.
   * @details On input, a_hybridDivergence holds kappa*D. On output, it holds the hybrid divergence in the cut cells,
   * and a_massDifference the mass that it misses.
   * @param[in,out] a_hybridDivergence          On input kappa*D, on output the hybrid divergence.
   * @param[out]    a_massDifference            Mass that the hybrid divergence misses.
   * @param[in]     a_nonConservativeDivergence Non-conservative divergence.
   * @param[in]     a_level                     Grid level
   */
  void
  hybridDivergence(LevelData<EBCellFAB>&             a_hybridDivergence,
                   LevelData<BaseIVFAB<Real>>&       a_massDifference,
                   const LevelData<BaseIVFAB<Real>>& a_nonConservativeDivergence,
                   const int                         a_level) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
