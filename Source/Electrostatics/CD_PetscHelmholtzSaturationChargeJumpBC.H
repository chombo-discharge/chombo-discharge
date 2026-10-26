/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzSaturationChargeJumpBC.H
 * @brief  Declaration of the saturation charge condition on the interface between two phases, for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZSATURATIONCHARGEJUMPBC_H
#define CD_PETSCHELMHOLTZSATURATIONCHARGEJUMPBC_H

// Our includes
#include <CD_PetscHelmholtzJumpBC.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief The saturation charge condition D_p = 0 on the interface, on the side of one phase p, for PetscHelmholtzOp.
 * @details The closure eliminates phi_B with phi_B = -sum_i w_p,i phi_p,i / w_B,p, so the flux into side p is zero, and
 * the flux into the other side couples the rows of both phases. There is no jump data. The surface charge density is
 * an output: updateJumpData computes it from the solution, sigma = b_0 * D_0 + b_1 * D_1, and keeps it in m_sigma, as
 * MFHelmholtzSaturationChargeJumpBC does.
 */
class PetscHelmholtzSaturationChargeJumpBC : public PetscHelmholtzJumpBC
{
public:
  /**
   * @brief Disallowed default constructor.
   */
  PetscHelmholtzSaturationChargeJumpBC() = delete;

  /**
   * @brief Full constructor. Must subsequently set the order and the weight, and call define.
   * @param[in] a_phase Phase on whose side the normal derivative is zero
   */
  PetscHelmholtzSaturationChargeJumpBC(const phase::which_phase a_phase) noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzSaturationChargeJumpBC() noexcept;

  /**
   * @brief Release the computed surface charge density, the stored stencils and the grids.
   */
  virtual void
  clear() noexcept override;

  /**
   * @brief Get the flux into one side of an interface cell, b_k * D_k, as a linear form.
   * @param[out] a_dataWeight   Weight w_g of the jump data, which is zero
   * @param[out] a_stencil      Stencil over PETSc rows of both phases
   * @param[in]  a_phase        Phase k of the side
   * @param[in]  a_level        Grid level
   * @param[in]  a_din          Grid index on the level
   * @param[in]  a_iv           Interface cell
   * @param[in]  a_bCoefPhase0  b-coefficient on the interface, on the side of phase 0
   * @param[in]  a_bCoefPhase1  b-coefficient on the interface, on the side of phase 1
   */
  virtual void
  getInterfaceFluxStencil(Real&            a_dataWeight,
                          PetscStencil&    a_stencil,
                          const int        a_phase,
                          const int        a_level,
                          const DataIndex& a_din,
                          const IntVect&   a_iv,
                          const Real       a_bCoefPhase0,
                          const Real       a_bCoefPhase1) const noexcept override;

  /**
   * @brief Get the jump data, which is zero: this closure has none.
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_iv    Interface cell
   * @param[in] a_sigma Surface charge density on the interface
   * @return Zero.
   */
  virtual Real
  getJumpData(const int          a_level,
              const DataIndex&   a_din,
              const IntVect&     a_iv,
              const EBAMRIVData& a_sigma) const noexcept override;

  /**
   * @brief Compute the surface charge density from the solution, sigma = b_0 * D_0 + b_1 * D_1, into m_sigma.
   * @param[in]     a_phi         Solution, on the rows of the PETSc grid
   * @param[in]     a_bCoefIrreg  b-coefficient on the EB faces
   */
  virtual void
  updateJumpData(const Vec& a_phi, const MFAMRIVData& a_bCoefIrreg) noexcept override;

  /**
   * @brief Get the surface charge density computed by the last call to updateJumpData.
   * @details In the units of the jump condition, i.e. b * dphi/dn. Defined on the cut cells of the first phase, and
   * zero on those that are not interface cells. Empty until updateJumpData has been called.
   * @return m_sigma
   */
  virtual const EBAMRIVData&
  getSigma() const noexcept;

protected:
  /**
   * @brief Phase on whose side the normal derivative is zero.
   */
  phase::which_phase m_phase;

  /**
   * @brief The surface charge density computed from the solution by updateJumpData.
   * @details Defined on the interface cells, per level.
   */
  EBAMRIVData m_sigma;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
