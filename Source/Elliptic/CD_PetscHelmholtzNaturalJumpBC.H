/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzNaturalJumpBC.H
 * @brief  Declaration of the natural jump condition on the interface between two phases, for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZNATURALJUMPBC_H
#define CD_PETSCHELMHOLTZNATURALJUMPBC_H

// Our includes
#include <CD_PetscHelmholtzJumpBC.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief The natural jump condition b_0 * D_0 + b_1 * D_1 = sigma on the interface between two phases.
 * @details The closure eliminates phi_B with phi_B = (sigma - b_0 sum_i w_0,i phi_0,i - b_1 sum_i w_1,i phi_1,i) /
 * (b_0 w_B,0 + b_1 w_B,1), so the flux into each side couples the rows of both phases, and g = sigma. The jump data is
 * an input, so updateJumpData does nothing.
 */
class PetscHelmholtzNaturalJumpBC : public PetscHelmholtzJumpBC
{
public:
  /**
   * @brief Default constructor. Must subsequently set the order and the weight, and call define.
   */
  PetscHelmholtzNaturalJumpBC() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzNaturalJumpBC() noexcept;

  /**
   * @brief Get the flux into one side of an interface cell, b_k * D_k, as a linear form.
   * @param[out] a_dataWeight   Weight w_g of the jump data
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
   * @brief Get the jump data, i.e. the surface charge density.
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_iv    Interface cell
   * @param[in] a_sigma Surface charge density on the interface
   * @return sigma at the interface cell.
   */
  virtual Real
  getJumpData(const int          a_level,
              const DataIndex&   a_din,
              const IntVect&     a_iv,
              const EBAMRIVData& a_sigma) const noexcept override;

  /**
   * @brief Does nothing: the surface charge density is an input of this closure, and nothing is kept.
   * @param[in]     a_phi         Solution, on the rows of the PETSc grid
   * @param[in]     a_bCoefIrreg  b-coefficient on the EB faces
   */
  virtual void
  updateJumpData(const Vec& a_phi, const MFAMRIVData& a_bCoefIrreg) noexcept override;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
