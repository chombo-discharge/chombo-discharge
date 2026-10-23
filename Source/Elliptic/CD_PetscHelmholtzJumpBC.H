/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzJumpBC.H
 * @brief  Declaration of a base class for the interface (jump) conditions of PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZJUMPBC_H
#define CD_PETSCHELMHOLTZJUMPBC_H

// Petsc includes
#include <petsc.h>

// Our includes
#include <CD_PetscHelmholtzEBStencils.H>
#include <CD_EBAMRData.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Base class for the conditions on the interface between two phases, for PetscHelmholtzOp.
 * @details An interface cell is a cut cell where both phases have a VoF. On each side k of the interface, the one-sided
 * normal derivative from PetscHelmholtzEBStencils is D_k = w_B,k * phi_B + sum_i w_k,i * phi_k,i over the valid cells
 * of phase k. A subclass supplies the closure: one equation in the D_k that eliminates phi_B, so that the flux into
 * each side is a linear form over the PETSc rows of both phases and the jump data g,
 *
 *   b_k * D_k = w_g,k * g + sum_j w_k,j * phi_j.
 *
 * Unlike the electrode and domain boundary conditions, this flux contains the b-coefficient: the closure combines the
 * b-coefficients of both phases, so the weights depend on them. PetscHelmholtzOp passes them in, and multiplies the
 * flux by the face area only. After a change of the b-coefficients, the fluxes must be fetched again; the least-squares
 * stencils are stored and are not recomputed.
 */
class PetscHelmholtzJumpBC : public PetscHelmholtzEBStencils
{
public:
  /**
   * @brief Default constructor. Must subsequently set the order and the weight, and call define.
   */
  PetscHelmholtzJumpBC() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzJumpBC() noexcept;

  /**
   * @brief Get the flux into one side of an interface cell, b_k * D_k, as a linear form.
   * @param[out] a_dataWeight   Weight w_g of the jump data, see getJumpData
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
                          const Real       a_bCoefPhase1) const noexcept = 0;

  /**
   * @brief Get the jump data g at an interface cell.
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_iv    Interface cell
   * @param[in] a_sigma Surface charge density on the interface
   * @return The jump data g.
   */
  virtual Real
  getJumpData(const int          a_level,
              const DataIndex&   a_din,
              const IntVect&     a_iv,
              const EBAMRIVData& a_sigma) const noexcept = 0;

  /**
   * @brief Update the jump data from the solution, after the solve.
   * @details For closures where the jump data is an output (e.g. the saturation charge), this computes it from the
   * solution with the same stencils and keeps it in the boundary condition object, as MFHelmholtzJumpBC does; the
   * solver's surface charge density is not changed. For closures where the jump data is an input, this does nothing.
   * @param[in]     a_phi         Solution, on the rows of the PETSc grid
   * @param[in]     a_bCoefIrreg  b-coefficient on the EB faces
   */
  virtual void
  updateJumpData(const Vec& a_phi, const MFAMRIVData& a_bCoefIrreg) noexcept = 0;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
