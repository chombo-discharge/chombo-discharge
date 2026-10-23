/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzEBBC.H
 * @brief  Declaration of a base class for electrode (EB) boundary conditions of PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZEBBC_H
#define CD_PETSCHELMHOLTZEBBC_H

// Our includes
#include <CD_PetscHelmholtzEBStencils.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Base class for the boundary conditions on the EB faces of cut cells that are not interface cells (electrodes).
 * @details The flux through the EB face of a cut cell is returned as a linear form in the cell values and the boundary
 * data g,
 *
 *   F = w_B * g + sum_j w_j * phi_j,
 *
 * where F is the normal derivative of phi with n pointing out of the cell, and the stencil is over PETSc rows, as in
 * PetscHelmholtzEBStencils. F does not contain the b-coefficient or the face area, which PetscHelmholtzOp multiplies
 * in. The stencil goes into the matrix, and w_B * g into the right-hand side. This mirrors PetscHelmholtzDomainBC for
 * the domain faces.
 */
class PetscHelmholtzEBBC : public PetscHelmholtzEBStencils
{
public:
  /**
   * @brief Default constructor. Must subsequently set the order and the weight, and call define.
   */
  PetscHelmholtzEBBC() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzEBBC() noexcept;

  /**
   * @brief Get the flux through the EB face of a cut cell, without the b-coefficient and the face area.
   * @param[out] a_boundaryWeight Weight w_B of the boundary data, see getBoundaryData
   * @param[out] a_stencil        Stencil over PETSc rows
   * @param[in]  a_phase          Phase
   * @param[in]  a_level          Grid level
   * @param[in]  a_din            Grid index on the level
   * @param[in]  a_vof            Cut cell
   */
  virtual void
  getEBFluxStencil(Real&            a_boundaryWeight,
                   PetscStencil&    a_stencil,
                   const int        a_phase,
                   const int        a_level,
                   const DataIndex& a_din,
                   const VolIndex&  a_vof) const noexcept = 0;

  /**
   * @brief Get the boundary data g at the EB face of a cut cell.
   * @details The flux depends on g through w_B * g, with w_B from getEBFluxStencil. g is evaluated whenever the
   * right-hand side is built.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cut cell
   * @return The boundary data g.
   */
  virtual Real
  getBoundaryData(const int        a_phase,
                  const int        a_level,
                  const DataIndex& a_din,
                  const VolIndex&  a_vof) const noexcept = 0;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
