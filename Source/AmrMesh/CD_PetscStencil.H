/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscStencil.H
 * @brief  Declaration of a simple stencil class for PETSc.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCSTENCIL_H
#define CD_PETSCSTENCIL_H

// Petsc includes
#include "petsc.h"

// Chombo includes
#include <REAL.H>
#include <Vector.H>

// Our includes
#include <CD_NamespaceHeader.H>

class PetscStencil
{
public:
  /**
   * @brief Default constructor. Creates an empty stencil.
   */
  PetscStencil();

  /**
   * @brief Destructor.
   */
  ~PetscStencil();

  /**
   * @brief Clear the stencil
   */
  void
  clear();

  /**
   * @brief Get the size of the stencil
   */
  int
  size() const noexcept;

  /**
   * @brief Get the PETSc row for the input index.
   * @param[in] a_entry Index in the m_rows vector.
   */
  PetscInt
  row(const int a_entry) const noexcept;

  /**
   * @brief Get the PETSc row for the input index.
   * @param[in] a_entry Index in the m_rows vector.
   */
  PetscScalar
  weight(const int a_entry) const noexcept;

  /**
   * @brief Add a row to this stencil
   * @param[in] a_row PETSc row
   * @param[in] a_weight PETSc weight
   */
  void
  add(const PetscInt& a_row, const PetscScalar& a_weight) noexcept;

  /**
   * @brief Add another stencil to this stencil.
   * @param[in] a_otherStencil Other PETSc stencil
   */
  PetscStencil&
  operator+=(const PetscStencil& a_otherStencil) noexcept;

  /**
   * @brief Scale this stencil.
   * @param[in] a_scalingFactor Scaling factor
   * @note This does not scale the constant term.
   */
  PetscStencil&
  operator*=(const Real& a_scalingFactor) noexcept;

protected:
  /**
   * @brief PETSc rows
   */
  Vector<PetscInt> m_rows;

  /**
   * @brief PETSc weights.
   */
  Vector<PetscScalar> m_weights;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
