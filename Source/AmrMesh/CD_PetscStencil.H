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

/**
 * @brief A linear stencil over global PETSc rows, i.e. a list of (row, weight) pairs.
 * @details This plays the role that VoFStencil plays for Chombo data, but each entry is a global PETSc row rather than
 * a VolIndex. A VoFStencil entry is only meaningful together with the level and grid patch it was built on, whereas
 * a PetscStencil entry can be passed directly to MatSetValues. Rows are unique within a stencil: adding an existing
 * row accumulates its weight.
 */
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
   * @brief Remove all entries from the stencil.
   */
  void
  clear();

  /**
   * @brief Get the number of entries in the stencil.
   * @return Number of (row, weight) entries.
   */
  int
  size() const noexcept;

  /**
   * @brief Get the PETSc row of a stencil entry.
   * @param[in] a_entry Entry index, in the range [0, size()).
   * @return Global PETSc row of the entry.
   */
  PetscInt
  row(const int a_entry) const noexcept;

  /**
   * @brief Get the weight of a stencil entry.
   * @param[in] a_entry Entry index, in the range [0, size()).
   * @return Weight of the entry.
   */
  PetscScalar
  weight(const int a_entry) const noexcept;

  /**
   * @brief Add a weight to a row of this stencil.
   * @details If the row is already in the stencil its weight is incremented, otherwise a new entry is appended.
   * Finding the row is a linear scan, so building a stencil of N entries costs O(N^2). This is deliberate: the
   * stencils assembled here have tens of entries, and at most on the order of a hundred (a 3D least-squares
   * neighbourhood expanded by coarse-fine interpolation stencils). At that size a scan over contiguous storage is
   * cheaper than hashing or sorting.
   * @param[in] a_row    Global PETSc row
   * @param[in] a_weight Weight to add to the row
   */
  void
  add(const PetscInt& a_row, const PetscScalar& a_weight) noexcept;

  /**
   * @brief Add another stencil to this stencil, entry by entry, using add().
   * @param[in] a_otherStencil Other PETSc stencil
   * @return This stencil.
   */
  PetscStencil&
  operator+=(const PetscStencil& a_otherStencil) noexcept;

  /**
   * @brief Scale all weights in this stencil.
   * @param[in] a_scalingFactor Scaling factor
   * @return This stencil.
   */
  PetscStencil&
  operator*=(const Real& a_scalingFactor) noexcept;

protected:
  /**
   * @brief Global PETSc rows. Unique within the stencil; m_rows[i] pairs with m_weights[i].
   */
  Vector<PetscInt> m_rows;

  /**
   * @brief Weights, one per entry in m_rows.
   */
  Vector<PetscScalar> m_weights;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
