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

// Std includes
#include <utility>
#include <vector>

// Petsc includes
#include "petsc.h"

// Our includes
#include <CD_NamespaceHeader.H>

/**
 * @brief A linear stencil over global PETSc rows, i.e. a list of (row, weight) pairs.
 * @details This plays the role that VoFStencil plays for Chombo data, but each entry is a global PETSc row rather than
 * a VolIndex. A VoFStencil entry is only meaningful together with the level and grid patch it was built on, whereas
 * a PetscStencil entry can be passed directly to MatSetValues.
 *
 * Entries are appended by add() and operator+= without looking for an existing entry with the same row, so a row
 * can appear more than once until compress() is called. compress() sorts the entries by row and merges duplicates,
 * which makes the rows unique and ascending. Duplicates are summed in the order they were added, so the result is the
 * same as merging on every add, bit for bit. Call compress() before passing the entries to MatSetValues: duplicate
 * rows would be summed with ADD_VALUES but would overwrite each other with INSERT_VALUES. Building a stencil of N
 * entries this way costs O(N log N), compared with O(N^2) for merging on every add; at a hundred entries the difference
 * is roughly a factor of three to five.
 */
class PetscStencil
{
public:
  /**
   * @brief Default constructor. Creates an empty stencil.
   */
  PetscStencil() = default;

  /**
   * @brief Destructor.
   */
  ~PetscStencil() = default;

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
   * @brief Append an entry to this stencil.
   * @details The entry is not merged with an existing entry for the same row; call compress() for that.
   * @param[in] a_row    Global PETSc row
   * @param[in] a_weight Weight of the entry
   */
  void
  add(const PetscInt& a_row, const PetscScalar& a_weight) noexcept;

  /**
   * @brief Append all entries of another stencil to this stencil.
   * @details Entries are not merged with existing entries for the same row; call compress() for that.
   * @param[in] a_otherStencil Other PETSc stencil
   * @return This stencil.
   */
  PetscStencil&
  operator+=(const PetscStencil& a_otherStencil) noexcept;

  /**
   * @brief Sort the entries by row and merge entries with the same row by summing their weights.
   * @details After this the rows are unique and in ascending order, which is the cheapest order for MatSetValues. The
   * sort is stable, so the weights of a row are summed in the order they were added.
   */
  void
  compress() noexcept;

  /**
   * @brief Scale all weights in this stencil.
   * @param[in] a_scalingFactor Scaling factor
   * @return This stencil.
   */
  PetscStencil&
  operator*=(const PetscScalar& a_scalingFactor) noexcept;

protected:
  /**
   * @brief Stencil entries as (global PETSc row, weight). Rows are unique and ascending only after compress().
   */
  std::vector<std::pair<PetscInt, PetscScalar>> m_entries;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
