/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzOp.H
 * @brief  Declaration of a Helmholtz operator assembled as a PETSc matrix.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZOP_H
#define CD_PETSCHELMHOLTZOP_H

// Petsc includes
#include <petsc.h>

// Chombo includes
#include <RefCountedPtr.H>
#include <RealVect.H>

// Our includes
#include <CD_PetscGrid.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_EBAMRData.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Helmholtz operator alpha*a*phi + beta*div(b*grad(phi)), assembled as a PETSc matrix over the rows of a
 * PetscGrid.
 * @details Unlike EBHelmholtzOp, which applies the operator to data, this class assembles the operator once as a matrix
 * so that it can be solved with PETSc. The domain boundary conditions enter through PetscHelmholtzDomainBC: the
 * phi-dependent part of the boundary flux goes into the matrix and the constant part into the right-hand side, see
 * addBoundaryConditionsToRHS.
 *
 * The operator is composite: it spans all AMR levels, and the rows are the valid cells of each level. Near a
 * refinement boundary:
 * - A fine-level face whose neighbor is a coarse-fine ghost cell uses the ghost cell's interpolation stencil from
 *   EBMultigridInterpolator::getGhostStencil, which is over level l and level l-1 rows.
 * - A coarse cell does not use its faces that open onto cells covered by the finer level. Instead, every fine face on
 *   the refinement boundary adds its flux to the coarse cell on the other side, scaled by the face area ratio. This is
 *   refluxing, built into the matrix, and makes the operator conservative across levels.
 *
 * Supported so far: one phase, no embedded boundaries, no periodic directions. Anything else aborts.
 */
class PetscHelmholtzOp
{
public:
  /**
   * @brief Default constructor. Must subsequently call define.
   */
  PetscHelmholtzOp() noexcept;

  /**
   * @brief Disallowed copy constructor.
   * @param[in] a_other Other operator
   */
  PetscHelmholtzOp(const PetscHelmholtzOp& a_other) = delete;

  /**
   * @brief Disallowed move constructor.
   * @param[in] a_other Other operator
   */
  PetscHelmholtzOp(const PetscHelmholtzOp&& a_other) = delete;

  /**
   * @brief Disallowed copy assignment operator.
   * @param[in] a_other Other operator
   * @return This operator.
   */
  PetscHelmholtzOp&
  operator=(const PetscHelmholtzOp& a_other) = delete;

  /**
   * @brief Disallowed move assignment operator.
   * @param[in] a_other Other operator
   * @return This operator.
   */
  PetscHelmholtzOp&
  operator=(const PetscHelmholtzOp&& a_other) = delete;

  /**
   * @brief Destructor. Destroys the matrix.
   */
  virtual ~PetscHelmholtzOp() noexcept;

  /**
   * @brief Define the operator.
   * @details This defines the domain BC on all levels and assembles the matrix.
   * @param[in] a_petscGrid     PETSc grid (rows)
   * @param[in] a_interpolators Coarse-fine ghost cell interpolators on each level. Entry 0 is not used.
   * @param[in] a_domainBc      Domain boundary conditions
   * @param[in] a_aCoef     Cell-centered a-coefficient
   * @param[in] a_bCoef     Face-centered b-coefficient
   * @param[in] a_probLo    Lower-left corner of the computational domain
   * @param[in] a_dx        Grid resolution on each level
   * @param[in] a_alpha     Operator alpha
   * @param[in] a_beta      Operator beta
   */
  virtual void
  define(const RefCountedPtr<PetscGrid>&                       a_petscGrid,
         const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_interpolators,
         const RefCountedPtr<PetscHelmholtzDomainBC>&          a_domainBc,
         const MFAMRCellData&                                  a_aCoef,
         const MFAMRFluxData&                                  a_bCoef,
         const RealVect&                                       a_probLo,
         const Vector<Real>&                                   a_dx,
         const Real                                            a_alpha,
         const Real                                            a_beta) noexcept;

  /**
   * @brief Get the assembled matrix.
   * @return m_matrix
   */
  virtual const Mat&
  getMatrix() const noexcept;

  /**
   * @brief Move the constant part of the domain boundary fluxes into a right-hand side.
   * @details The operator is A*phi + g, where A is the matrix and g holds the constant parts of the boundary fluxes.
   * Solving A*phi + g = rhs is the same as solving A*phi = rhs - g, so this subtracts g from a_rhs.
   * @param[in,out] a_rhs Right-hand side, on the rows of the PETSc grid
   */
  virtual void
  addBoundaryConditionsToRHS(Vec& a_rhs) const noexcept;

protected:
  /**
   * @brief Is defined or not.
   */
  bool m_isDefined;

  /**
   * @brief Verbose output or not.
   */
  bool m_verbose;

  /**
   * @brief PETSc grid.
   */
  RefCountedPtr<PetscGrid> m_petscGrid;

  /**
   * @brief Coarse-fine ghost cell interpolators on each level.
   */
  Vector<RefCountedPtr<EBMultigridInterpolator>> m_interpolators;

  /**
   * @brief Domain boundary conditions.
   */
  RefCountedPtr<PetscHelmholtzDomainBC> m_domainBc;

  /**
   * @brief Cell-centered a-coefficient.
   */
  MFAMRCellData m_aCoef;

  /**
   * @brief Face-centered b-coefficient.
   */
  MFAMRFluxData m_bCoef;

  /**
   * @brief Lower-left corner of the computational domain.
   */
  RealVect m_probLo;

  /**
   * @brief Grid resolution on each level.
   */
  Vector<Real> m_dx;

  /**
   * @brief Operator alpha.
   */
  Real m_alpha;

  /**
   * @brief Operator beta.
   */
  Real m_beta;

  /**
   * @brief The assembled matrix.
   */
  Mat m_matrix;

  /**
   * @brief Abort unless the grids are supported by this operator: one phase, no embedded boundaries, no periodic
   * directions.
   */
  virtual void
  checkSupported() const noexcept;

  /**
   * @brief Assemble the matrix.
   */
  virtual void
  assembleMatrix() noexcept;

  /**
   * @brief Add a stencil to a row of the matrix.
   * @param[in] a_row     Global PETSc row
   * @param[in] a_stencil Stencil over global PETSc rows
   */
  virtual void
  addToMatrix(const PetscInt a_row, const PetscStencil& a_stencil) noexcept;

  /**
   * @brief Get the b-coefficient on a face.
   * @param[in] a_iv    Cell
   * @param[in] a_dir   Coordinate direction of the face
   * @param[in] a_side  Side of the cell the face is on
   * @param[in] a_level Level
   * @param[in] a_din   Grid index
   * @param[in] a_phase Phase
   * @return The b-coefficient on the face.
   */
  virtual Real
  getBcoef(const IntVect&       a_iv,
           const int            a_dir,
           const Side::LoHiSide a_side,
           const int            a_level,
           const DataIndex&     a_din,
           const int            a_phase) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
