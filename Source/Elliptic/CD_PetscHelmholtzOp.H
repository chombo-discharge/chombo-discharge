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
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_EBAMRData.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Helmholtz operator alpha*a*phi + beta*div(b*grad(phi)), assembled as a PETSc matrix over the rows of a
 * PetscGrid.
 * @details Unlike EBHelmholtzOp, which applies the operator to data, this class assembles the operator once as a matrix
 * so that it can be solved with PETSc. The domain boundary conditions enter through PetscHelmholtzDomainBC, which gives
 * the flux through a domain face as w_B * g + sum_i w_i * phi_i with g the boundary data. The w_i go into the matrix.
 * The w_B are kept per domain face when the matrix is assembled, and w_B * g goes into the right-hand side, see
 * addBoundaryConditionsToRHS. The boundary data can thus change without the matrix being reassembled.
 *
 * Supported so far: a single level, no embedded boundaries, no periodic directions. Anything else aborts.
 *
 * The operator is assembled finite-volume style: the row of a cell is alpha*a*phi plus beta/dx times the sum, over the
 * faces of the cell, of the face area fraction times b times the face flux, where each face flux is a stencil over
 * PETSc rows. Cells covered by a finer level have no row. Multilevel assembly is not implemented yet; it follows this
 * scheme:
 *
 * - Face between two valid cells on the same level: the two-point flux, as on a single level.
 * - Domain face: the stencil and constant from PetscHelmholtzDomainBC. The stencil has parts on the next coarser
 *   level, this level, and the next finer level. The part on this level may contain coarse-fine ghost cells, which
 *   are replaced as on the fine side of a refinement boundary below. Any other entry that is not a PETSc row, e.g. a
 *   cell covered by a finer grid, aborts, as do entries outside the ghost cells of the PetscGrid maps.
 * - Fine side of a refinement boundary: the two-point flux between the fine cell and the ghost cell across the face.
 *   The ghost cell is not a PETSc row, so it is replaced by its interpolation stencil from
 *   EBMultigridInterpolator::getGhostStencil, which is over cells on the fine level and the next coarser level, and
 *   mapped with PetscGrid::getPetscStencilWithCoarserLevel.
 * - Coarse side of a refinement boundary: the flux through the coarse face is the area-weighted average of the fine
 *   face fluxes that cover it, using the same fine-face stencils as the fine rows. The coarse and fine rows thus see
 *   the same flux through the interface, which makes the composite operator conservative.
 *
 * The coefficients alpha, beta, a and b are multiplied into the matrix entries at assembly, so changing any of them
 * means refilling the matrix (beta alone could be applied with MatScale). Refilling is not implemented yet; it follows
 * this plan:
 *
 * - Refill the existing matrix: MatZeroEntries and the same assembly, which keeps the preallocation and the sparsity
 *   pattern. Every stencil is recomputed, so the result is the same as a fresh assembly. The AMG preconditioner keeps
 *   its interpolation (PCGAMGSetReuseInterpolation) and only recomputes its coarse operators.
 * - If recomputing the stencils dominates, store the face flux stencils without beta and b, one per face, for the faces
 *   whose stencils are expensive: EB faces, domain faces, and the fine faces on refinement boundaries (with the ghost
 *   cells already replaced). The coarse side of a refinement boundary reuses the stencils of the fine faces that
 *   cover it, each multiplied by its own fine-face b. Regular faces between valid cells on one level are two-point
 *   differences and are recomputed.
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
   * @details This defines the domain BC on the level and assembles the matrix.
   * @param[in] a_petscGrid PETSc grid (rows)
   * @param[in] a_domainBc  Domain boundary conditions
   * @param[in] a_aCoef     Cell-centered a-coefficient
   * @param[in] a_bCoef     Face-centered b-coefficient
   * @param[in] a_probLo    Lower-left corner of the computational domain
   * @param[in] a_dx        Grid resolution on each level
   * @param[in] a_alpha     Operator alpha
   * @param[in] a_beta      Operator beta
   */
  virtual void
  define(const RefCountedPtr<PetscGrid>&              a_petscGrid,
         const RefCountedPtr<PetscHelmholtzDomainBC>& a_domainBc,
         const MFAMRCellData&                         a_aCoef,
         const MFAMRFluxData&                         a_bCoef,
         const RealVect&                              a_probLo,
         const Vector<Real>&                          a_dx,
         const Real                                   a_alpha,
         const Real                                   a_beta) noexcept;

  /**
   * @brief Get the assembled matrix.
   * @return m_matrix
   */
  virtual const Mat&
  getMatrix() const noexcept;

  /**
   * @brief Move the boundary data part of the domain boundary fluxes into a right-hand side.
   * @details The operator is A*phi + c, where A is the matrix and c holds the w_B * g parts of the boundary fluxes,
   * scaled like the rest of the face flux. Solving A*phi + c = rhs is the same as solving A*phi = rhs - c, so this
   * subtracts c from a_rhs. It evaluates the boundary data but not the stencils, using the w_B kept from assembly.
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
   * @brief Weight w_B of the boundary data in the flux through each domain face.
   * @details From PetscHelmholtzDomainBC::getFaceFluxStencil, filled when the matrix is assembled. Defined on the phase
   * that is assembled.
   */
  EBAMRIFData m_boundaryWeights;

  /**
   * @brief Abort unless the grids are supported by this operator: one level, no embedded boundaries.
   */
  virtual void
  checkSupported() const noexcept;

  /**
   * @brief Assemble the matrix.
   */
  virtual void
  assembleMatrix() noexcept;

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
