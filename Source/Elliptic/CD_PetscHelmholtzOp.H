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

// Std includes
#include <vector>

// Chombo includes
#include <RefCountedPtr.H>
#include <RealVect.H>
#include <LayoutData.H>
#include <BaseFab.H>
#include <BaseIVFAB.H>

// Our includes
#include <CD_PetscGrid.H>
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_EBMultigridInterpolator.H>
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
 * The operator is composite: it spans all AMR levels, and its rows are the valid cells of every level. Supported so
 * far: one phase, no embedded boundaries, no periodic directions. Anything else aborts.
 *
 * The operator is assembled finite-volume style. The row of a cell is its finite-volume balance divided by dx^SpaceDim
 * and kappa-weighted, as in FieldSolverGMG: kappa*alpha*a*phi plus beta/dx times the sum, over the faces of the cell
 * (including the EB face), of the face area fraction (the boundary area on the EB face) times b times s*F. Here F is
 * the face flux, the derivative of phi in the +dir direction on the face, which is a stencil over PETSc rows, and
 * s = -1 on the low side of the cell and +1 on the high side. The right-hand side must use the same weighting: a volume
 * source is multiplied by kappa, and every factor that multiplies the stencil part of a face flux in the matrix also
 * multiplies its boundary data part w_B * g on the right-hand side, see addBoundaryConditionsToRHS. Without embedded
 * boundaries kappa and the area fractions are 1, and are not applied yet. Cells covered by a finer level have no row.
 * The faces are assembled as follows:
 *
 * - Face between two valid cells on the same level: the two-point flux, as on a single level.
 * - Domain face: the stencils and the boundary weight w_B from PetscHelmholtzDomainBC::getFaceFluxStencil. The
 *   stencil has parts on the next coarser level, this level, and the next finer level. The part on this level may
 *   contain coarse-fine ghost cells, which are replaced as on the fine side of a refinement boundary below. Any other
 *   entry that is not a PETSc row, e.g. a cell covered by a finer grid, aborts, as do entries outside the ghost cells
 *   of the PetscGrid maps.
 * - Fine side of a refinement boundary: the two-point flux between the fine cell and the ghost cell across the face.
 *   The ghost cell is not a PETSc row, so it is replaced by its interpolation stencil from
 *   EBMultigridInterpolator::getGhostStencil, which is over cells on the fine level and the next coarser level, and
 *   mapped with PetscGrid::getPetscStencilWithCoarserLevel.
 * - Coarse side of a refinement boundary: the flux through the coarse face is the area-weighted average of the fine
 *   face fluxes that cover it, using the same fine-face stencils as the fine rows. The coarse and fine rows thus see
 *   the same flux through the interface, which makes the composite operator conservative. Each fine face adds its flux
 *   to the coarse cell across it, through the coarsened fine grids (PetscGrid::getAMRToPetscCoFi), with the weight
 *   -1/refRat^SpaceDim relative to the fine row. The coarse cell skips its own faces onto cells covered by the finer
 *   level. The coarse row may be on another rank.
 *
 * The coefficients alpha, beta, a and b are multiplied into the matrix entries, so changing any of them means refilling
 * the matrix. The a- and b-coefficients are held by reference, so the caller can change their values in place;
 * setCoefficients points the operator to other data, and setAlphaAndBeta changes alpha and beta. Neither touches the
 * matrix: fillMatrix is the only function that does. It refills the existing matrix, keeping its preallocation and
 * sparsity pattern, so the AMG preconditioner can keep its interpolation (PCGAMGSetReuseInterpolation) and only
 * recompute its coarse operators.
 *
 * Faces between two valid cells on one level have two-point stencils, which fillMatrix computes as it goes. The other
 * faces have stencils that are expensive to compute, so they are computed once per define and stored without beta and
 * b; fillMatrix then only multiplies them by beta and the b-coefficient on the face. The stored stencils are over PETSc
 * rows, with ghost cells already replaced, and are held per category so that fillMatrix iterates over each category
 * on its own:
 *
 * - Domain faces: m_domainFluxStencils, with their boundary weights in m_boundaryWeights.
 * - EB faces: m_ebFluxStencils.
 * - Fine side of a refinement boundary: m_fineCFFluxStencils, multiplied by the fine-face b.
 * - Coarse side of a refinement boundary: m_coarCFFluxStencils, one stencil per fine face under the coarse face.
 *   Each is multiplied by the b-coefficient on its fine face, which fillMatrix first copies from the fine level into
 *   the index space of the coarse level (m_bCoefFiCo).
 *
 * So far the stencils of the domain faces and of both sides of the refinement boundaries are computed in fillMatrix
 * rather than stored; the holders, and the functions that define them, are placeholders.
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
   * @details This defines the domain BC, creates the matrix, and fills it.
   * The a- and b-coefficients are held by reference, see setCoefficients.
   * @param[in] a_petscGrid     PETSc grid (rows)
   * @param[in] a_interpolators Coarse-fine ghost cell interpolators on each level. Entry 0 is not used.
   * @param[in] a_domainBc      Domain boundary conditions
   * @param[in] a_aCoef         Cell-centered a-coefficient
   * @param[in] a_bCoef         Face-centered b-coefficient
   * @param[in] a_probLo        Lower-left corner of the computational domain
   * @param[in] a_dx            Grid resolution on each level
   * @param[in] a_alpha         Operator alpha
   * @param[in] a_beta          Operator beta
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
   * @brief Set alpha and beta. The matrix is not changed until fillMatrix is called.
   * @param[in] a_alpha Operator alpha
   * @param[in] a_beta  Operator beta
   */
  virtual void
  setAlphaAndBeta(const Real a_alpha, const Real a_beta) noexcept;

  /**
   * @brief Point the operator to other a- and b-coefficients. The matrix is not changed until fillMatrix is called.
   * @details The coefficients are held by reference, so later changes to their values are seen by fillMatrix.
   * @param[in] a_aCoef Cell-centered a-coefficient
   * @param[in] a_bCoef Face-centered b-coefficient
   */
  virtual void
  setCoefficients(const MFAMRCellData& a_aCoef, const MFAMRFluxData& a_bCoef) noexcept;

  /**
   * @brief Fill the matrix from the current alpha, beta, and a- and b-coefficients.
   * @details This refills the existing matrix, keeping its preallocation and sparsity pattern. It also stores the
   * boundary weights used by addBoundaryConditionsToRHS.
   */
  virtual void
  fillMatrix() noexcept;

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
   * subtracts c from a_rhs. It evaluates the boundary data but not the stencils, using the w_B kept from assembly. Only
   * boundary terms are added: the domain boundary data now, and the EB boundary data and surface charge later. The
   * caller supplies the volume source in a_rhs, kappa-weighted as described for the class.
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
   * @brief Coarse-fine ghost cell interpolators on each level. Entry 0 is not used.
   */
  Vector<RefCountedPtr<EBMultigridInterpolator>> m_interpolators;

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
   * @details From PetscHelmholtzDomainBC::getFaceFluxStencil, filled by fillMatrix. Defined on the phase that is
   * assembled.
   */
  EBAMRIFData m_boundaryWeights;

  /**
   * @brief Flux stencils through the domain faces, without beta and b, over PETSc rows.
   * @details Indexed as [dir][side][level], with Side::Lo = 0 and Side::Hi = 1. Each BaseFab is defined over the cells
   * of the box that abut the domain side, and holds the stencil of the flux through the domain face of each cell.
   * Placeholder: not defined or used yet.
   */
  Vector<RefCountedPtr<LayoutData<BaseFab<PetscStencil>>>> m_domainFluxStencils[SpaceDim][2];

  /**
   * @brief Flux stencils through the EB faces, without beta and b, over PETSc rows.
   * @details Indexed by level, with one stencil per cut cell. Placeholder: not defined or used yet.
   */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<PetscStencil>>>> m_ebFluxStencils;

  /**
   * @brief Flux stencils through the fine faces on a refinement boundary, seen from the fine side.
   * @details Without beta and b, over PETSc rows, with the coarse-fine ghost cells replaced. Indexed as
   * [dir][side][level]. Each BaseFab is defined over the cells of the box that abut the box side, and holds the stencil
   * of the flux through that side for the cells whose neighbor across it is a coarse-fine ghost cell; other entries are
   * empty. Placeholder: not defined or used yet.
   */
  Vector<RefCountedPtr<LayoutData<BaseFab<PetscStencil>>>> m_fineCFFluxStencils[SpaceDim][2];

  /**
   * @brief Flux stencils through the fine faces on a refinement boundary, seen from the coarse side.
   * @details Without beta and b, over PETSc rows. Indexed as [dir][side][level], for the levels below the finest. Each
   * BaseFab is defined over the cells of the coarse box. For a coarse cell whose neighbor across the side is covered by
   * the finer level, the entry holds the stencils of the fine faces under the coarse face, in a fixed order that
   * matches the fine faces in m_bCoefFiCo; other entries are empty. Placeholder: not defined or used yet.
   */
  Vector<RefCountedPtr<LayoutData<BaseFab<std::vector<PetscStencil>>>>> m_coarCFFluxStencils[SpaceDim][2];

  /**
   * @brief The b-coefficient of the finer level, in the index space of the coarse level.
   * @details Entry lvl is defined on the grids of level lvl refined to level lvl + 1, i.e.
   * PetscGrid::getMFLevelGridsFiCo()[lvl + 1], which share DataIndex with level lvl, so that the coarse side of a
   * refinement boundary finds the fine-face b under its own grid index. Copied from the finer level by fillMatrix.
   * Placeholder: not defined or used yet.
   */
  EBAMRFluxData m_bCoefFiCo;

  /**
   * @brief Abort unless the grids are supported by this operator: an interpolator on every level above the coarsest, no
   * periodic directions, no embedded boundaries, and no cells in any phase other than the gas phase.
   */
  virtual void
  checkSupported() const noexcept;

  /**
   * @brief Add a stencil to a row of the matrix, with ADD_VALUES.
   * @details Thread safe: MatSetValues is called in a named critical section. The row may be on another rank.
   * @param[in] a_row     Global PETSc row
   * @param[in] a_stencil Stencil over global PETSc rows
   */
  virtual void
  addToMatrix(const PetscInt a_row, const PetscStencil& a_stencil) noexcept;

  /**
   * @brief Replace the coarse-fine ghost cells in a stencil over level a_level cells by their interpolation stencils.
   * @details Each ghost cell entry is removed from a_stencil, and its stencil from
   * EBMultigridInterpolator::getGhostStencil, times the entry's weight, is added to a_stencil (the part on level
   * a_level) and a_stencilCoar (the part on level a_level - 1, in the index space of grid a_din coarsened). Entries
   * that are not ghost cells are kept.
   * @param[in,out] a_stencil     Stencil over cells on level a_level
   * @param[in,out] a_stencilCoar Stencil over cells on level a_level - 1
   * @param[in]     a_level       Grid level
   * @param[in]     a_din         Grid index on the level
   */
  virtual void
  replaceGhostCells(VoFStencil&      a_stencil,
                    VoFStencil&      a_stencilCoar,
                    const int        a_level,
                    const DataIndex& a_din) const noexcept;

  /**
   * @brief Create the matrix, with its preallocation, and allocate the boundary weights.
   */
  virtual void
  createMatrix() noexcept;

  /**
   * @brief Define the flux stencils through the domain faces, m_domainFluxStencils.
   * @details Placeholder: not implemented yet, and aborts.
   */
  virtual void
  defineDomainFluxStencils() noexcept;

  /**
   * @brief Define the flux stencils through the EB faces, m_ebFluxStencils.
   * @details Placeholder: not implemented yet, and aborts.
   */
  virtual void
  defineEBFluxStencils() noexcept;

  /**
   * @brief Define the flux stencils through the refinement boundaries, m_fineCFFluxStencils and m_coarCFFluxStencils.
   * @details Placeholder: not implemented yet, and aborts.
   */
  virtual void
  defineCoarseFineFluxStencils() noexcept;

  /**
   * @brief Copy the b-coefficient of each finer level into the index space of the coarse level, m_bCoefFiCo.
   * @details Placeholder: not implemented yet, and aborts.
   */
  virtual void
  copyFineBcoefToCoarse() noexcept;

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
