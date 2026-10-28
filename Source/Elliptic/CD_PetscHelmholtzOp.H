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

// Our includes
#include <CD_PetscGrid.H>
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_PetscHelmholtzEBBC.H>
#include <CD_PetscHelmholtzJumpBC.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_MFMultigridInterpolator.H>
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
 * The operator is composite: it spans all AMR levels, and its rows are the valid cells of every level, in each phase.
 * Embedded boundaries (electrodes) and the interface between two phases are supported; periodic directions are not,
 * and abort.
 *
 * The operator is assembled finite-volume style. The row of a cell is its finite-volume balance divided by dx^SpaceDim
 * and kappa-weighted, as in FieldSolverGMG, and then multiplied by the row weight (dx/dx_0)^SpaceDim of its level (see
 * getRowWeight): kappa*alpha*a*phi plus beta/dx times the sum, over the faces of the cell (including the EB face), of
 * the face area fraction (the boundary area on the EB face) times b times s*F. Here F is the face flux, the derivative
 * of phi in the +dir direction on the face, which is a stencil over PETSc rows, and s = -1 on the low side of the cell
 * and +1 on the high side. The right-hand side must use the same weighting: a volume source is multiplied by kappa,
 * every factor that multiplies the stencil part of a face flux in the matrix also multiplies its boundary data part
 * w_B * g on the right-hand side (see addBoundaryConditionsToRHS), and the whole right-hand side is then multiplied by
 * the row weights (scaleRows; addBoundaryConditionsToRHS does not apply them). Cells covered by a finer level have no
 * row. The faces are assembled as follows:
 *
 * - Face between two valid cells on the same level: the two-point flux, as on a single level.
 * - Face of a cut cell: as in EBHelmholtzOp, the two-point fluxes b * dphi/dx through the face and its neighbouring
 *   face are interpolated to the face centroid (EBArith::getInterpStencil), with b on each face. The interpolation
 *   falls back to the face itself if it would use a cell covered by the finer level. Coarse-fine ghost cells in it are
 *   replaced by their interpolation stencils.
 * - Domain face of a cut cell: the domain BC fluxes of the face and its neighbouring domain face, interpolated to the
 *   face centroid in the same way, as EBHelmholtzDomainBC does. On the right-hand side the boundary data of each of
 *   these faces enters with its own w_B.
 * - EB face: the flux from PetscHelmholtzEBBC::getEBFluxStencil, times the boundary area and the b-coefficient on the
 *   EB, and w_B * g on the right-hand side. On an interface cell, which is cut in both phases, the flux comes from
 *   PetscHelmholtzJumpBC::getInterfaceFluxStencil instead; it already contains the b-coefficients, so it is only
 *   multiplied by the boundary area, and w_g * g goes on the right-hand side.
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
 *   -1/refRat^SpaceDim relative to the fine row before the row weights, and -1 after them. The coarse cell skips
 *   its own faces onto cells covered by the finer level. The coarse row may be on another rank.
 *
 * The coefficients alpha, beta, a and b are multiplied into the matrix entries, so changing any of them means refilling
 * the matrix. The a- and b-coefficients are held by reference, so the caller can change their values in place;
 * setCoefficients points the operator to other data, and setAlphaAndBeta changes alpha and beta. Neither touches the
 * matrix: fillMatrix is the only function that does. It refills the existing matrix, keeping its preallocation and
 * sparsity pattern, so the AMG preconditioner can keep its interpolation (PCGAMGSetReuseInterpolation) and only
 * recompute its coarse operators.
 *
 * After any change of alpha, beta, a or b, the caller must call fillMatrix and then rebuild the right-hand side
 * (addBoundaryConditionsToRHS) before solving. The stencils can be reused, but the weights must be redone, and both the
 * matrix and the right-hand side carry them. addBoundaryConditionsToRHS uses the current beta and b, so a right-hand
 * side built after a coefficient change, against a matrix filled before it, gives a matrix and a right-hand side that
 * do not belong together. That is a user error, and it is not detected.
 *
 * Faces between two valid cells on one level have two-point stencils, which fillMatrix computes as it goes. The other
 * faces have stencils that are expensive to compute, so they are computed once per define and stored without beta and
 * b; fillMatrix then only multiplies them by beta and the b-coefficient on the face. The stored stencils are over PETSc
 * rows, with ghost cells already replaced, and are held per category so that fillMatrix iterates over each category
 * on its own:
 *
 * - Domain faces: m_domainFluxStencils, with their boundary weights in m_boundaryWeights.
 * - EB faces: in the EB and jump boundary conditions, see below.
 * - Fine side of a refinement boundary: m_fineCFFluxStencils, multiplied by the fine-face b.
 * - Coarse side of a refinement boundary: m_coarCFFluxStencils, one stencil per fine face under the coarse face.
 *   Each is multiplied by the b-coefficient on its fine face, which fillMatrix first copies from the fine level into
 *   the index space of the coarse level (m_bCoefFiCo).
 *
 * So far the stencils of the domain faces and of both sides of the refinement boundaries are computed in fillMatrix
 * rather than stored; the holders, and the functions that define them, are placeholders.
 *
 * Embedded boundaries and interfaces follow this design:
 *
 * - Rows: every cut cell of each phase is a row. Multi-valued cells are a hard error (PetscGrid aborts on them).
 * - The flux through the EB face of a cut cell is F = w_B * g + sum_j w_j * phi_j, where the stencil is a PetscStencil
 *   over PETSc rows of any phase and any level. The shared primitive, in PetscHelmholtzEBStencils, is the one-sided
 *   normal derivative of phi over the cells of the cell's phase, dphi/dn = w_B * phi_B + sum_i w_i * phi_i, from a
 *   least-squares fit (LeastSquares::computeDualLevelStencils) over the cell's level and the next finer one.
 *   Coarse-fine ghost cells in it are replaced by their interpolation stencils, as for the domain faces, and cells
 *   covered by the finer level by the finer cells under them. It searches the half-circle, quadrant and radius
 *   neighbourhoods, each from the requested order down to 1, and aborts if none gives a stencil. One order and one
 *   weight, set by the solver, apply to all EB face stencils.
 * - Electrode faces: PetscHelmholtzEBBC, with Dirichlet (g = phi_B) and Neumann (g = dphi/dn) subclasses, and
 *   PetscHelmholtzElectrostaticEBBC (Source/Electrostatics), which takes g from the electrodes.
 * - Interface cells, where both phases have a VoF: PetscHelmholtzJumpBC. The one-sided derivatives of both phases are
 *   combined by a closure that eliminates phi_B, so that the flux into each side is a linear form over the rows of
 *   both phases plus w_g * g. The closure is the subclass: the natural jump b_1 D_1 + b_2 D_2 = sigma (g = sigma), and
 *   the saturation charge jump D_gas = 0 (Source/Electrostatics), for which sigma is an output. A virtual hook computes
 *   the jump data from the solution after the solve.
 * - The operator: kappa * alpha * a on the diagonal, each face flux times its area fraction, and the EB face flux
 *   times the boundary area and the b-coefficient on the EB. The flux through a cut face, and through the domain face
 *   of a cut cell, is the face-center fluxes of the neighbouring faces interpolated to the face centroid, as in
 *   EBHelmholtzOp. The b-coefficient on the EB is a third coefficient in define and setCoefficients.
 * - The one-sided least-squares stencils are computed once per define and stored in the boundary condition objects
 *   (PetscHelmholtzEBStencils). fillMatrix multiplies the electrode fluxes by beta and the b-coefficient on the EB, and
 *   fetches the interface fluxes, whose weights depend on the b-coefficients of both phases, with those coefficients.
 *   addBoundaryConditionsToRHS adds w_B * g for the electrode faces and w_g * g for the interface faces, and takes
 *   sigma as an argument.
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
   * @brief Release the matrix, the boundary conditions, the coefficients and all per-level data, so that the operator
   * can be defined again, e.g. after a regrid.
   */
  virtual void
  clear() noexcept;

  /**
   * @brief Define the operator.
   * @details This clears any previous definition, defines the domain, EB and jump BCs, creates the matrix, and fills
   * it. The a- and b-coefficients are held by reference, see setCoefficients.
   * @param[in] a_petscGrid     PETSc grid (rows)
   * @param[in] a_interpolators Coarse-fine ghost cell interpolators of each phase, on each level. Entry 0 is not used.
   * @param[in] a_domainBc      Domain boundary conditions
   * @param[in] a_ebBc          Boundary conditions on the EB faces of cut cells that are not interface cells
   * @param[in] a_jumpBc        Conditions on the EB faces of interface cells
   * @param[in] a_aCoef         Cell-centered a-coefficient
   * @param[in] a_bCoef         Face-centered b-coefficient
   * @param[in] a_bCoefIrreg    b-coefficient on the EB faces
   * @param[in] a_probLo        Lower-left corner of the computational domain
   * @param[in] a_dx            Grid resolution on each level
   * @param[in] a_alpha         Operator alpha
   * @param[in] a_beta          Operator beta
   */
  virtual void
  define(const RefCountedPtr<PetscGrid>&              a_petscGrid,
         const Vector<MFMultigridInterpolator>&       a_interpolators,
         const RefCountedPtr<PetscHelmholtzDomainBC>& a_domainBc,
         const RefCountedPtr<PetscHelmholtzEBBC>&     a_ebBc,
         const RefCountedPtr<PetscHelmholtzJumpBC>&   a_jumpBc,
         const MFAMRCellData&                         a_aCoef,
         const MFAMRFluxData&                         a_bCoef,
         const MFAMRIVData&                           a_bCoefIrreg,
         const RealVect&                              a_probLo,
         const Vector<Real>&                          a_dx,
         const Real                                   a_alpha,
         const Real                                   a_beta) noexcept;

  /**
   * @brief Set alpha and beta. The matrix is not changed until fillMatrix is called.
   * @details Before solving, the caller must call fillMatrix and rebuild the right-hand side, see the class
   * documentation.
   * @param[in] a_alpha Operator alpha
   * @param[in] a_beta  Operator beta
   */
  virtual void
  setAlphaAndBeta(const Real a_alpha, const Real a_beta) noexcept;

  /**
   * @brief Point the operator to other a- and b-coefficients. The matrix is not changed until fillMatrix is called.
   * @details The coefficients are held by reference, so later changes to their values are seen by fillMatrix. After
   * this call, or after changing the values in place, the caller must call fillMatrix and rebuild the right-hand side
   * before solving, see the class documentation. The coefficients must be defined on the grids of the PetscGrid;
   * fillMatrix aborts otherwise.
   * @param[in] a_aCoef      Cell-centered a-coefficient
   * @param[in] a_bCoef      Face-centered b-coefficient
   * @param[in] a_bCoefIrreg b-coefficient on the EB faces
   */
  virtual void
  setCoefficients(const MFAMRCellData& a_aCoef, const MFAMRFluxData& a_bCoef, const MFAMRIVData& a_bCoefIrreg) noexcept;

  /**
   * @brief Fill the matrix from the current alpha, beta, and a- and b-coefficients.
   * @details This refills the existing matrix, keeping its preallocation and sparsity pattern. It also stores the
   * boundary weights used by addBoundaryConditionsToRHS. Call it after every change of alpha, beta, a or b, and rebuild
   * the right-hand side afterwards. Aborts if the a- or b-coefficient is not defined on the grids of the PetscGrid.
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
   * @brief Get the weight that the rows of a level are multiplied by.
   * @details The rows are volume-weighted, (dx/dx_0)^SpaceDim times the kappa-weighted divergence, so that a flux
   * through a refinement boundary enters the coarse and fine rows with the same weight, and the row scale does not grow
   * by refRat^2 per level. The weight is 1 on level 0. The right-hand side must carry the same weights, see scaleRows.
   * @param[in] a_level Grid level
   * @return (dx[a_level]/dx[0])^SpaceDim
   */
  virtual Real
  getRowWeight(const int a_level) const noexcept;

  /**
   * @brief Multiply, or divide, each local row of a vector by its row weight.
   * @details Use it to weight a right-hand side built from cell data and addBoundaryConditionsToRHS before the solve,
   * and to remove the weight from a residual afterwards.
   * @param[in,out] a_vec     Vector over the PETSc rows
   * @param[in]     a_inverse If true, divide by the weights instead
   */
  virtual void
  scaleRows(Vec& a_vec, const bool a_inverse) const noexcept;

  /**
   * @brief Move the boundary data part of the domain boundary fluxes into a right-hand side.
   * @details The operator is A*phi + c, where A is the matrix and c holds the w_B * g parts of the boundary fluxes,
   * scaled like the rest of the face flux. Solving A*phi + c = rhs is the same as solving A*phi = rhs - c, so this
   * subtracts c from a_rhs. It evaluates the boundary data but not the stencils, using the w_B kept from assembly. Only
   * boundary terms are added: the domain boundary data, the electrode data, and the jump data (e.g. the surface charge
   * density). The caller supplies the volume source in a_rhs, kappa-weighted as described for the class. The scaling
   * uses the current beta and b, so call this after the fillMatrix that matches them.
   * @param[in,out] a_rhs   Right-hand side, on the rows of the PETSc grid
   * @param[in]     a_sigma Surface charge density on the interface between the phases, on the cut cells of the first
   * phase, in the units of the jump condition b_0 * D_0 + b_1 * D_1 = sigma
   */
  virtual void
  addBoundaryConditionsToRHS(Vec& a_rhs, const EBAMRIVData& a_sigma) const noexcept;

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
   * @brief Boundary conditions on the EB faces of cut cells that are not interface cells.
   */
  RefCountedPtr<PetscHelmholtzEBBC> m_ebBc;

  /**
   * @brief Conditions on the EB faces of interface cells.
   */
  RefCountedPtr<PetscHelmholtzJumpBC> m_jumpBc;

  /**
   * @brief Coarse-fine ghost cell interpolators of each phase, on each level. Entry 0 is not used.
   */
  Vector<MFMultigridInterpolator> m_interpolators;

  /**
   * @brief Cell-centered a-coefficient.
   */
  MFAMRCellData m_aCoef;

  /**
   * @brief Face-centered b-coefficient.
   */
  MFAMRFluxData m_bCoef;

  /**
   * @brief b-coefficient on the EB faces.
   */
  MFAMRIVData m_bCoefIrreg;

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
  EBAMRIFData m_boundaryWeights[2];

  /**
   * @brief Flux stencils through the domain faces, without beta and b, over PETSc rows.
   * @details Indexed as [dir][side][level], with Side::Lo = 0 and Side::Hi = 1. Each BaseFab is defined over the cells
   * of the box that abut the domain side, and holds the stencil of the flux through the domain face of each cell.
   * Placeholder: not defined or used yet.
   */
  Vector<RefCountedPtr<LayoutData<BaseFab<PetscStencil>>>> m_domainFluxStencils[SpaceDim][2];

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
   * @brief Abort unless the grids are supported by this operator: an interpolator for each phase on every level above
   * the coarsest, and no periodic directions.
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
   * @brief Allocate the boundary weights, and create the matrix with an exact preallocation.
   * @details The preallocation comes from assembling the operator once into a MATPREALLOCATOR, with fillMatrix. The
   * matrix entries are zero afterwards; fillMatrix must be called to fill them.
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
   * @brief Get the coarse-fine interpolator of a level and phase.
   * @param[in] a_level Grid level
   * @param[in] a_phase Phase
   * @return The interpolator, or null on level 0.
   */
  virtual RefCountedPtr<EBMultigridInterpolator>
  getInterpolator(const int a_level, const int a_phase) const noexcept;

  /**
   * @brief Get the cells around a grid that are covered by the finer level.
   * @details These are the cells that the interpolation of a cut-face flux to its centroid must not use. Empty on the
   * finest level and on grids without cut cells.
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_phase Phase
   * @return The covered cells within two cells of the grid.
   */
  virtual IntVectSet
  getCoveredCells(const int a_level, const DataIndex& a_din, const int a_phase) const noexcept;

  /**
   * @brief Check if a cell is an interface cell, i.e. cut in both phases.
   * @param[in] a_iv    Cell
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @return True if both phases are cut in a_iv.
   */
  virtual bool
  isInterfaceCell(const IntVect& a_iv, const int a_level, const DataIndex& a_din) const noexcept;

  /**
   * @brief Get the flux b * dphi/dx through the centroid of an interior cut face, as in EBHelmholtzOp.
   * @details The two-point fluxes through the face and its neighbouring face, with b on each, interpolated to the face
   * centroid. The stencil is over cells on the level, and may contain coarse-fine ghost cells.
   * @param[in] a_face         Face
   * @param[in] a_coveredCells Cells that the interpolation must not use, see getCoveredCells
   * @param[in] a_level        Grid level
   * @param[in] a_din          Grid index on the level
   * @param[in] a_phase        Phase
   * @return The flux stencil, including the b-coefficients but not the area fraction.
   */
  virtual VoFStencil
  getCutFaceFluxStencil(const FaceIndex&  a_face,
                        const IntVectSet& a_coveredCells,
                        const int         a_level,
                        const DataIndex&  a_din,
                        const int         a_phase) const noexcept;

  /**
   * @brief Get the row of the coarse cell across a refinement boundary from a fine cell.
   * @details Aborts if the coarse cell is not a row.
   * @param[in] a_ivGhost Coarse-fine ghost cell on level a_level, across the refinement boundary
   * @param[in] a_refRat  Refinement ratio between level a_level - 1 and level a_level
   * @param[in] a_level   Grid level of the fine cell
   * @param[in] a_din     Grid index on level a_level
   * @param[in] a_phase   Phase
   * @return The global PETSc row of the coarse cell under a_ivGhost.
   */
  virtual PetscInt
  getCoarseRow(const IntVect&   a_ivGhost,
               const int        a_refRat,
               const int        a_level,
               const DataIndex& a_din,
               const int        a_phase) const noexcept;

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
