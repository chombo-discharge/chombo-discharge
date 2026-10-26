/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzEBStencils.H
 * @brief  Declaration of the shared EB face stencils for the EB and jump boundary conditions of PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZEBSTENCILS_H
#define CD_PETSCHELMHOLTZEBSTENCILS_H

// Chombo includes
#include <RefCountedPtr.H>
#include <RealVect.H>
#include <Vector.H>
#include <LayoutData.H>
#include <BaseIVFAB.H>

// Our includes
#include <CD_PetscGrid.H>
#include <CD_PetscStencil.H>
#include <CD_MFMultigridInterpolator.H>
#include <CD_VofUtils.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Shared one-sided EB face stencils for the EB and jump boundary conditions of PetscHelmholtzOp.
 * @details For each cut cell of each phase, this class holds the one-sided normal derivative of phi at the EB face,
 *
 *   dphi/dn = w_B * phi_B + sum_i w_i * phi_i,
 *
 * where phi_B is the value on the EB face, n points out of the cell (into the boundary, i.e. minus the EB normal), and
 * the stencil is over PETSc rows of the cell's phase. phi is cell-centered.
 *
 * The stencil is a weighted least-squares fit at the EB face centroid, with phi_B known, over the cells in a
 * neighbourhood of the cut cell on its own level, within a radius equal to the order:
 *
 * - Cells that are PETSc rows are used as they are.
 * - Coarse-fine ghost cells are used as well, and are replaced by their interpolation stencils from the
 *   EBMultigridInterpolator (PetscGrid::getPetscStencilWithGhostCells), which are over rows on this level and the next
 *   coarser one. The interpolator fills the ghost cells within its radius of a cut cell, so the order must not exceed
 *   that radius.
 * - Cells covered by the next finer level are replaced by the finer cells under them that are PETSc rows.
 *
 * Every cell of the neighbourhood enters the fit, as in FieldSolverGMG. The search tries the half-circle, quadrant and
 * radius neighbourhoods, each from the order down to 1. A fit is accepted only if w_B > 0, i.e. if the normal
 * derivative increases with the boundary value, as in FieldSolverGMG. If no fit is accepted, or if a cut cell has no EB
 * normal, the cell has no flux through its EB face: its stencil is empty with w_B = 0. The stencils are computed once,
 * in define, and contain neither the b-coefficient nor the face area.
 *
 * PetscHelmholtzEBBC (electrode faces) and PetscHelmholtzJumpBC (interface faces) derive from this class and build
 * their fluxes from these stencils.
 */
class PetscHelmholtzEBStencils
{
public:
  /**
   * @brief Default constructor. Must subsequently set the order and the weight, and call define.
   */
  PetscHelmholtzEBStencils() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzEBStencils() noexcept;

  /**
   * @brief Disallowed copy constructor.
   * @param[in] a_other Other object
   */
  PetscHelmholtzEBStencils(const PetscHelmholtzEBStencils& a_other) = delete;

  /**
   * @brief Disallowed move constructor.
   * @param[in] a_other Other object
   */
  PetscHelmholtzEBStencils(const PetscHelmholtzEBStencils&& a_other) = delete;

  /**
   * @brief Disallowed copy assignment operator.
   * @param[in] a_other Other object
   * @return This object.
   */
  PetscHelmholtzEBStencils&
  operator=(const PetscHelmholtzEBStencils& a_other) = delete;

  /**
   * @brief Disallowed move assignment operator.
   * @param[in] a_other Other object
   * @return This object.
   */
  PetscHelmholtzEBStencils&
  operator=(const PetscHelmholtzEBStencils&& a_other) = delete;

  /**
   * @brief Set the order of the stencils.
   * @param[in] a_order Highest order of the stencils. The search falls back to lower orders, down to 1.
   */
  virtual void
  setOrder(const int a_order) noexcept;

  /**
   * @brief Set the least-squares weighting.
   * @param[in] a_weight Exponent of the inverse-distance weighting in the least-squares fit.
   */
  virtual void
  setWeight(const int a_weight) noexcept;

  /**
   * @brief Release the stored stencils, the grids and the interpolators, so that the object can be defined again, e.g.
   * after a regrid.
   * @details The order and the weight are settings, and are kept.
   */
  virtual void
  clear() noexcept;

  /**
   * @brief Define over the AMR hierarchy, and compute the one-sided normal derivative stencils of every cut cell.
   * @details The order and the weight must be set first. Clears any previous definition first. Aborts if the order
   * exceeds the radius within which the interpolators fill the ghost cells.
   * @param[in] a_petscGrid     PETSc grid. Holds the grids on every level, their coarsenings and refinements, and the
   * per-cell AMR information.
   * @param[in] a_interpolators Coarse-fine ghost cell interpolators of both phases, on each level. Entry 0 is not used.
   * @param[in] a_probLo        Lower-left corner of the computational domain
   * @param[in] a_dx            Grid resolution on each level
   */
  virtual void
  define(const RefCountedPtr<PetscGrid>&        a_petscGrid,
         const Vector<MFMultigridInterpolator>& a_interpolators,
         const RealVect&                        a_probLo,
         const Vector<Real>&                    a_dx) noexcept;

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
   * @brief Highest order of the stencils, or -1 if not set.
   */
  int m_order;

  /**
   * @brief Least-squares weighting, or -1 if not set.
   */
  int m_weight;

  /**
   * @brief PETSc grid.
   */
  RefCountedPtr<PetscGrid> m_petscGrid;

  /**
   * @brief Coarse-fine ghost cell interpolators of both phases, on each level.
   */
  Vector<MFMultigridInterpolator> m_interpolators;

  /**
   * @brief Lower-left corner of the computational domain.
   */
  RealVect m_probLo;

  /**
   * @brief Grid resolution on each level.
   */
  Vector<Real> m_dx;

  /**
   * @brief The stencils of the one-sided normal derivatives, over PETSc rows.
   * @details Indexed as [phase][level], with one stencil per cut cell of the phase.
   */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<PetscStencil>>>> m_derivStencils[2];

  /**
   * @brief The weights w_B of the boundary value in the one-sided normal derivatives.
   * @details Indexed as [phase][level], with one weight per cut cell of the phase.
   */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<Real>>>> m_derivWeights[2];

  /**
   * @brief Get the stored one-sided normal derivative at the EB face of a cut cell.
   * @param[out] a_boundaryWeight Weight w_B of the boundary value phi_B
   * @param[out] a_stencil        Stencil over PETSc rows of the cell's phase
   * @param[in]  a_phase          Phase
   * @param[in]  a_level          Grid level
   * @param[in]  a_din            Grid index on the level
   * @param[in]  a_vof            Cut cell
   */
  virtual void
  getNormalDerivativeStencil(Real&            a_boundaryWeight,
                             PetscStencil&    a_stencil,
                             const int        a_phase,
                             const int        a_level,
                             const DataIndex& a_din,
                             const VolIndex&  a_vof) const noexcept;

  /**
   * @brief Compute the one-sided normal derivative at the EB face of a cut cell.
   * @details Runs the neighbourhood and order search described for the class. If no stencil is found, w_B = 0 and the
   * stencil is empty.
   * @param[out] a_boundaryWeight Weight w_B of the boundary value phi_B
   * @param[out] a_stencil        Stencil over PETSc rows of the cell's phase
   * @param[in]  a_phase          Phase
   * @param[in]  a_level          Grid level
   * @param[in]  a_din            Grid index on the level
   * @param[in]  a_vof            Cut cell
   */
  virtual void
  computeNormalDerivativeStencil(Real&            a_boundaryWeight,
                                 PetscStencil&    a_stencil,
                                 const int        a_phase,
                                 const int        a_level,
                                 const DataIndex& a_din,
                                 const VolIndex&  a_vof) const noexcept;

  /**
   * @brief Try to compute the one-sided normal derivative at the EB face of a cut cell, for one neighbourhood and one
   * order.
   * @details Gathers the cells of the neighbourhood as described for the class, and fits them by least squares.
   * @param[out] a_boundaryWeight Weight w_B of the boundary value phi_B
   * @param[out] a_stencil        Stencil over PETSc rows of the cell's phase
   * @param[in]  a_phase          Phase
   * @param[in]  a_level          Grid level
   * @param[in]  a_din            Grid index on the level
   * @param[in]  a_vof            Cut cell
   * @param[in]  a_neighborhood   Neighbourhood of the cut cell to take cells from
   * @param[in]  a_order          Order of the fit. Also the radius of the neighbourhood.
   * @return True if there were enough cells for the fit, it succeeded and w_B > 0, in which case a_boundaryWeight and
   * a_stencil are set.
   */
  virtual bool
  getLeastSquaresStencil(Real&                        a_boundaryWeight,
                         PetscStencil&                a_stencil,
                         const int                    a_phase,
                         const int                    a_level,
                         const DataIndex&             a_din,
                         const VolIndex&              a_vof,
                         const VofUtils::Neighborhood a_neighborhood,
                         const int                    a_order) const noexcept;

  /**
   * @brief Get the physical position of the centroid of the EB face of a cut cell.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cut cell
   * @return Position of the EB face centroid.
   */
  virtual RealVect
  getBoundaryPosition(const int        a_phase,
                      const int        a_level,
                      const DataIndex& a_din,
                      const VolIndex&  a_vof) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
