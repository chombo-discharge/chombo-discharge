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
#include <CD_NamespaceHeader.H>

/**
 * @brief Shared one-sided EB face stencils for the EB and jump boundary conditions of PetscHelmholtzOp.
 * @details For each cut cell of each phase, this class holds the one-sided normal derivative of phi at the EB face,
 *
 *   dphi/dn = w_B * phi_B + sum_i w_i * phi_i,
 *
 * where phi_B is the value on the EB face, n points out of the cell (into the boundary), and the stencil is over PETSc
 * rows of the cell's phase. The stencil uses valid cells only: no coarse-fine ghost cells and no cells covered by a
 * finer grid. It takes cells from the cell's level and at most one neighbouring level (the coarser level on the fine
 * side of a refinement boundary, the finer level next to covered cells). The search tries the half-circle, quadrant and
 * radius neighbourhoods, each from the order down to 1, and aborts if none gives a stencil. The stencils are computed
 * once, in define, and do not contain the b-coefficient or the face area.
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
   * @brief Define over the AMR hierarchy, and compute the one-sided normal derivative stencils of every cut cell.
   * @details The order and the weight must be set first.
   * @param[in] a_petscGrid PETSc grid. Holds the grids on every level, their coarsenings and refinements, and the
   * per-cell AMR information.
   * @param[in] a_probLo    Lower-left corner of the computational domain
   * @param[in] a_dx        Grid resolution on each level
   */
  virtual void
  define(const RefCountedPtr<PetscGrid>& a_petscGrid, const RealVect& a_probLo, const Vector<Real>& a_dx) noexcept;

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
   * @brief Compute the one-sided normal derivative at the EB face of a cut cell, from valid cells only.
   * @details Runs the neighbourhood and order search described for the class, and aborts if no stencil is found.
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
