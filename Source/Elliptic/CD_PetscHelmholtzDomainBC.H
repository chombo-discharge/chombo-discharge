/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzDomainBC.H
 * @brief  Declaration of a domain boundary condition class for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZDOMAINBC_H
#define CD_PETSCHELMHOLTZDOMAINBC_H

// Chombo includes
#include <EBLevelGrid.H>
#include <Stencils.H>
#include <RealVect.H>
#include <LoHiSide.H>
#include <RefCountedPtr.H>
#include <Vector.H>

// Our includes
#include <CD_PetscGrid.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Base class for domain boundary conditions of PetscHelmholtzOp.
 * @details This mirrors EBHelmholtzDomainBC, but where that class computes the flux through a domain face from data,
 * this class returns the flux as a linear form in the cell values and the boundary data g, i.e.
 *
 *   F = w_B * g + sum_i w_i * phi_i
 *
 * g is what the user specifies on the boundary, e.g. phi for a Dirichlet BC and dphi/dn for a Neumann BC. The weights
 * w_i and w_B depend only on the grid and the discretization, and come from getFaceFluxStencil: the w_i go into the
 * matrix, and w_B is kept by the operator. g comes from getBoundaryData, and w_B * g goes into the right-hand side, so
 * that g can change between solves without the matrix being reassembled. F is the derivative of phi in the +a_dir
 * direction on the face, i.e. it does not include the b-coefficient, which PetscHelmholtzOp multiplies in on the
 * face. The discretizations in the subclasses are the same as in the EBHelmholtz domain BCs, so that
 * PetscHelmholtzOp and EBHelmholtzOp discretize the same problem.
 *
 * Like the PETSc matrix, the BC is composite: it is defined once over the whole AMR hierarchy, and each query names the
 * phase, level, grid index and cell it is for. Through the PetscGrid a subclass has the grids on every level, their
 * coarsenings and refinements (which share DataIndex with the level they came from), and knows for every cell whether
 * it is covered by a finer grid. This is what a stencil needs in order to reach into a coarser or finer level.
 */
class PetscHelmholtzDomainBC
{
public:
  /**
   * @brief Default constructor. Must subsequently call define.
   */
  PetscHelmholtzDomainBC();

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzDomainBC();

  /**
   * @brief Disallowed copy constructor.
   * @param[in] a_other Other BC
   */
  PetscHelmholtzDomainBC(const PetscHelmholtzDomainBC& a_other) = delete;

  /**
   * @brief Disallowed move constructor.
   * @param[in] a_other Other BC
   */
  PetscHelmholtzDomainBC(const PetscHelmholtzDomainBC&& a_other) = delete;

  /**
   * @brief Disallowed copy assignment operator.
   * @param[in] a_other Other BC
   * @return This BC.
   */
  PetscHelmholtzDomainBC&
  operator=(const PetscHelmholtzDomainBC& a_other) = delete;

  /**
   * @brief Disallowed move assignment operator.
   * @param[in] a_other Other BC
   * @return This BC.
   */
  PetscHelmholtzDomainBC&
  operator=(const PetscHelmholtzDomainBC&& a_other) = delete;

  /**
   * @brief Define the BC over the AMR hierarchy.
   * @param[in] a_petscGrid PETSc grid. Holds the grids on every level, their coarsenings and refinements, and the
   * per-cell AMR information (e.g. whether a cell is covered by a finer grid).
   * @param[in] a_probLo    Lower-left corner of the computational domain
   * @param[in] a_dx        Grid resolution on each level
   */
  virtual void
  define(const RefCountedPtr<PetscGrid>& a_petscGrid, const RealVect& a_probLo, const Vector<Real>& a_dx);

  /**
   * @brief Get the stencil of the face flux through a domain face, and the weight of the boundary data.
   * @details The stencil gives the weights w_i of the flux, excluding the b-coefficient. It is returned
   * in three parts, one per level, which are summed: a_stencil over cells on level a_level in the index space of grid
   * a_din, a_stencilCoar over cells on level a_level - 1 in the index space of grid a_din coarsened
   * (PetscGrid::getMFLevelGridsCoFi()[a_level - 1]), and a_stencilFine over cells on level a_level + 1 in the index
   * space of grid a_din refined (PetscGrid::getMFLevelGridsFiCo()[a_level + 1]). Previous contents of the three
   * stencils are discarded.
   *
   * This function does not require the entries to be PETSc rows: they may be e.g. coarse-fine ghost cells or cells
   * covered by a finer grid. It is up to the operator that uses the stencil to make it consistent with the PetscGrid
   * before it goes into the matrix, by replacing the entries it knows how to replace and failing on the rest. See
   * PetscHelmholtzOp for what that operator accepts.
   * @param[out] a_boundaryWeight Weight w_B of the boundary data in the flux, see getBoundaryData
   * @param[out] a_stencilCoar Stencil over cells on level a_level - 1
   * @param[out] a_stencil     Stencil over cells on level a_level
   * @param[out] a_stencilFine Stencil over cells on level a_level + 1
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cell adjacent to the domain face
   * @param[in] a_dir   Coordinate direction of the face
   * @param[in] a_side  Domain side of the face
   */
  virtual void
  getFaceFluxStencil(Real&                a_boundaryWeight,
                     VoFStencil&          a_stencilCoar,
                     VoFStencil&          a_stencil,
                     VoFStencil&          a_stencilFine,
                     const int            a_phase,
                     const int            a_level,
                     const DataIndex&     a_din,
                     const VolIndex&      a_vof,
                     const int            a_dir,
                     const Side::LoHiSide a_side) const = 0;

  /**
   * @brief Get the boundary data g at a domain face.
   * @details The flux depends on g through w_B * g, with w_B from getFaceFluxStencil. g is evaluated whenever the
   * right-hand side is built.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cell adjacent to the domain face
   * @param[in] a_dir   Coordinate direction of the face
   * @param[in] a_side  Domain side of the face
   * @return The boundary data g.
   */
  virtual Real
  getBoundaryData(const int            a_phase,
                  const int            a_level,
                  const DataIndex&     a_din,
                  const VolIndex&      a_vof,
                  const int            a_dir,
                  const Side::LoHiSide a_side) const = 0;

protected:
  /**
   * @brief Is defined or not.
   */
  bool m_isDefined;

  /**
   * @brief PETSc grid.
   * @details For grids on level lvl, m_petscGrid->getMFLevelGridsCoFi()[lvl - 1] holds them coarsened to level lvl - 1,
   * and m_petscGrid->getMFLevelGridsFiCo()[lvl + 1] holds them refined to level lvl + 1. Both share DataIndex with
   * level lvl. m_petscGrid->getAMRToPetsc()[lvl] has the per-cell AMR information.
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
   * @brief Get the physical position of the centroid of a domain face.
   * @details The face is the face of a_vof on the domain side, in the EB graph of the phase. On a face that is not cut
   * by the EB, the centroid is the face center. A cell with more than one face on the domain side is not supported.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cell adjacent to the domain face
   * @param[in] a_dir   Coordinate direction of the face
   * @param[in] a_side  Domain side of the face
   * @return Position of the face centroid.
   */
  RealVect
  getBoundaryPosition(const int            a_phase,
                      const int            a_level,
                      const DataIndex&     a_din,
                      const VolIndex&      a_vof,
                      const int            a_dir,
                      const Side::LoHiSide a_side) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
