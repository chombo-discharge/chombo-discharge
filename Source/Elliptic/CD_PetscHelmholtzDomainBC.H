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

// Std includes
#include <utility>

// Chombo includes
#include <EBLevelGrid.H>
#include <Stencils.H>
#include <RealVect.H>
#include <LoHiSide.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Base class for domain boundary conditions of PetscHelmholtzOp.
 * @details This mirrors EBHelmholtzDomainBC, but where that class computes the flux through a domain face from data,
 * this class returns the flux as a linear form in the cell values, i.e.
 *
 *   F = sum_i w_i * phi_i + c
 *
 * The weights w_i go into the matrix and the constant c into the right-hand side. F is the derivative of phi in the
 * +a_dir direction on the face, i.e. it does not include the b-coefficient, which PetscHelmholtzOp multiplies in on
 * the face. The discretizations in the subclasses are the same as in the EBHelmholtz domain BCs, so that
 * PetscHelmholtzOp and EBHelmholtzOp discretize the same problem.
 *
 * The stencil is returned as two parts, (level, level+1), the form taken by PetscGrid::getPetscStencilWithFinerLevel.
 * The level+1 part is reserved for boundary stencils that reach under a finer level; it is empty for now.
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
   * @brief Define the BC on a level.
   * @param[in] a_dataLocation Data location. Only cell centers are supported.
   * @param[in] a_eblg         Grids on the level
   * @param[in] a_probLo       Lower-left corner of the computational domain
   * @param[in] a_dx           Grid resolution on the level
   */
  virtual void
  define(const Location::Cell a_dataLocation, const EBLevelGrid& a_eblg, const RealVect& a_probLo, const Real a_dx);

  /**
   * @brief Get the stencil of the face flux through a domain face.
   * @details The stencil gives the part of the flux that depends on phi, excluding the b-coefficient.
   * @param[in] a_vof  Cell adjacent to the domain face
   * @param[in] a_dir  Coordinate direction of the face
   * @param[in] a_side Domain side of the face
   * @param[in] a_din  Grid index of the cell
   * @return The stencil as (cells on this level, cells on the next finer level).
   */
  virtual std::pair<VoFStencil, VoFStencil>
  getFaceFluxStencil(const VolIndex&      a_vof,
                     const int            a_dir,
                     const Side::LoHiSide a_side,
                     const DataIndex&     a_din) const = 0;

  /**
   * @brief Get the constant (inhomogeneous) part of the face flux through a domain face.
   * @details Excludes the b-coefficient, like getFaceFluxStencil.
   * @param[in] a_vof  Cell adjacent to the domain face
   * @param[in] a_dir  Coordinate direction of the face
   * @param[in] a_side Domain side of the face
   * @param[in] a_din  Grid index of the cell
   * @return The constant part of the flux.
   */
  virtual Real
  getFaceFluxConstant(const VolIndex&      a_vof,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const DataIndex&     a_din) const = 0;

protected:
  /**
   * @brief Is defined or not.
   */
  bool m_isDefined;

  /**
   * @brief Data location.
   */
  Location::Cell m_dataLocation;

  /**
   * @brief Grids on the level.
   */
  EBLevelGrid m_eblg;

  /**
   * @brief Lower-left corner of the computational domain.
   */
  RealVect m_probLo;

  /**
   * @brief Grid resolution on the level.
   */
  Real m_dx;

  /**
   * @brief Get the physical position of the center of a domain face.
   * @param[in] a_iv   Cell adjacent to the domain face
   * @param[in] a_dir  Coordinate direction of the face
   * @param[in] a_side Domain side of the face
   * @return Position of the face center.
   */
  RealVect
  getBoundaryPosition(const IntVect& a_iv, const int a_dir, const Side::LoHiSide a_side) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
