/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzDirichletDomainBC.H
 * @brief  Declaration of a Dirichlet domain boundary condition class for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZDIRICHLETDOMAINBC_H
#define CD_PETSCHELMHOLTZDIRICHLETDOMAINBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Dirichlet domain boundary condition for PetscHelmholtzOp.
 * @details The face derivative is the one-sided difference between the cell center and the face,
 * (phi_B - phi)/(dx/2) on the high side and (phi - phi_B)/(dx/2) on the low side, as in EBHelmholtzDirichletDomainBC.
 */
class PetscHelmholtzDirichletDomainBC : public PetscHelmholtzDomainBC
{
public:
  /**
   * @brief Default constructor. Must subsequently call setValue.
   */
  PetscHelmholtzDirichletDomainBC();

  /**
   * @brief Constructor with a constant value.
   * @param[in] a_value Value of phi on the boundary
   */
  PetscHelmholtzDirichletDomainBC(const Real a_value);

  /**
   * @brief Constructor with a spatially varying value.
   * @param[in] a_value Value of phi on the boundary
   */
  PetscHelmholtzDirichletDomainBC(const std::function<Real(const RealVect& a_pos)>& a_value);

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzDirichletDomainBC();

  /**
   * @brief Set a constant value on the boundary.
   * @param[in] a_value Value of phi on the boundary
   */
  void
  setValue(const Real a_value);

  /**
   * @brief Set a spatially varying value on the boundary.
   * @param[in] a_value Value of phi on the boundary
   */
  void
  setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /**
   * @brief Get the stencil of the face flux through a domain face.
   * @param[in] a_vof  Cell adjacent to the domain face
   * @param[in] a_dir  Coordinate direction of the face
   * @param[in] a_side Domain side of the face
   * @param[in] a_level Level of the cell
   * @param[in] a_din  Grid index of the cell
   * @return The stencil as (cells on this level, cells on the next finer level).
   */
  std::pair<VoFStencil, VoFStencil>
  getFaceFluxStencil(const VolIndex&      a_vof,
                     const int            a_dir,
                     const Side::LoHiSide a_side,
                     const int            a_level,
                     const DataIndex&     a_din) const override;

  /**
   * @brief Get the constant part of the face flux through a domain face.
   * @param[in] a_vof  Cell adjacent to the domain face
   * @param[in] a_dir  Coordinate direction of the face
   * @param[in] a_side Domain side of the face
   * @param[in] a_level Level of the cell
   * @param[in] a_din  Grid index of the cell
   * @return The constant part of the flux.
   */
  Real
  getFaceFluxConstant(const VolIndex&      a_vof,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const int            a_level,
                      const DataIndex&     a_din) const override;

protected:
  /**
   * @brief Use a constant value or not.
   */
  bool m_useConstant;

  /**
   * @brief Use a spatially varying value or not.
   */
  bool m_useFunction;

  /**
   * @brief Constant value.
   */
  Real m_constantValue;

  /**
   * @brief Spatially varying value.
   */
  std::function<Real(const RealVect& a_pos)> m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
