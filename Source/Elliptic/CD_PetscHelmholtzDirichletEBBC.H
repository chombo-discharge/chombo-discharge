/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzDirichletEBBC.H
 * @brief  Declaration of a Dirichlet boundary condition on the EB faces for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZDIRICHLETEBBC_H
#define CD_PETSCHELMHOLTZDIRICHLETEBBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_PetscHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Dirichlet boundary condition phi = phi_B on the EB faces, for PetscHelmholtzOp.
 * @details The flux is the one-sided normal derivative from PetscHelmholtzEBStencils, F = w_B * phi_B + sum_i w_i *
 * phi_i, so g = phi_B.
 */
class PetscHelmholtzDirichletEBBC : public PetscHelmholtzEBBC
{
public:
  /**
   * @brief Default constructor. Must subsequently set the value, the order and the weight, and call define.
   */
  PetscHelmholtzDirichletEBBC() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzDirichletEBBC() noexcept;

  /**
   * @brief Set a constant value on the EB.
   * @param[in] a_value Value of phi on the EB
   */
  virtual void
  setValue(const Real a_value) noexcept;

  /**
   * @brief Set a spatially varying value on the EB.
   * @param[in] a_value Value of phi on the EB, as a function of position
   */
  virtual void
  setValue(const std::function<Real(const RealVect& a_pos)>& a_value) noexcept;

  /**
   * @brief Get the flux through the EB face of a cut cell, without the b-coefficient and the face area.
   * @param[out] a_boundaryWeight Weight w_B of the boundary data
   * @param[out] a_stencil        Stencil over PETSc rows
   * @param[in]  a_phase          Phase
   * @param[in]  a_level          Grid level
   * @param[in]  a_din            Grid index on the level
   * @param[in]  a_vof            Cut cell
   */
  virtual void
  getEBFluxStencil(Real&            a_boundaryWeight,
                   PetscStencil&    a_stencil,
                   const int        a_phase,
                   const int        a_level,
                   const DataIndex& a_din,
                   const VolIndex&  a_vof) const noexcept override;

  /**
   * @brief Get the boundary data, i.e. phi on the EB.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cut cell
   * @return phi on the EB.
   */
  virtual Real
  getBoundaryData(const int        a_phase,
                  const int        a_level,
                  const DataIndex& a_din,
                  const VolIndex&  a_vof) const noexcept override;

protected:
  /**
   * @brief Use the constant value or not.
   */
  bool m_useConstant;

  /**
   * @brief Use the spatially varying value or not.
   */
  bool m_useFunction;

  /**
   * @brief Constant value on the EB.
   */
  Real m_constantValue;

  /**
   * @brief Spatially varying value on the EB.
   */
  std::function<Real(const RealVect& a_pos)> m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
