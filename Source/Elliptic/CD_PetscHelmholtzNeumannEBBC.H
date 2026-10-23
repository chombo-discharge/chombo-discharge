/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzNeumannEBBC.H
 * @brief  Declaration of a Neumann boundary condition on the EB faces for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZNEUMANNEBBC_H
#define CD_PETSCHELMHOLTZNEUMANNEBBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_PetscHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Neumann boundary condition dphi/dn = g on the EB faces, for PetscHelmholtzOp.
 * @details n points out of the cell, into the EB, as in EBHelmholtzNeumannEBBC. The flux is the given normal
 * derivative, so the stencil is empty and w_B = 1.
 */
class PetscHelmholtzNeumannEBBC : public PetscHelmholtzEBBC
{
public:
  /**
   * @brief Default constructor. Must subsequently set the normal derivative, and call define.
   */
  PetscHelmholtzNeumannEBBC() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzNeumannEBBC() noexcept;

  /**
   * @brief Set a constant normal derivative on the EB.
   * @param[in] a_DphiDn Normal derivative of phi on the EB
   */
  virtual void
  setDphiDn(const Real a_DphiDn) noexcept;

  /**
   * @brief Set a spatially varying normal derivative on the EB.
   * @param[in] a_DphiDn Normal derivative of phi on the EB, as a function of position
   */
  virtual void
  setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn) noexcept;

  /**
   * @brief Get the flux through the EB face of a cut cell. The stencil is empty and w_B = 1.
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
   * @brief Get the boundary data, i.e. dphi/dn on the EB.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cut cell
   * @return dphi/dn on the EB.
   */
  virtual Real
  getBoundaryData(const int        a_phase,
                  const int        a_level,
                  const DataIndex& a_din,
                  const VolIndex&  a_vof) const noexcept override;

protected:
  /**
   * @brief Use the constant normal derivative or not.
   */
  bool m_useConstant;

  /**
   * @brief Use the spatially varying normal derivative or not.
   */
  bool m_useFunction;

  /**
   * @brief Constant normal derivative on the EB.
   */
  Real m_constantDphiDn;

  /**
   * @brief Spatially varying normal derivative on the EB.
   */
  std::function<Real(const RealVect& a_pos)> m_functionDphiDn;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
