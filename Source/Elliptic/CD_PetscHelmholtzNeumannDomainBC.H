/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzNeumannDomainBC.H
 * @brief  Declaration of a Neumann domain boundary condition class for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZNEUMANNDOMAINBC_H
#define CD_PETSCHELMHOLTZNEUMANNDOMAINBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Neumann domain boundary condition for PetscHelmholtzOp.
 * @details As in EBHelmholtzNeumannDomainBC, a positive dphi/dn is a flux into the domain. The flux has no
 * phi-dependent part, so the stencil is empty.
 */
class PetscHelmholtzNeumannDomainBC : public PetscHelmholtzDomainBC
{
public:
  /**
   * @brief Default constructor. Must subsequently call setDphiDn.
   */
  PetscHelmholtzNeumannDomainBC();

  /**
   * @brief Constructor with a constant normal derivative.
   * @param[in] a_DphiDn Value of dphi/dn on the boundary
   */
  PetscHelmholtzNeumannDomainBC(const Real a_DphiDn);

  /**
   * @brief Constructor with a spatially varying normal derivative.
   * @param[in] a_DphiDn Value of dphi/dn on the boundary
   */
  PetscHelmholtzNeumannDomainBC(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzNeumannDomainBC();

  /**
   * @brief Set a constant normal derivative on the boundary.
   * @param[in] a_DphiDn Value of dphi/dn on the boundary
   */
  void
  setDphiDn(const Real a_DphiDn);

  /**
   * @brief Set a spatially varying normal derivative on the boundary.
   * @param[in] a_DphiDn Value of dphi/dn on the boundary
   */
  void
  setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /**
   * @brief Get the stencil of the face flux through a domain face. The stencil is empty and w_B = -s, with s = -1 on
   * the low side and +1 on the high side: a positive dphi/dn is a flux into the domain, as in
   * EBHelmholtzNeumannDomainBC.
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
  void
  getFaceFluxStencil(Real&                a_boundaryWeight,
                     VoFStencil&          a_stencilCoar,
                     VoFStencil&          a_stencil,
                     VoFStencil&          a_stencilFine,
                     const int            a_phase,
                     const int            a_level,
                     const DataIndex&     a_din,
                     const VolIndex&      a_vof,
                     const int            a_dir,
                     const Side::LoHiSide a_side) const override;

  /**
   * @brief Get the boundary data, i.e. dphi/dn on the boundary.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cell adjacent to the domain face
   * @param[in] a_dir   Coordinate direction of the face
   * @param[in] a_side  Domain side of the face
   * @return dphi/dn on the boundary.
   */
  Real
  getBoundaryData(const int            a_phase,
                  const int            a_level,
                  const DataIndex&     a_din,
                  const VolIndex&      a_vof,
                  const int            a_dir,
                  const Side::LoHiSide a_side) const override;

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
   * @brief Constant dphi/dn.
   */
  Real m_constantDphiDn;

  /**
   * @brief Spatially varying dphi/dn.
   */
  std::function<Real(const RealVect& a_pos)> m_functionDphiDn;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
