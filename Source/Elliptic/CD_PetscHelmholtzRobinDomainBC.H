/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzRobinDomainBC.H
 * @brief  Declaration of a Robin domain boundary condition class for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZROBINDOMAINBC_H
#define CD_PETSCHELMHOLTZROBINDOMAINBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Robin domain boundary condition A*phi + B*dphi/dn = C for PetscHelmholtzOp.
 * @details As in EBHelmholtzRobinDomainBC, phi is extrapolated linearly to the face from the cell and its interior
 * neighbor, and the flux is C/B + s*A*phi_face/B with s = -1 on the low side and +1 on the high side.
 */
class PetscHelmholtzRobinDomainBC : public PetscHelmholtzDomainBC
{
public:
  /**
   * @brief Default constructor. Must subsequently call setCoefficients.
   */
  PetscHelmholtzRobinDomainBC();

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzRobinDomainBC();

  /**
   * @brief Set constant coefficients.
   * @param[in] a_A Coefficient on phi
   * @param[in] a_B Coefficient on dphi/dn
   * @param[in] a_C Right-hand side
   */
  void
  setCoefficients(const Real a_A, const Real a_B, const Real a_C);

  /**
   * @brief Set spatially varying coefficients.
   * @param[in] a_A Coefficient on phi
   * @param[in] a_B Coefficient on dphi/dn
   * @param[in] a_C Right-hand side
   */
  void
  setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                  const std::function<Real(const RealVect& a_pos)>& a_B,
                  const std::function<Real(const RealVect& a_pos)>& a_C);

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
   * @return The constant part of the flux, C/B.
   */
  Real
  getFaceFluxConstant(const VolIndex&      a_vof,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const int            a_level,
                      const DataIndex&     a_din) const override;

protected:
  /**
   * @brief Use constant coefficients or not.
   */
  bool m_useConstant;

  /**
   * @brief Use spatially varying coefficients or not.
   */
  bool m_useFunction;

  /**
   * @brief Constant coefficient on phi.
   */
  Real m_constantA;

  /**
   * @brief Constant coefficient on dphi/dn.
   */
  Real m_constantB;

  /**
   * @brief Constant right-hand side.
   */
  Real m_constantC;

  /**
   * @brief Spatially varying coefficient on phi.
   */
  std::function<Real(const RealVect& a_pos)> m_functionA;

  /**
   * @brief Spatially varying coefficient on dphi/dn.
   */
  std::function<Real(const RealVect& a_pos)> m_functionB;

  /**
   * @brief Spatially varying right-hand side.
   */
  std::function<Real(const RealVect& a_pos)> m_functionC;

  /**
   * @brief Get the coefficients at a domain face.
   * @param[out] a_A    Coefficient on phi
   * @param[out] a_B    Coefficient on dphi/dn
   * @param[out] a_C    Right-hand side
   * @param[in]  a_iv   Cell adjacent to the domain face
   * @param[in]  a_dir  Coordinate direction of the face
   * @param[in]  a_side Domain side of the face
   * @param[in]  a_level Level of the cell
   */
  void
  getCoefficients(Real&                a_A,
                  Real&                a_B,
                  Real&                a_C,
                  const IntVect&       a_iv,
                  const int            a_dir,
                  const Side::LoHiSide a_side,
                  const int            a_level) const;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
