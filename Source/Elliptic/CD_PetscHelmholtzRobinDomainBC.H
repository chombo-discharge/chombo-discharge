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
 * neighbor, and the flux is C/B + s*A*phi_face/B with s = -1 on the low side and +1 on the high side. The interior
 * neighbor must be a cell on the same level, i.e. not covered by the EB or by a finer grid. Where it is not, the
 * extrapolation needs a least-squares reconstruction over the valid cells on this level and the next finer level,
 * which is not implemented yet and aborts.
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
   * @brief Get the boundary data, i.e. the coefficient C.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cell adjacent to the domain face
   * @param[in] a_dir   Coordinate direction of the face
   * @param[in] a_side  Domain side of the face
   * @return C at the face.
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
   * @param[in]  a_phase Phase
   * @param[in]  a_level Grid level
   * @param[in]  a_din   Grid index on the level
   * @param[in]  a_vof   Cell adjacent to the domain face
   * @param[in]  a_dir   Coordinate direction of the face
   * @param[in]  a_side  Domain side of the face
   */
  void
  getCoefficients(Real&                a_A,
                  Real&                a_B,
                  Real&                a_C,
                  const int            a_phase,
                  const int            a_level,
                  const DataIndex&     a_din,
                  const VolIndex&      a_vof,
                  const int            a_dir,
                  const Side::LoHiSide a_side) const;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
