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
   * @brief Get the boundary data, i.e. phi on the boundary.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cell adjacent to the domain face
   * @param[in] a_dir   Coordinate direction of the face
   * @param[in] a_side  Domain side of the face
   * @return phi on the boundary.
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
