/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzElectrostaticDomainBC.H
 * @brief  Declaration of an electrostatic domain boundary condition class for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZELECTROSTATICDOMAINBC_H
#define CD_PETSCHELMHOLTZELECTROSTATICDOMAINBC_H

// Std includes
#include <map>
#include <memory>

// Our includes
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_ElectrostaticDomainBc.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Domain boundary conditions for PetscHelmholtzOp, taken from ElectrostaticDomainBc.
 * @details This mirrors EBHelmholtzElectrostaticDomainBC: it holds a Dirichlet or Neumann BC object on each domain side
 * and routes the calls to it. The BC values are evaluated whenever getFaceFluxConstant is called, so time-dependent
 * values (FieldSolver captures its time by reference) enter through the right-hand side.
 */
class PetscHelmholtzElectrostaticDomainBC : public PetscHelmholtzDomainBC
{
public:
  /**
   * @brief Disallowed default constructor.
   */
  PetscHelmholtzElectrostaticDomainBC() = delete;

  /**
   * @brief Constructor. Creates a Dirichlet or Neumann BC object on each domain side.
   * @param[in] a_electrostaticBCs Electrostatic domain boundary conditions
   */
  PetscHelmholtzElectrostaticDomainBC(const ElectrostaticDomainBc& a_electrostaticBCs);

  /**
   * @brief Disallowed copy constructor.
   * @param[in] a_other Other BC
   */
  PetscHelmholtzElectrostaticDomainBC(const PetscHelmholtzElectrostaticDomainBC& a_other) = delete;

  /**
   * @brief Disallowed move constructor.
   * @param[in] a_other Other BC
   */
  PetscHelmholtzElectrostaticDomainBC(const PetscHelmholtzElectrostaticDomainBC&& a_other) = delete;

  /**
   * @brief Destructor.
   */
  ~PetscHelmholtzElectrostaticDomainBC();

  /**
   * @brief Disallowed copy assignment operator.
   * @param[in] a_other Other BC
   * @return This BC.
   */
  PetscHelmholtzElectrostaticDomainBC&
  operator=(const PetscHelmholtzElectrostaticDomainBC& a_other) = delete;

  /**
   * @brief Disallowed move assignment operator.
   * @param[in] a_other Other BC
   * @return This BC.
   */
  PetscHelmholtzElectrostaticDomainBC&
  operator=(const PetscHelmholtzElectrostaticDomainBC&& a_other) = delete;

  /**
   * @brief Define the BC objects on every domain side.
   * @param[in] a_dataLocation Data location
   * @param[in] a_eblg         Grids on each level
   * @param[in] a_probLo       Lower-left corner of the computational domain
   * @param[in] a_dx           Grid resolution on each level
   */
  void
  define(const Location::Cell       a_dataLocation,
         const Vector<EBLevelGrid>& a_eblg,
         const RealVect&            a_probLo,
         const Vector<Real>&        a_dx) override final;

  /**
   * @brief Get the face flux stencil from the BC object on the domain side.
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
                     const DataIndex&     a_din) const override final;

  /**
   * @brief Get the constant part of the face flux from the BC object on the domain side.
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
                      const DataIndex&     a_din) const override final;

protected:
  /**
   * @brief Electrostatic domain boundary conditions.
   */
  ElectrostaticDomainBc m_electrostaticBCs;

  /**
   * @brief BC objects on each domain side, created in the constructor.
   */
  std::map<ElectrostaticDomainBc::DomainSide, std::shared_ptr<PetscHelmholtzDomainBC>> m_bcObjects;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
