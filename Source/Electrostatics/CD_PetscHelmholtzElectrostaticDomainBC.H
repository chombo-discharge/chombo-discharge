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
#include <utility>

// Our includes
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_ElectrostaticDomainBc.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Domain boundary conditions for PetscHelmholtzOp, taken from ElectrostaticDomainBc.
 * @details This mirrors EBHelmholtzElectrostaticDomainBC: it holds a Dirichlet or Neumann BC object on each domain side
 * and routes the calls to it. The BC values are evaluated whenever getBoundaryData is called, so time-dependent
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
   * @param[in] a_petscGrid PETSc grid
   * @param[in] a_probLo    Lower-left corner of the computational domain
   * @param[in] a_dx        Grid resolution on each level
   */
  void
  define(const RefCountedPtr<PetscGrid>& a_petscGrid,
         const RealVect&                 a_probLo,
         const Vector<Real>&             a_dx) override final;

  /**
   * @brief Get the face flux stencil from the BC object on the domain side.
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
                     const Side::LoHiSide a_side) const override final;

  /**
   * @brief Get the boundary data from the BC object on the domain side.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cell adjacent to the domain face
   * @param[in] a_dir   Coordinate direction of the face
   * @param[in] a_side  Domain side of the face
   * @return The boundary data.
   */
  Real
  getBoundaryData(const int            a_phase,
                  const int            a_level,
                  const DataIndex&     a_din,
                  const VolIndex&      a_vof,
                  const int            a_dir,
                  const Side::LoHiSide a_side) const override final;

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
