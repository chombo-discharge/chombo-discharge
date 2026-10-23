/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzElectrostaticEBBC.H
 * @brief  Declaration of the electrode boundary conditions on the EB faces for PetscHelmholtzOp.
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_PETSCHELMHOLTZELECTROSTATICEBBC_H
#define CD_PETSCHELMHOLTZELECTROSTATICEBBC_H

// Our includes
#include <CD_PetscHelmholtzDirichletEBBC.H>
#include <CD_ElectrostaticEbBc.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Electrode boundary conditions on the EB faces for PetscHelmholtzOp, taken from ElectrostaticEbBc.
 * @details This is a Dirichlet condition whose value on an EB face is the potential of the closest electrode, as in
 * MFHelmholtzElectrostaticEBBC. The electrode potentials capture the time by reference, so time-dependent potentials
 * enter through the right-hand side.
 */
class PetscHelmholtzElectrostaticEBBC : public PetscHelmholtzDirichletEBBC
{
public:
  /**
   * @brief Disallowed default constructor.
   */
  PetscHelmholtzElectrostaticEBBC() = delete;

  /**
   * @brief Full constructor. Must subsequently set the order and the weight, and call define.
   * @param[in] a_electrostaticBCs Electrode boundary conditions
   */
  PetscHelmholtzElectrostaticEBBC(const ElectrostaticEbBc& a_electrostaticBCs) noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~PetscHelmholtzElectrostaticEBBC() noexcept;

  /**
   * @brief Get the boundary data, i.e. the potential of the closest electrode.
   * @param[in] a_phase Phase
   * @param[in] a_level Grid level
   * @param[in] a_din   Grid index on the level
   * @param[in] a_vof   Cut cell
   * @return The potential of the closest electrode, at the EB face centroid.
   */
  virtual Real
  getBoundaryData(const int        a_phase,
                  const int        a_level,
                  const DataIndex& a_din,
                  const VolIndex&  a_vof) const noexcept override;

protected:
  /**
   * @brief Electrode boundary conditions.
   */
  ElectrostaticEbBc m_electrostaticBCs;

  /**
   * @brief Get the potential of the electrode closest to a position.
   * @param[in] a_pos Position
   * @return The potential of the closest electrode.
   */
  virtual Real
  getElectrodePotential(const RealVect& a_pos) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
