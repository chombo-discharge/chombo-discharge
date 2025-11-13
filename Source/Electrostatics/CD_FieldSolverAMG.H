/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldSolverAMG.H
  @brief  Declaration of FieldSolverAMG
  @author Robert Marskar
*/

#ifndef CD_FieldSolverAMG_H
#define CD_FieldSolverAMG_H

// Our includes
#include <CD_MFHelmholtzPETSc.H>
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Implementation of FieldSolver which uses PETSc for solving the resulting system of equations. 
*/
class FieldSolverAMG : public FieldSolver
{
public:
  /*!
    @brief Weak constructor. 
    @details This sets default options. 
  */
  FieldSolverAMG();

  /*!
    @brief Disallowed copy constructor. 
    @param[in] a_other Other solver
  */
  FieldSolverAMG(const FieldSolverAMG& a_other) = delete;

  /*!
    @brief Disallowed move constructor.
    @param[in] a_other Other solver
  */
  FieldSolverAMG(const FieldSolverAMG&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment operator.
    @param[in] a_other Other solver
  */
  FieldSolverAMG&
  operator=(const FieldSolverAMG& a_other) = delete;

  /*!
    @brief Disallowed move assignment operator.
    @param[in] a_other Other solver
  */
  FieldSolverAMG&
  operator=(const FieldSolverAMG&& a_other) = delete;

  /*!
    @brief Constructor (does nothing)
  */
  virtual ~FieldSolverAMG();

  /*!
    @brief Solves Poisson equation onto a_phi using a_rho and a_sigma as right-hand sides. 
    @param[inout] a_potential Potential
    @param[in] a_rho Space charge density
    @param[in] a_sigma Surface charge density. 
    @param[in] a_zeroPhi Set a_potential to zero first. 
    @return True if we found a solution and false otherwise. 
    @note a_sigma must be defined on the gas phase. 
  */
  virtual bool
  solve(MFAMRCellData&       a_potential,
        const MFAMRCellData& a_rho,
        const EBAMRIVData&   a_sigma,
        const bool           a_zeroPhi = false) override;

  /*!
    @brief Parse all class options from command-line or input script. 
  */
  virtual void
  parseOptions() override;

  /*!
    @brief Parse runtime options from command line or input script. 
  */
  virtual void
  parseRuntimeOptions() override;

  /*!
    @brief Compute the cell-centered electric field. 
    @param[out] a_electricField Cell-centered electric field
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential) const override;

  /*!
    @brief Compute the face-centered electric field. 
    @param[out] a_electricField Face-centered electric field
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(MFAMRFluxData& a_E, const MFAMRCellData& a_potential) const override;

  /*!
    @brief Compute the cell-centered electric field on a specific phase. 
    @param[out] a_electricField Cell-centered electric field
    @param[in]  a_phase         Phase
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(EBAMRCellData&           a_E,
                       const phase::which_phase a_phase,
                       const MFAMRCellData&     a_potential) const override;

  /*!
    @brief Compute the face-centered electric field on a specific phase. 
    @param[out] a_electricField Face-centered electric field
    @param[in]  a_phase         Phase
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(EBAMRFluxData&           a_E,
                       const phase::which_phase a_phase,
                       const MFAMRCellData&     a_potential) const override;

  /*!
    @brief Cache state before regridding
    @param[in] a_lbase          Coarsest level which changes during regrid. 
    @param[in] a_oldFinestLevel Finest AMR level before regrid. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

  /*!
    @brief Regrid method. 
    @details This calls the parent regrid function first. We only need the overwrite because we have to signal that multigrid must be set up again. 
    @param[in] a_lmin           Coarsest level allowed to change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level before the regrid. 
    @details This linearly interpolates (with limiters) m_potential to the new grids and recomputes the electric field (from the interpolated potential). 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

  /*!
    @brief Registers required operators for AMR for this class. 
  */
  virtual void
  registerOperators() override;

  /*!
    @brief Set up multigrid solver
  */
  virtual void
  setupSolver() override;

  /*!
    @brief Set new permittivities for the multigrid solver. 
    @param[in] a_permittivityCell Permittivity on cell center
    @param[in] a_permittivityFace Permittivity on faces
    @param[in] a_permittivityEB   Permittivity on EB
  */
  virtual void
  setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                          const MFAMRFluxData& a_permittivityFace,
                          const MFAMRIVData&   a_permittivityEB) override;

protected:

  /*!
    @brief Set the permittivities
    @details This sets m_permittivityCell, m_permittivityFace, and m_permittivityEB and fills
    necessary ghost faces/cells. 
  */
  virtual void
  setPermittivities() override;
};

#include <CD_NamespaceFooter.H>

#endif
