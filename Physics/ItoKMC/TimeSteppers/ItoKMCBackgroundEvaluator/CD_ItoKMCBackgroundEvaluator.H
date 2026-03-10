/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCBackgroundEvaluator.H
  @brief  File containing a subclass of ItoKMCGodunovStepper.
  @author Robert Marskar
*/

#ifndef CD_ItoKMCBackgroundEvaluator_H
#define CD_ItoKMCBackgroundEvaluator_H

// Our includes
#include <CD_ItoKMCGodunovStepper.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Implementation of ItoKMCGodunovStepper that also evaluates the background field at every regrid/time step.
    */
    template <typename I = ItoSolver, typename C = CdrCTU, typename R = McPhoto, typename F = FieldSolverGMG>
    class ItoKMCBackgroundEvaluator : public ItoKMCGodunovStepper<I, C, R, F>
    {
    public:
      /*!
	@brief Constructor.
	@param[in] a_physics Physics kernel that will be used.
      */
      ItoKMCBackgroundEvaluator(RefCountedPtr<ItoKMCPhysics>& a_physics) noexcept;

      /*!
	@brief Post-initialization operations. Computes the background field.
      */
      virtual void
      postInitialize() noexcept override final;

      /*!
	@brief Post-regrid operations. Computes the background field.
      */
      virtual void
      postRegrid() noexcept override final;

      /*!
	@brief Post-checkpoint operations. Computes the background field.
      */
      virtual void
      postCheckpointSetup() noexcept override final;

      /*!
	@brief Print a new step report. This is the old one plus the maximum field change.
      */
      virtual void
      printStepReport() noexcept override final;

      /*!
	@brief Overriden advance method.
	@param[in] a_dt Time step
      */
      virtual Real
      advance(const Real a_dt) override;

      /*!
	@brief Compute a time step used for the advance method
	@note This function returns a zero value if space charge effects have set in, which lets Driver put plot files in the crash folder.
      */
      virtual Real
      computeDt() override;

    protected:
      /*!
	@brief Exit criterion for maximum field -- if the background field changes by this much the simulation is terminated.
      */
      Real m_maxFieldExitCrit;

      /*!
	@brief Exit criterion for relative field -- if the relative change in the background field in any cell changes by this much the simulation is terminated.
      */
      Real m_relFieldExitCrit;

      /*!
        @brief Initial time step size.
        Introduced to limit the initial particles from moving to far in the initial timestep.

        Default constructed to 0.0, which indicates normal operation.
        */
      Real m_maxInitialTimeStep;

      /*!
	@brief Change in max field.
      */
      Real m_maxFieldChange;

      /*!
	@brief Change in relative field.
      */
      Real m_relFieldChange;

      /*!
        @brief Total charge on the electrode, given by \int(eps*E dA)
      */
      Real m_electrodeCharge;

      /*!
        @brief Total charge that left through the electrode in the last time step.
      */
      Real m_ohmicCharge;

      /*!
        @brief Integrated optical emissions (phi)
      */
      Real m_sumOpticalPhi;

      /*!
        @brief Integrated optical emissions (source term)
      */
      Real m_sumOpticalSrc;

      /*!
        @brief Name of the optical solver
      */
      std::string m_opticalSolver;

      /*!
	@brief For storing the background field.
	@details This field is space- and surface-charge free and is filled by computeBackgroundField
      */
      MFAMRCellData m_backgroundField;

      /*!
	@brief Computes the background field. I.e., it allocates and fills m_backgroundField.
      */
      virtual void
      computeBackgroundField() noexcept;

      /*!
	@brief Evaluate the maximum change in the background field.
	@return Returns a pair: relative change, maximum change
      */
      [[nodiscard]]
      virtual std::pair<Real, Real>
      evaluateSpaceChargeEffects() noexcept;

      /*!
        @brief Compute total charge on electrode
        @details This returns \int(eps * E)*dA
      */
      [[nodiscard]]
      virtual Real
      integrateElectrodeSurfaceCharge() const noexcept;

      /*!
        @brief Compute the total ohmic current through the electrode.
        @details This returns the total charge that left the electrode in the last time step.
      */
      [[nodiscard]]
      virtual Real
      integrateOhmicCharge() const noexcept;

      /*!
        @brief Integrate the sum of the optical excitations.
        @details This returns the sum of the emissions, and the sum of the excitations.
      */
      [[nodiscard]]
      virtual std::pair<Real, Real>
      integrateOpticalExcitations() const noexcept;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCBackgroundEvaluatorImplem.H>

#endif
