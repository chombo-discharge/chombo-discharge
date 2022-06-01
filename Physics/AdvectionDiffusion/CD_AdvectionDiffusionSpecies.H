/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AdvectionDiffusionSpecies.H
  @brief  Simple species for AdvectionDiffusion test problem
  @author Robert Marskar
*/

#ifndef CD_AdvectionDiffusionSpecies_H
#define CD_AdvectionDiffusionSpecies_H

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace AdvectionDiffusion {

    /*!
      @brief Advection and diffused species for AdvectionDiffusionStepper
      @details The initial conditions and mobility/diffusion settings are initialized upon construction.
    */
    class AdvectionDiffusionSpecies : public CdrSpecies
    {
    public:
      /*!
	@brief Constructor
      */
      AdvectionDiffusionSpecies();

      /*!
	@brief Destructor
      */
      virtual ~AdvectionDiffusionSpecies();

      /*!
	@brief Initial data. 
	@details This returns a super-Gaussian blob. 
	@param[in] a_pos  Position.
	@param[in] a_time Time
      */
      Real
      initialData(const RealVect a_pos, const Real a_time) const override;

    protected:
      /*!
	@brief Gaussian blob radius
      */
      Real m_blobRadius;

      /*!
	@brief Gaussian blob amplitude
      */
      Real m_blobAmplitude;

      /*!
	@brief Gaussian blob center
      */
      RealVect m_blobCenter;
    };
  } // namespace AdvectionDiffusion
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
