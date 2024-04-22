/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DischargeInceptionSpecies.H
  @brief  Simple species for DischargeInception test problem
  @author Robert Marskar
*/

#ifndef CD_DischargeInceptionSpecies_H
#define CD_DischargeInceptionSpecies_H

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace DischargeInception {

    /*!
      @brief Advection and diffused species for DischargeInceptionStepper
      @details The initial conditions and mobility/diffusion settings are initialized upon construction.
    */
    class DischargeInceptionSpecies : public CdrSpecies
    {
    public:
      /*!
	@brief Constructor
	@param[in] a_initialData Initial data
	@param[in] a_mobile      Mobile or not
	@param[in] a_diffusive   Diffusive or not
      */
      DischargeInceptionSpecies(const std::function<Real(const RealVect& a_pos)>& a_initialData,
                                const bool                                        a_mobile,
                                const bool                                        a_diffusive)
      {
        m_initialData  = a_initialData;
        m_isMobile     = a_mobile;
        m_isDiffusive  = a_diffusive;
        m_chargeNumber = 0;
        m_name         = "Ion distribution";
      }

      /*!
	@brief Destructor
      */
      virtual ~DischargeInceptionSpecies()
      {}

      /*!
	@brief Initial data. 
	@details This returns a super-Gaussian blob. 
	@param[in] a_pos  Position.
	@param[in] a_time Time
      */
      Real
      initialData(const RealVect a_pos, const Real a_time) const override
      {
        return m_initialData(a_pos);
      }

    protected:
      /*!
	@brief Initial ion density. 
      */
      std::function<Real(const RealVect& a_pos)> m_initialData;
    };
  } // namespace DischargeInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
