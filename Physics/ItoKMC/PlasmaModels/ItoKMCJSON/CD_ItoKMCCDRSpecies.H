/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCCdrSpecies.H
  @brief  Declaration of a simple CdrSpecies derived class for usage in ItoKMC
  @author Robert Marskar
*/

#ifndef CD_ItoKMCCDRSpecies_H
#define CD_ItoKMCCDRSpecies_H

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Simple CdrSpecies class for usage with ItoKMC
    */
    class ItoKMCCDRSpecies : public CdrSpecies
    {
    public:
      /*!
	@brief Full constructor
	@param[in] a_name         Name of the species
	@param[in] a_chargeNumber Charge number
	@param[in] a_isMbile      Mobile or not
	@param[in] a_isDiffusive  Diffusive or not
      */
      ItoKMCCDRSpecies(const std::string a_name,
                       const int         a_chargeNumber,
                       const bool        a_isMobile,
                       const bool        a_isDiffusive)
      {
        m_name         = a_name;
        m_chargeNumber = a_chargeNumber;
        m_isMobile     = a_isMobile;
        m_isDiffusive  = a_isDiffusive;

        m_initialParticles.clear();
      }

      /*!
	@brief Destructor. Does nothing
      */
      virtual ~ItoKMCCDRSpecies() noexcept
      {}

      /*!
	@brief Initial data
      */
      Real
      initialData(const RealVect a_pos, const Real a_time) const override final
      {
        return 0.0;
      }
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
