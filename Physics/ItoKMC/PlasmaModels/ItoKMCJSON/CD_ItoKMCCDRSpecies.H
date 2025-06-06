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

// Std includes
#include <functional>

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
	@brief Function alias for e.g. initial data.
	@param[in] a_position Physical coordinates
	@param[in] a_time     Time
      */
      using FunctionXt = std::function<Real(const RealVect a_position, const Real a_time)>;

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

        m_initialData = [&](const RealVect a_pos, const Real a_time) -> Real {
          return 0.0;
        };
      }

      /*!
	@brief Destructor. Does nothing
      */
      virtual ~ItoKMCCDRSpecies() noexcept
      {}

      /*!
	@brief Returns initial data for this species.
	@param[in] a_pos Physical position
	@param[in] a_time Time
      */
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const override final
      {
        return m_initialData(a_pos, a_time);
      }

      /*!
	@brief Set the initial data function for this species
      */
      virtual void
      setInitialData(const FunctionXt& a_initialData)
      {
        m_initialData = a_initialData;
      }

    protected:
      /*!
	@brief Initial data function
      */
      FunctionXt m_initialData;
    };

  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
