/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCPhotonSpecies.H
  @brief  Declaration of a simple PhotonSpecies derived class for usage in ItoKMC
  @author Robert Marskar
*/

#ifndef CD_ItoKMCPhotonSpecies_H
#define CD_ItoKMCPhotonSpecies_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Simple RtSpecies class for usage with ItoKMCJSON
    */
    class ItoKMCPhotonSpecies : public RtSpecies
    {
    public:
      /*!
	@brief Full constructor
	@param[in] a_name          Name of the species
	@param[in] a_kappaFunction Absorption function
      */
      ItoKMCPhotonSpecies(const std::string a_name, const std::function<Real(const RealVect& a_pos)>& a_kappaFunction)
      {
        m_name  = a_name;
        m_kappa = a_kappaFunction;
      }

      /*!
	@brief Destructor. Does nothing
      */
      virtual ~ItoKMCPhotonSpecies() noexcept {}

      /*!
	@brief Initial data
      */
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const override final
      {
        return m_kappa(a_pos);
      }

    protected:
      /*!
	@brief Absorption function
      */
      std::function<Real(const RealVect& a_position)> m_kappa;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
