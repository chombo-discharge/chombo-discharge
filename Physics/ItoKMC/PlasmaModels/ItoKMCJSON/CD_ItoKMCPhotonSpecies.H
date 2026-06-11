/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCPhotonSpecies.H
  @brief  Declaration of a simple PhotonSpecies derived class for usage in ItoKMC
  @author Robert Marskar
*/

#ifndef CD_ITOKMCPHOTONSPECIES_H
#define CD_ITOKMCPHOTONSPECIES_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /**
      @brief Simple RtSpecies class for usage with ItoKMCJSON
    */
    class ItoKMCPhotonSpecies : public RtSpecies
    {
    public:
      /**
	@brief Full constructor
	@param[in] a_name          Name of the species
	@param[in] a_kappaFunction Absorption function
      */
      ItoKMCPhotonSpecies(const std::string& a_name, const std::function<Real(const RealVect& a_pos)>& a_kappaFunction)
        : m_kappa(a_kappaFunction)
      {
        m_name = a_name;
      }

      /**
	@brief Destructor. Does nothing
      */
      ~ItoKMCPhotonSpecies() noexcept override = default;

      ItoKMCPhotonSpecies(const ItoKMCPhotonSpecies&) = default;
      ItoKMCPhotonSpecies&
      operator=(const ItoKMCPhotonSpecies&)      = default;
      ItoKMCPhotonSpecies(ItoKMCPhotonSpecies&&) = default;
      ItoKMCPhotonSpecies&
      operator=(ItoKMCPhotonSpecies&&) = default;

      /**
	@brief Get absorption coefficient at a physical position.
	@param[in] a_pos Physical position.
	@return Absorption coefficient (inverse absorption length) at a_pos.
      */
      Real
      getAbsorptionCoefficient(const RealVect& a_pos) const final
      {
        return m_kappa(a_pos);
      }

    protected:
      /**
	@brief Absorption function
      */
      std::function<Real(const RealVect& a_position)> m_kappa;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
