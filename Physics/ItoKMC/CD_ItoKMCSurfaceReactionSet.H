/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCSurfaceReactionSet.H
  @brief  Declaration of a class that encapsulates a set of secondary emission reactions
  @author Robert Marskar
*/

#ifndef CD_ITOKMCSURFACEREACTIONSET_H
#define CD_ITOKMCSURFACEREACTIONSET_H

// Our includes
#include <CD_ItoKMCSurfaceReactions.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /**
      @brief Class for holding a set of surface reactions
    */
    class ItoKMCSurfaceReactionSet
    {
    public:
      /**
	@brief Local enum for specifying surface type
      */
      enum class Surface
      {
        Dielectric,
        Electrode
      };

      /**
	@brief Local enum for specifying species type
      */
      enum class Species
      {
        Plasma,
        Photon
      };
      /**
	@brief Allowed, but you need to call the define function afterwards. 
      */
      ItoKMCSurfaceReactionSet() = default;

      /**
	@brief Destructor. 
      */
      virtual inline ~ItoKMCSurfaceReactionSet() noexcept;

      /**
	@brief Add a surface reaction to the set.
	@param[in] a_speciesIndex Index of the outgoing/impinging species that triggers the reaction.
	@param[in] a_reaction     Reaction to add.
	@param[in] a_surface      Whether this is a dielectric or electrode surface.
	@param[in] a_speciesType  Whether the outgoing species is a plasma or photon species.
      */
      inline void
      add(const size_t                            a_speciesIndex,
          const ItoKMCSurfaceReaction&            a_reaction,
          const ItoKMCSurfaceReactionSet::Surface a_surface,
          const ItoKMCSurfaceReactionSet::Species a_speciesType) noexcept;

      /**
	@brief Get the dielectric surface reactions where the outgoing species is a plasma species.
	@return Map from species index to reaction list for dielectric plasma reactions.
      */
      const std::map<size_t, ItoKMCSurfaceReactions>&
      getDielectricPlasmaReactions() const noexcept;

      /**
	@brief Get the dielectric surface reactions where the outgoing species is a photon species.
	@return Map from species index to reaction list for dielectric photon reactions.
      */
      const std::map<size_t, ItoKMCSurfaceReactions>&
      getDielectricPhotonReactions() const noexcept;

      /**
	@brief Get the electrode surface reactions where the outgoing species is a plasma species.
	@return Map from species index to reaction list for electrode plasma reactions.
      */
      const std::map<size_t, ItoKMCSurfaceReactions>&
      getElectrodePlasmaReactions() const noexcept;

      /**
	@brief Get the electrode surface reactions where the outgoing species is a photon species.
	@return Map from species index to reaction list for electrode photon reactions.
      */
      const std::map<size_t, ItoKMCSurfaceReactions>&
      getElectrodePhotonReactions() const noexcept;

    protected:
      /**
	@brief Dielectric reaction set where outgoing species is a plasma species
      */
      std::map<size_t, ItoKMCSurfaceReactions> m_dielectricPlasmaReactions;

      /**
	@brief Dielectric reaction set where outgoing species is a photon species
      */
      std::map<size_t, ItoKMCSurfaceReactions> m_dielectricPhotonReactions;

      /**
	@brief Electrode reaction set where outgoing species is a plasma species
      */
      std::map<size_t, ItoKMCSurfaceReactions> m_electrodePlasmaReactions;

      /**
	@brief Electrode reaction set where outgoing species is a photon species
      */
      std::map<size_t, ItoKMCSurfaceReactions> m_electrodePhotonReactions;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCSurfaceReactionSetImplem.H>

#endif
