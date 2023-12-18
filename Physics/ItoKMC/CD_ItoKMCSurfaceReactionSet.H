/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCSurfaceReactionSet.H
  @brief  Declaration of a class that encapsulates a set of secondary emission reactions
  @author Robert Marskar
*/

#ifndef CD_ItoKMCSurfaceReactionSet_H
#define CD_ItoKMCSurfaceReactionSet_H

// Our includes
#include <CD_ItoKMCSurfaceReaction.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Class for holding a set of surface reactions
    */
    class ItoKMCSurfaceReactionSet
    {
    public:
      /*!
	@brief Local enum for specifying surface type
      */
      enum class Surface
      {
        Dielectric,
        Electrode
      };

      /*!
	@brief Local enum for specifying species type
      */
      enum class Species
      {
        Plasma,
        Photon
      };
      /*!
	@brief Allowed, but you need to call the define function afterwards. 
      */
      ItoKMCSurfaceReactionSet() = default;

      /*!
	@brief Destructor. 
      */
      virtual inline ~ItoKMCSurfaceReactionSet() noexcept;

      /*!
	@brief Add a reaction
	@param[in] a_reaction Reaction to add
	@param[in] a_surface  Specification of surface type
	@param[in] a_species  Specification of species type
      */
      inline void
      addSurfaceReaction(const ItoKMCSurfaceReaction&            a_reaction,
                         const ItoKMCSurfaceReactionSet::Surface a_surface,
                         const ItoKMCSurfaceReactionSet::Species a_species) noexcept;

      /*!
	@brief Get the dielectric surface reactions where outgoing species is a plasma species
      */
      const std::vector<ItoKMCSurfaceReaction>&
      getDielectricPlasmaReactions() const noexcept;

      /*!
	@brief Get the dielectric surface reactions where outgoing species is a photon species
      */
      const std::vector<ItoKMCSurfaceReaction>&
      getDielectricPhotonReactions() const noexcept;

      /*!
	@brief Get the electrode surface reactions where outgoing species is a plasma species
      */
      const std::vector<ItoKMCSurfaceReaction>&
      getElectrodePlasmaReactions() const noexcept;

      /*!
	@brief Get the electrode surface reactions where outgoing species is a plasma species
      */
      const std::vector<ItoKMCSurfaceReaction>&
      getElectrodePhotonReactions() const noexcept;

    protected:
      /*!
	@brief Dielectric reaction set where outgoing species is a plasma species
      */
      std::vector<ItoKMCSurfaceReaction> m_dielectricPlasmaReactions;

      /*!
	@brief Dielectric reaction set where outgoing species is a photon species
      */
      std::vector<ItoKMCSurfaceReaction> m_dielectricPhotonReactions;

      /*!
	@brief Electrode reaction set where outgoing species is a plasma species
      */
      std::vector<ItoKMCSurfaceReaction> m_electrodePlasmaReactions;

      /*!
	@brief Electrode reaction set where outgoing species is a photon species
      */
      std::vector<ItoKMCSurfaceReaction> m_electrodePhotonReactions;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCSurfaceReactionSetImplem.H>

#endif
