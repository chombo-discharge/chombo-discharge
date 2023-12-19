/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCSurfaceReactions.H
  @brief  Declaration of a class that encapsulates a set secondary emission at EBs and domain edges/faces
  @author Robert Marskar
*/

#ifndef CD_ItoKMCSurfaceReactions_H
#define CD_ItoKMCSurfaceReactions_H

// Std includes
#include <random>

// Our includes
#include <CD_ItoKMCSurfaceReaction.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Reaction class for describing multiple pathways of secondary emissions due to a single outgoing particle in ItoKMCPhysics. 
      @details This class supports multiple emissions due to single impinging particles, and encapsulates all reactions A -> B, A -> B + C, etc.
    */
    class ItoKMCSurfaceReactions
    {
    public:
      /*!
	@brief Alias for list of products
      */
      using Products = std::list<size_t>;

      /*!
	@brief Allowed, but you need to call the define function afterwards. 
      */
      inline ItoKMCSurfaceReactions() noexcept;

      /*!
	@brief Destructor. 
      */
      virtual inline ~ItoKMCSurfaceReactions() noexcept;

      /*!
	@brief Add a reaction
      */
      inline void
      add(const ItoKMCSurfaceReaction& a_surfaceReaction) noexcept;

      /*!
	@brief Reset method for clearing class
      */
      inline void
      reset() noexcept;

      /*!
	@brief Get the reactant
      */
      inline const size_t&
      getReactant() const noexcept;

      /*!
	@brief Get the plasma products for each reaction
      */
      inline const std::vector<Products>&
      getProducts() const noexcept;

      /*!
	@brief Get the reaction efficiency
      */
      inline const std::vector<Real>&
      getEfficiencies() const noexcept;

      /*!
	@brief Get the discrete distribution for the reaction efficiencies
      */
      inline std::discrete_distribution<size_t>&
      getDistribution() const noexcept;

    protected:

      /*!
	@brief Is defined or not
      */
      bool m_isDefined;
      
      /*!
	@brief Reactant index.
      */
      size_t m_reactant;

      /*!
	@brief Plasma product indices
      */
      std::vector<Products> m_products;

      /*!
	@brief Reaction efficiencies
      */
      std::vector<Real> m_efficiencies;

      /*!
	@brief Discrete distribution -- useful when drawing a reaction.
      */
      mutable std::discrete_distribution<size_t> m_distribution;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCSurfaceReactionsImplem.H>

#endif
