/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCSurfaceReaction.H
  @brief  Declaration of a class that encapsulates a set secondary emission at EBs and domain edges/faces
  @author Robert Marskar
*/

#ifndef CD_ItoKMCSurfaceReaction_H
#define CD_ItoKMCSurfaceReaction_H

// Std includes
#include <list>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Reaction class for describing secondary emissions due to a single outgoing particle in ItoKMCPhysics. 
      @details This class supports multiple emissions due to single impinging particles, and encapsulates all reactions A -> B, A -> B + C, etc.
    */
    class ItoKMCSurfaceReaction
    {
    public:
      /*!
	@brief Allowed, but you need to call the define function afterwards. 
      */
      ItoKMCSurfaceReaction() = default;

      /*!
	@brief Full constructor. Calls the define function. 
	@param[in] a_reactant   Index in the list of plasma or photon species that triggers the reaction.
	@param[in] a_products   List of plasma species potentially generated by the reaction.
	@param[in] a_efficiency Reaction efficiency
      */
      inline ItoKMCSurfaceReaction(const size_t            a_reactant,
                                   const std::list<size_t> a_products,
                                   const Real              a_efficiency) noexcept;

      /*!
	@brief Destructor. 
      */
      virtual inline ~ItoKMCSurfaceReaction() noexcept;

      /*!
	@brief Full define function.
	@param[in] a_reactant   Index in the list of plasma or photon species that triggers the reaction.
	@param[in] a_products   List of plasma species potentially generated by the reaction.
	@param[in] a_efficiency Reaction efficiency
      */
      inline void
      define(const size_t a_reactant, const std::list<size_t> a_products, const Real a_efficiency) noexcept;

      /*!
	@brief Get the reactant
      */
      inline const size_t&
      getReactant() const noexcept;

      /*!
	@brief Get the reactant
      */
      inline size_t&
      getReactant() noexcept;

      /*!
	@brief Get the plasma products.
      */
      inline const std::list<size_t>&
      getProducts() const noexcept;

      /*!
	@brief Get the plasma products.
      */
      inline std::list<size_t>&
      getProducts() noexcept;

      /*!
	@brief Get the reaction efficiency
      */
      inline const Real&
      getEfficiency() const noexcept;

      /*!
	@brief Get the reaction efficiencies. 
      */
      inline Real&
      getEfficiency() noexcept;

    protected:
      /*!
	@brief Source reactant index.
      */
      size_t m_reactant;

      /*!
	@brief Plasma product indices
      */
      std::list<size_t> m_products;

      /*!
	@brief Reaction efficiency
      */
      Real m_efficiency;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCSurfaceReactionImplem.H>

#endif
