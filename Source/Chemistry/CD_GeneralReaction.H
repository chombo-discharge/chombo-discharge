/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_GeneralReaction.H
  @brief  Declaration of a "general" reaction type. 
  @author Robert Marskar
*/

#ifndef CD_GeneralReaction_H
#define CD_GeneralReaction_H

// Std includes
#include <map>
#include <vector>
#include <list>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Generic class for encapsulating reactive processes of the type A + B -> C 
  @todo Need a template type for indicating how we should piggyback ride other quantities (such as energy) on the reaction state. 
*/
template <size_t RTypes, size_t Ptypes, typename VType = std::vector, typename UType = long long>
class GeneralReaction
{
public:
  /*!
    @brief Disallowed constructor. Use the full constructor. 
  */
  GeneralReaction() = delete;

  /*!
    @brief Copy constructor (uses default). 
  */
  GeneralReaction(const Reaction&) = default;

  /*!
    @brief Disallowed move constructor
  */
  GeneralReaction(const Reaction&&) = delete;

  /*!
    @brief Full constructor
  */
  inline GeneralReaction(const std::array<std::list<size_t>, RTypes>& a_reactants,
                         const std::array<std::list<size_t>, PTypes>& a_products) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~GeneralReaction();

  /*!
    @brief Copy assignment operator
  */
  GeneralReaction&
  operator=(const GeneralReaction&) = default;

  /*!
    @brief Disallowed move operator
  */
  GeneralReaction&
  operator=(const GeneralReaction&&) = delete;

  /*!
    @brief Compute the propensity function for this reaction type
  */
  virtual inline Real
  propensity(const std::array<VType<UType>, RTypes>& a_state) const noexcept;

protected:
  /*!
    @brief Left-hand side reactants
  */
  std::array<std::list<size_t>, Rtypes> m_reactants;

  /*!
    @brief Left-hand side reactants
  */
  std::array<std::list<size_t>, Ptypes> m_products;

  /*!
    @brief State change due to firing of exactly one reaction. 
    @details First index is species and second index is change (e.g., lose one particle per reaction);
  */
  std::map<size_t, int> m_stateChange;

  /*!
    @brief Compute state change
  */
  inline void
  computeStateChange() noexcept;
};

#include <CD_NamespaceFooter.H>
