/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_GeneralReactionImplem.H
  @brief  Implementation of CD_GeneralReaction.H
  @author Robert Marskar
*/

#ifndef CD_GeneralReactionImplem_H
#define CD_GeneralReactionImplem_H

// Our includes
#include <CD_GeneralReaction.H>
#include <CD_NamespaceHeader.H>

template <size_t RTypes, size_t PTypes>
inline
GeneralReaction<RTypes, PTypes>::GeneralReaction(const std::array<std::list<size_t>, RTypes>& a_reactants,
						 const std::array<std::list<size_t>, PTypes>& a_products) noexcept {

  m_reactants = a_reactants;
  m_products = a_products;
}

template <size_t RTypes, size_t PTypes>
inline
GeneralReaction<RTypes, PTypes>::~GeneralReaction() {

}

#include <CD_NamespaceFooter.H>
