/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaPhotoReactionJSON.H
  @brief  Declaration of a class that encapsulates reactions for usage with CdrPlasmaJSON
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaPhotoReactionJSON_H
#define CD_CdrPlasmaPhotoReactionJSON_H

// Std includes
#include <map>
#include <vector>
#include <list>

// Our includes
#include <CD_CdrPlasmaPhysics.H>
#include <CD_NamespaceHeader.H>

namespace KMC {

/*!
  @brief Generic class for encapsulating reactive processes of the type A + B + Y -> C + D + X where A/C are "plasma species", B and D
  are "neutral species" and Y and X are photons involved in the reaction. In general, the left hand side is called "reactants" and
  the right hand side are "products". 
*/
  class Reaction
  {
  public:
    /*!
      @brief Disallowed constructor
    */
    Reaction() = delete;

    /*!
      @brief Copy constructor
    */
    Reaction(const Reaction&) = default;

    /*!
      @brief Disallowed move constructor
    */
    Reaction(const Reaction&&) = delete;

    /*!
      @brief Destructor
    */
    virtual ~Reaction();

    /*!
      @brief 
    */
  protected:

    std::list<int> m_plasmaReactants;
  };
} // namespace KMC

#include <CD_NamespaceFooter.H>
