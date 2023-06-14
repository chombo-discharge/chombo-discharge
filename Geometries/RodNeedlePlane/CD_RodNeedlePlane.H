/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodNeedlePlane.H
  @brief  Rod plane electrode geometry, with a needle protruding out of the rod.
  @author Fanny Skirbekk
*/

#ifndef CD_RodNeedlePlane_H
#define CD_RodNeedlePlane_H

//Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A class that defines a rod electrode with a needle out of the tip, and a plane (electrode) a set distance from the needle tip
*/
class RodNeedlePlane : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor. Puts object in usable state based on the input given in the .options-file
  */
  RodNeedlePlane();
  /*!
    @brief Destructor. Destroys the object and cleans up when the object goes out of scope.
  */
  virtual ~RodNeedlePlane();

protected:
  /*!
    @brief Defines a Rod electrode with a needle out of the tip as one electrode
  */
  void
  defineRodWNeedle();
  /*!
    @brief Defines a Rod electrode
  */
  void
  defineRod();
  /*!
    @brief Defines a Needle electrode
  */
  void
  defineNeedle();
  /*!
    @brief Defines a Plane electrode
  */
  void
  definePlane();
};

#include <CD_NamespaceFooter.H>

#endif
