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
class RodNeedlePlane : public ComputationalGeometry {
public:
  RodNeedlePlane();
  ~RodNeedlePlane();
protected:
  void defineRodWNeedle();
  void defineRod();
  void defineNeedle();
  void definePlane();
};

#include <CD_NamespaceFooter.H>

#endif
