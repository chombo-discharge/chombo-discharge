/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Copyright © 2023 SINTEF Energy Research
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodNeedleDisk.H
  @brief  Rod plane electrode geometry, with a needle protruding out of the rod.
  @author Fanny Skirbekk
  @author Robert Marskar
*/

#ifndef CD_RodNeedleDisk_H
#define CD_RodNeedleDisk_H

//Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A class that defines a rod electrode with a needle out of the tip, and a plane (electrode) a set distance from the needle tip.
  @details This geometry is an axi-symmetric geometry oriented along one of the axis. All input points are in 1D along the specified axis.
*/
class RodNeedleDisk : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor. Puts object in usable state based on the input given in the .options-file
  */
  RodNeedleDisk() noexcept;
  /*!
    @brief Destructor. Destroys the object and cleans up when the object goes out of scope.
  */
  virtual ~RodNeedleDisk() noexcept;

protected:
  /*!
    @brief Define the rod+needle electrode
  */
  void
  defineRodNeedle() noexcept;

  /*!
    @brief Defines a Plane electrode
  */
  void
  defineDisk();
};

#include <CD_NamespaceFooter.H>

#endif
