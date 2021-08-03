/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Location.H
  @brief  Declaration of cell positions.
  @author Robert Marskar
*/

#ifndef CD_Location_H
#define CD_Location_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  Enums for specying where cell/face data is located
*/
namespace Location {
  enum class Cell{
    Center,
    Centroid,
    Boundary
  };

  enum class Face {
    Center,
    Centroid
  };
};

#include <CD_NamespaceFooter.H>

#endif
