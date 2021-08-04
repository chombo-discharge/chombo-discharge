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

// Chombo includes
#include <EBISBox.H>

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

  /*!
    @brief Compute the position (ignoring the "origin) of a Vof. 
    @param[in] a_location  Position within the grid cell (e.g. center, centroid, boundary).
    @param[in] a_vof       VolIndex to compute the position of
    @param[in] a_ebisbox   EBISBox
    @param[in] a_dx        Grid resolution
    @return Position of a face, ignoring the "origin" of the coordinate system since this is only used together with LeastSquares::displacement(...)
  */
  static RealVect position(const Location::Cell a_location,
			   const VolIndex&      a_vof,
			   const EBISBox&       a_ebisbox,
			   const Real&          a_dx);

  /*!
    @brief Compute the position (ignoring the "origin) of a FaceIndex
    @param[in] a_location  Position within the grid cell (e.g. face center, face centroid)
    @param[in] a_face      Face
    @param[in] a_ebisbox   EBISBox
    @param[in] a_dx        Grid resolution
    @return Position of a face, ignoring the "origin" of the coordinate system since this is only used together with LeastSquares::displacement(...)
  */
  static RealVect position(const Location::Face a_location,
			   const FaceIndex&     a_face,
			   const EBISBox&       a_ebisbox,
			   const Real&          a_dx);
};

#include <CD_NamespaceFooter.H>

#include <CD_LocationImplem.H>

#endif
