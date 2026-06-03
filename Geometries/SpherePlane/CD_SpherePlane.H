/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SpherePlane.H
  @brief  Sphere-plane gap geometry
  @author Robert Marskar
*/

#ifndef CD_SpherePlane_H
#define CD_SpherePlane_H

// our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Sphere-plane gap geometry.
  @details The geometry is always oriented along +y, in both 2D and 3D. 
*/
class SpherePlane : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor
  */
  SpherePlane();

  /*!
    @brief Destructor
  */
  virtual ~SpherePlane();
};

#include <CD_NamespaceFooter.H>

#endif
