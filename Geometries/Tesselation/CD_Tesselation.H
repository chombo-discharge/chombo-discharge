/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Tesselation.H
  @brief  Geometry for surface tesselation
  @author Robert Marskar
*/

#ifndef CD_Tesselation_H
#define CD_Tesselation_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A geometry consisting of a single polygon grid.
*/
class Tesselation : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor
  */
  Tesselation();

  /*!
    @brief Destructor
  */
  ~Tesselation();
};

#include <CD_NamespaceFooter.H>

#endif
