/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NoisePlane.H
  @brief  Rough sphere geometry
  @author Robert Marskar
*/

#ifndef CD_NoisePlane_H
#define CD_NoisePlane_H

#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief This is a class that a rough sphere (either electrode or dielectric)
*/
class NoisePlane : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor.
  */
  NoisePlane();

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~NoisePlane();
};

#include <CD_NamespaceFooter.H>

#endif
