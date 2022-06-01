/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Aerosol.H
  @brief  Declaration of an aerosol geometry
  @author Robert Marskar
*/

#ifndef CD_Aerosol_H
#define CD_Aerosol_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Simple class for making a bunch of (manually constructored) dielectric sphere
*/
class Aerosol : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor
  */
  Aerosol();

  /*!
    @brief Destructor
  */
  ~Aerosol();
};

#include <CD_NamespaceFooter.H>

#endif
