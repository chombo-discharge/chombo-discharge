/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Aerosol.H
  @brief  Aerosol geometry
  @author Robert Marskar
*/

#ifndef CD_Aerosol_H
#define CD_Aerosol_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

class Aerosol : public ComputationalGeometry {
public:

  Aerosol();
  ~Aerosol();
};

#include <CD_NamespaceFooter.H>

#endif
