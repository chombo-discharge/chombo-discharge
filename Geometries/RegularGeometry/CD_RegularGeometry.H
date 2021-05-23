/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RegularGeometry.H
  @brief  Regular geometry (no EBs)
  @author Robert Marskar
*/

#ifndef CD_RegularGeometry_H
#define CD_RegularGeometry_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief This class implements an empty geometry
*/
class RegularGeometry : public ComputationalGeometry {
public:

  RegularGeometry();
  ~RegularGeometry();
};

#include <CD_NamespaceFooter.H>

#endif
