/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ElectrodeArray.H
  @brief  Electrode array geometry.
  @author Robert Marskar
*/

#ifndef CD_ElectrodeArray_H
#define CD_ElectrodeArray_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief This is a class that defines an array of electrodes.
*/
class ElectrodeArray : public ComputationalGeometry {
public:

  ElectrodeArray();
  ~ElectrodeArray();
};

#include <CD_NamespaceFooter.H>

#endif
