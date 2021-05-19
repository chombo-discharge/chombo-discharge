/*!
  @file   electrode_array.H
  @brief  electrode array geometry.
  @author Robert Marskar
  @date   Dec. 2018
*/

#ifndef _ELECTRODE_ARRAY_GEOM_
#define _ELECTRODE_ARRAY_GEOM_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief This is a class that defines an array of electrodes.
*/
class Electrode_array : public ComputationalGeometry {
public:

  electrode_array();
  ~electrode_array();

protected:


};
#include "CD_NamespaceFooter.H"
#endif
