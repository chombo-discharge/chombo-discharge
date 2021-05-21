/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_WireWire.H
  @brief  Declaration of a Wire-wire geometry
  @author Robert Marskar
*/

#ifndef CD_WireWire_H
#define CD_WireWire_H

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief This is a class that defines a wire-wire geometry
*/
class WireWire : public ComputationalGeometry {
public:
  
  WireWire();
  ~WireWire();

protected:

  void addWire(ParmParse& a_pp);
};

#include <CD_NamespaceFooter.H>

#endif
