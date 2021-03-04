/*!
  @file   wire_wire.H
  @brief  Declaration of a Wire-wire geometry
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _WIRE_WIRE_
#define _WIRE_WIRE_

#include "computational_geometry.H"

#include <ParmParse.H>

/*!
  @brief This is a class that defines a wire-wire geometry
*/
class wire_wire : public computational_geometry {
public:
  
  wire_wire();
  ~wire_wire();

protected:

  void addWire(ParmParse& a_pp);
};

#endif
