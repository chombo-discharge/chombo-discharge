/*!
  @file wire_wire_geometry.H
  @brief Sphere-sphere geometry (one dielectric sphere and one electrode)
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _WIRE_WIRE_GEOM_
#define _WIRE_WIRE_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines a wire-wire geometry
  
  @details 

  Class options
  -------------

      wire_wire_geometry.eps0           = 1.0
      wire_wire_geometry.turn_off_wire1 = false
      wire_wire_geometry.turn_off_wire2 = false
      wire_wire_geometry.wire1_center1  = -1.0 0.0 -100.0
      wire_wire_geometry.wire1_center2  = -1.0 0.0  100.0
      wire_wire_geometry.wire2_center1  = -1.0 0.0 -100.0
      wire_wire_geometry.wire2_center2  = -1.0 0.0  100.0
      wire_wire_geometry.wire1_radius   = 0.1
      wire_wire_geometry.wire2_radius   = 0.1
      wire_wire_geometry.wire1_live     = true
      wire_wire_geometry.wire2_live     = false


*/
class wire_wire_geometry : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  wire_wire_geometry();

  /*!
    @brief Destructor
  */
  ~wire_wire_geometry();

  /*!
    @brief Dump script for this class
  */
  static void dumpScript(std::string a_filename = "wire_wire_geometry.input");
};

#endif
