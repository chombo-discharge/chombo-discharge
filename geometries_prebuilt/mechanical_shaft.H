/*!
  @file   mechanical_shaft.H
  @brief  Declaration of the mechanical shaft geometry (ABB problem)
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _MECHANICAL_SHAFT_
#define _MECHANICAL_SHAFT_

#include "computational_geometry.H"

/*!
  @brief Mechanical shaft class

  @details This class consist of a hexagonal dielectric rod surrounded by a hollow cylinder. The 6 corners on the rod are rounded, 
  and so are the internal and external edges of the hollow cylinder. 

  Class options
  -------------

      mechanical_shaft.eps0                      = 1               # Background permittivity
      mechanical_shaft.corner_curvatures         = 250E-6          # Corner curvatures
      mechanical_shaft.turn_off_electrode        = false           # Turn on/off electrode
      mechanical_shaft.turn_off_dielectric       = false           # Turn on/off dielectric
      mechanical_shaft.electrode_live            = true            # Live electrode or not
      mechanical_shaft.electrode_inner_radius    = 7.E-3           # Electrode inner radius
      mechanical_shaft.electrode_outer_radius    = 1.E-2           # Electrode outer radius
      mechanical_shaft.electrode_height          = 1.E-2           # Electrode height
      mechanical_shaft.electrode_center          = 0.0 0.0 7.E-3   # Electrode center
      mechanical_shaft.dielectric_permittivity   = 4.0             # Dielectric permittivity
      mechanical_shaft.dielectric_rod_sides      = 6               # Number of sides
      mechanical_shaft.dielectric_rod_radius     = 5E-3            # Dielectric rod radius
      mechanical_shaft.dielectric_rod_length     = 1.0             # Dielectric rod length


  @note This class only works in 3D. 
*/
class mechanical_shaft : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  mechanical_shaft();

  /*!
    @brief Destructor
  */
  ~mechanical_shaft();
};

#endif
