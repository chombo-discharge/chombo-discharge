/*!
  @file rod_slab_geometry.H
  @brief Sphere-sphere geometry (one dielectric sphere and one electrode)
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _ROD_SLAB_GEOM_
#define _ROD_SLAB_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines one electrode sphere and a dielectric slab. Both can be turned off through class options.

  @details 
  
  Class options
  -------------

      rod_slab.eps0                      = 1                                # Background permittivity
      rod_slab.corner_curvatures         = 250E-6                           # Corner curvatures
      rod_slab.turn_off_electrode        = false                            # Turn on/off electrode
      rod_slab.turn_off_dielectric       = false                            # Turn on/off dielectric
      rod_slab.electrode_live            = true                             # Live electrode or not
      rod_slab.electrode_radius          = 5.E-3                            # Electrode inner radius
      rod_slab.electrode_center1         = 0.0 0.0 1E-2                     # Center 1
      rod_slab.electrode_center2         = 0.0 0.0 0E-2                     # Center 2
      rod_slab.dielectric_permittivity   = 4.0                              # Dielectric permittivity
      rod_slab.dielectric_corner_lo      = -1.0123E-2 -2.0123E-2 -2.0123E-2 # Low corner
      rod_slab.dielectric_corner_hi      =  1.0        2.0123E-2 -1.5123E-2 # High corner

*/
class rod_slab : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  rod_slab();

  /*!
    @brief Destructor
  */
  ~rod_slab();

protected:


};

#endif
