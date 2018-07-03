/*!
  @file   rod_sphere.H
  @brief  Rod-sphere geometry (one dielectric sphere and one rod electrode)
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _ROD_SPHERE_GEOM_
#define _ROD_SPHERE_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines one electrode sphere and a dielectric slab. Both can be turned off through class options.

  @details 
  
  Class options
  -------------

      rod_sphere.eps0                      = 1                                # Background permittivity
      rod_sphere.turn_off_electrode        = false                            # Turn on/off electrode
      rod_sphere.turn_off_dielectric       = false                            # Turn on/off dielectric
      rod_sphere.electrode_live            = true                             # Live electrode or not
      rod_sphere.electrode_radius          = 5.E-3                            # Electrode inner radius
      rod_sphere.electrode_center1         = 0.0 0.0 1E-2                     # Center 1
      rod_sphere.electrode_center2         = 0.0 0.0 0E-2                     # Center 2
      rod_sphere.dielectric_permittivity   = 4.0                              # Dielectric permittivity
      rod_sphere.dielectric_center         = 0.0 0.0 0.0                      # Dielectric center
      rod_sphere.dielectric_radius         = 1.0                              # Dielectric radius

*/
class rod_sphere : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  rod_sphere();

  /*!
    @brief Destructor
  */
  ~rod_sphere();

protected:


};

#endif
