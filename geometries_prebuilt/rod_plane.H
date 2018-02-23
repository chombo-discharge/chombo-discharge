/*!
  @file   rod_plane_geometry.H
  @brief  Rod-plane geometry. The dielectric plane can have surface roughness
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _ROD_PLANE_GEOM_
#define _ROD_PLANE_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines one electrode sphere and a dielectric slab. Both can be turned off through class options.

  @details 
  
  Class options
  -------------

      rod_plane.eps0                      = 1                                # Background permittivity
      rod_plane.turn_off_electrode        = false                            # Turn on/off electrode
      rod_plane.turn_off_dielectric       = false                            # Turn on/off dielectric

      rod_plane.electrode_live            = true                             # Live electrode or not
      rod_plane.electrode_radius          = 5.E-3                            # Electrode inner radius
      rod_plane.electrode_center1         = 0.0 0.0 1E-2                     # Center 1
      rod_plane.electrode_center2         = 0.0 0.0 0E-2                     # Center 2

      rod_plane.dielectric_permittivity   = 4.0                              # Dielectric permittivity
      rod_plane.plane_normal              = 0.0 1.0 0.0                      # Plane normal
      rod_plane.plane_point               = 0.0 0.0 0.0                      # Point on plane
      rod_plane.noise_amplitude           = 0.0                              # Surface roughness base amplitude
      rod_plane.noise_frequency           = 1.0 1.0 1.0                      # Surface roughness frequency
      rod_plane.noise_persistence         = 0.5                              # Drop-off factor between each octave
      rod_plane.noise_reseed              = false                            # Reseed the RNG

*/
class rod_plane : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  rod_plane();

  /*!
    @brief Destructor
  */
  ~rod_plane();

protected:


};

#endif
