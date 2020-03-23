/*!
  @file   rod_plane_profile.H
  @brief  Rod-plane geometry. The dielectric plane can have surface roughness
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _ROD_PLANE_PROFILE_GEOM_
#define _ROD_PLANE_PROFILE_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines one electrode sphere and a dielectric slab. Both can be turned off through class options.

  @details 
  
  Class options
  -------------

      rod_plane_profile.eps0                      = 1                                # Background permittivity
      rod_plane_profile.turn_off_electrode        = false                            # Turn on/off electrode
      rod_plane_profile.turn_off_dielectric       = false                            # Turn on/off dielectric

      rod_plane_profile.electrode_live            = true                             # Live electrode or not
      rod_plane_profile.electrode_radius          = 5.E-3                            # Electrode inner radius
      rod_plane_profile.electrode_center1         = 0.0 0.0 1E-2                     # Center 1
      rod_plane_profile.electrode_center2         = 0.0 0.0 0E-2                     # Center 2

      rod_plane_profile.dielectric_permittivity   = 4.0                              # Dielectric permittivity
      rod_plane_profile.plane_normal              = 0.0 1.0 0.0                      # Plane normal
      rod_plane_profile.plane_point               = 0.0 0.0 0.0                      # Point on plane
      rod_plane_profile.noise_amplitude           = 0.0                              # Surface roughness base amplitude
      rod_plane_profile.noise_frequency           = 1.0 1.0 1.0                      # Surface roughness frequency
      rod_plane_profile.noise_persistence         = 0.5                              # Drop-off factor between each octave
      rod_plane_profile.noise_reseed              = false                            # Reseed the RNG

*/
class rod_plane_profile : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  rod_plane_profile();

  /*!
    @brief Destructor
  */
  ~rod_plane_profile();

protected:


};

#endif
