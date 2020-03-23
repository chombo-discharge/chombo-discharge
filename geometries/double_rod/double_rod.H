/*!
  @file   double_rod.H
  @brief  Double rod geometry.
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _DOUBLE_ROD_GEOM_
#define _DOUBLE_ROD_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines one electrode sphere and a dielectric slab. Both can be turned off through class options.

  @details 
  
  Class options
  -------------

      double_rod.eps0                      = 1                                # Background permittivity
      double_rod.turn_off_rod1             = false                            # Turn on/off electrode
      double_rod.turn_off_rod2             = false                            # Turn on/off electrode

      double_rod.rod1_live                 = true                             # Live electrode or not
      double_rod.rod2_live                 = true                             # Live electrode or not
      double_rod.rod1_radius               = 5.E-3                            # Rod radius
      double_rod.rod2_radius               = 5.E-3                            # Rod radius
      double_rod.rod1_center1              = -1E-2 0.0 1E-2                   # Center 1
      double_rod.rod1_center2              = -1E-2 0.0 0E-2                   # Center 2
      double_rod.rod2_center1              =  1E-2 0.0 1E-2                   # Center 1
      double_rod.rod2_center2              =  1E-2 0.0 0E-2                   # Center 2
*/
class double_rod : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  double_rod();

  /*!
    @brief Destructor
  */
  ~double_rod();

protected:


};

#endif
