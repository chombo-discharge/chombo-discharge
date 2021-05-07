/*!
  @file   double_rod.H
  @brief  Double rod geometry.
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _DOUBLE_ROD_GEOM_
#define _DOUBLE_ROD_GEOM_

#include "computational_geometry.H"

namespace ChomboDischarge {
  
  /*!
    @brief This is a class that defines two dielectric rods.
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
}

#endif
