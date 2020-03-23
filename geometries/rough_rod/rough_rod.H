/*!
  @file   rough_rod.H
  @brief  Rod geometry
  @author Robert Marskar
  @date   Oct. 2018
*/

#ifndef _ROUGH_ROD_GEOM_
#define _ROUGH_ROD_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines a wire-wire geometry
*/
class rough_rod : public computational_geometry {
public:

  /*!
    @brief Constructor. First rod is the electrode. 
  */
  rough_rod();

  /*!
    @brief Destructor
  */
  ~rough_rod();
};

#endif
