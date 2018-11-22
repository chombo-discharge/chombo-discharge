/*!
  @file   coaxial_packed_bed.H
  @brief  Coaxial cable geometry
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _COAXIAL_PACKED_BED_
#define _COAXIAL_PACKED_BED_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines a coaxial cable geometry 
*/
class coaxial_packed_bed : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  coaxial_packed_bed();

  /*!
    @brief Destructor
  */
  ~coaxial_packed_bed();
};

#endif
