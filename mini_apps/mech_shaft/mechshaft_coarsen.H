/*!
  @file   mechshaft_coarsen.H
  @brief  Declaration of a geometric coarsening mesh optimization for the mechshaft geometry
  @author Robert Marskar
  @date   June 2018
*/

#ifndef _MECHSHAFT_COARSEN_
#define _MECHSHAFT_COARSEN_

#include "geo_coarsener.H"

/*!
  @brief Geometric coarsening class for the mechshaft geometry. Useful for mesh reduction.
*/
class mechshaft_coarsen : public geo_coarsener {
public:

  /*!
    @brief Constructor
  */
  mechshaft_coarsen();

  /*!
    @brief Destructor
  */
  virtual ~mechshaft_coarsen();
  
};
#endif
