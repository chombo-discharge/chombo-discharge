/*!
  @file   rod_coarsen.H
  @brief  Declaration of a geometric coarsening mesh optimization for the mechshaft geometry
  @author Robert Marskar
  @date   June 2018
*/

#ifndef _ROD_COARSEN_
#define _ROD_COARSEN_

#include "geo_coarsener.H"

/*!
  @brief Geometric coarsening class for the mechshaft geometry. Useful for mesh reduction.
*/
class rod_coarsen : public geo_coarsener {
public:

  /*!
    @brief Constructor
  */
  rod_coarsen();

  /*!
    @brief Destructor
  */
  virtual ~rod_coarsen();
  
};
#endif
