/*!
  @file   double_rod.H
  @brief  Double rod geometry.
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _DOUBLE_ROD_GEOM_
#define _DOUBLE_ROD_GEOM_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief This is a class that defines two dielectric rods.
*/
class double_rod : public ComputationalGeometry {
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
#include "CD_NamespaceFooter.H"

#endif
