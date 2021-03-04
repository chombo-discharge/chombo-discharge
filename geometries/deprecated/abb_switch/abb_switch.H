/*!
  @file   abb_switch.H
  @brief  ABB switch geometry
  @author Robert Marskar
  @date   Jan .2019
*/

#ifndef _ABB_SWITCH_
#define _ABB_SWITCH_

#include "computational_geometry.H"

/*!
  @brief This is a class that takes the switch PLY files and generates the geometry. 
  so I don't 

*/
class abb_switch : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  abb_switch();

  /*!
    @brief Destructor
  */
  ~abb_switch();

protected:


};

#endif
