/*!
  @file   dcel_triangle.cpp
  @brief  Declaration of a polygon class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_triangle.H"
#include "dcel_edge.H"
#include "dcel_vert.H"

#if 1
#include <ParmParse.H>
#endif

dcel_triangle::dcel_triangle(){

}

dcel_triangle::~dcel_triangle(){

}

Real dcel_triangle::signed_distance(const RealVect a_x0) const {

  
  for (edge_iterator iter(static_cast<const dcel_poly*> (this)); iter.ok(); ++iter){
    const dcel_edge* edge = iter();

  }


  return 0.0;
}
