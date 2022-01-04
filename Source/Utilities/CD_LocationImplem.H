/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LocationImplem.H
  @brief  Implementation of CD_Location.H
  @author Robert Marskar
*/

#ifndef CD_LocationImplem_H
#define CD_LocationImplem_H

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

inline
RealVect Location::position(const Location::Cell a_position,
			    const VolIndex&      a_vof,
			    const EBISBox&       a_ebisbox,
			    const Real&          a_dx){

  RealVect ret;
  
  switch(a_position){
  case Location::Cell::Center:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit;
    break;
  case Location::Cell::Centroid:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit + a_ebisbox.centroid(a_vof);
    break;
  case Location::Cell::Boundary:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit + a_ebisbox.bndryCentroid(a_vof);
    break;
  default:
    MayDay::Error("LeastSquares::position -- bad location input.");
    break;
  }

  ret *= a_dx;

  return ret;
}

inline
RealVect Location::position(const Location::Face a_position,
			    const FaceIndex&     a_face,
			    const EBISBox&       a_ebisbox,
			    const Real&          a_dx){
  RealVect ret;
  const IntVect iv = a_face.gridIndex(Side::Hi);
  
  switch(a_position){
  case Location::Face::Center:
    ret = (RealVect(iv) + 0.5*RealVect::Unit) - 0.5*BASISREALV(a_face.direction());
    break;
  case Location::Face::Centroid: {
    ret = (RealVect(iv) + 0.5*RealVect::Unit) - 0.5*BASISREALV(a_face.direction());

    const RealVect faceCentroid = a_ebisbox.centroid(a_face);
    
    for (int dir = 0; dir < SpaceDim; dir++){ // Direction a_face is undefined. 
      if(dir != a_face.direction()){
	ret[dir] += faceCentroid[dir];
      }
    }
    break;
  }
  default:
    MayDay::Error("LeastSquares::positon - bad location input.");
    break;
  }

  ret *= a_dx;

  return ret;
}

#include <CD_NamespaceFooter.H>

#endif
