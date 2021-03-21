/*!
  @file   dcel_BVH.H
  @brief  Declaration of some default lambdas for partitioning bounding volume hierarchies with dcel_face objects. 
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BVH_H_
#define _DCEL_BVH_H_

#include "BVH.H"
#include "dcel_face.H"

namespace dcel {
  
  template <class T, class BV>
  BVH::StopFunctionT<T, faceT<T>, BV>  defaultStopFunction;

  template <class T>
  BVH::PartitionFunctionT<faceT<T> >   defaultPartitionFunction;
  
  template <class T, class BV>
  BVH::PartitionFunctionT<faceT<T> >   partitionMinimumOverlap;

  template <class T, class BV>
  BVH::PartitionFunctionT<faceT<T> >   partitionSAH;

  template <class T, class BV>
  BVH::BVConstructorT<faceT<T>, BV>    defaultBVConstructor;
}

#include "dcel_BVHI.H"

#endif
