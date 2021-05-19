/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelBVH.H
  @brief  File which contains partitioners and lambdas for enclosing dcel_face in bounding volume heirarchies
  @author Robert Marskar
*/

#ifndef CD_DcelBVH_H
#define CD_DcelBVH_H

// Our includes
#include <CD_BVH.H>
#include <CD_DcelFace.H>
#include <CD_NamespaceHeader.H>
  
namespace Dcel {

  template <class T>
  using PrimitiveList = std::vector<std::shared_ptr<const Dcel::FaceT<T> > >;

  constexpr int primitivesPerLeafNode = 1;

  template <class T, class BV>
  BVH::StopFunctionT<T, FaceT<T>, BV> defaultStopFunction = [=](const BVH::NodeT<T, FaceT<T>, BV>& a_node){
    const auto& primitives = a_node.getPrimitives();
    const int depth        = a_node.getDepth();

    return primitives.size() <= primitivesPerLeafNode;
  };

  template <class T, class BV>
  BVH::BVConstructorT<FaceT<T>, BV> defaultBVConstructor = [](const PrimitiveList<T>& a_primitives){
    std::vector<Vec3T<T> > coordinates;

    for (const auto& f : a_primitives){
      const auto faceCoordinates = f->getAllVertexCoordinates();
      coordinates.insert(coordinates.end(), faceCoordinates.begin(), faceCoordinates.end());
    }

    return BV(coordinates);
  };
  
  template <class T>
  BVH::PartitionFunctionT<FaceT<T> > defaultPartitionFunction = [](const PrimitiveList<T>& a_primitives){

    auto lo = Vec3T<T>::max();
    auto hi = Vec3T<T>::min();
    for (const auto& p : a_primitives){
      lo = min(lo, p->getCentroid());
      hi = max(hi, p->getCentroid());
    }
    const auto delta   = (hi-lo);
    const int splitDir = delta.maxDir(true);

    PrimitiveList<T> sortedPrimitives(a_primitives);
  
    std::sort(sortedPrimitives.begin(), sortedPrimitives.end(),
	      [=](const std::shared_ptr<const FaceT<T> >& f1, const std::shared_ptr<const FaceT<T> >& f2) -> bool {
		return f1->getCentroid(splitDir) < f2->getCentroid(splitDir);
	      });

    const int splitIndex = (sortedPrimitives.size()-1)/2;

    PrimitiveList<T> lPrims(sortedPrimitives.begin(), sortedPrimitives.begin() + splitIndex+1);
    PrimitiveList<T> rPrims(sortedPrimitives.begin() + splitIndex + 1, sortedPrimitives.end());
  
    return std::make_pair(lPrims, rPrims);
  };

  template <class T, class BV>
  BVH::PartitionFunctionT<FaceT<T> > partitionMinimumOverlap = [](const PrimitiveList<T>& a_primitives){
    constexpr int DIM = 3;
    
    const int splitIndex    = (a_primitives.size() - 1)/2;

    T minOverlap = std::numeric_limits<T>::infinity();

    std::pair<PrimitiveList<T>, PrimitiveList<T> > ret;
  
    for (int dir = 0; dir < DIM; dir++){

      PrimitiveList<T> sortedPrims(a_primitives);
      std::sort(sortedPrims.begin(), sortedPrims.end(),
		[=](const std::shared_ptr<const FaceT<T> >& f1, const std::shared_ptr<const FaceT<T> >& f2){
		  return f1->getCentroid(dir) < f2->getCentroid(dir);
		});

      PrimitiveList<T> lPrims(sortedPrims.begin(), sortedPrims.begin() + splitIndex+1);
      PrimitiveList<T> rPrims(sortedPrims.begin() + splitIndex + 1, sortedPrims.end());

      const BV leftBV  = defaultBVConstructor<T, BV>(lPrims);
      const BV rightBV = defaultBVConstructor<T, BV>(rPrims);

      const T curOverlap = getOverlappingVolume(leftBV, rightBV);

      if (curOverlap < minOverlap){
	minOverlap = curOverlap;

	ret = std::make_pair(lPrims, rPrims);
      }
    }

    return ret;
  };

  template <class T, class BV>
  BVH::PartitionFunctionT<FaceT<T> > partitionSAH = [](const PrimitiveList<T>& a_primitives){
    constexpr int DIM   = 3; 
    constexpr int nBins = 16;
    constexpr T invBins = 1./nBins;
    constexpr T Ct      = 0.0;
    constexpr T Ci      = 1.0;

    const auto curBV   = defaultBVConstructor<T, BV>(a_primitives);
    const auto curArea = curBV.getArea();

    auto lo = Vec3T<T>::max();
    auto hi = Vec3T<T>::min();
    for (const auto& p : a_primitives){
      lo = min(lo, p->getCentroid());
      hi = max(hi, p->getCentroid());
    }
    const auto delta = (hi-lo)*invBins;

    T minCost = std::numeric_limits<T>::max();

    std::pair<PrimitiveList<T>, PrimitiveList<T> > ret;
  
    for (int dir = 0; dir < DIM; dir++){

      for (int ibin = 0; ibin <= nBins; ibin++){
	const Vec3T<T> pos = lo + T(1.0*ibin)*delta;

	PrimitiveList<T> lPrims;
	PrimitiveList<T> rPrims;
	
	for (const auto& p : a_primitives){
	  if(p->getCentroid()[dir] <= pos[dir]){
	    lPrims.emplace_back(p);
	  }
	  else{
	    rPrims.emplace_back(p);
	  }
	}

	const auto numLeft  = lPrims.size();
	const auto numRight = rPrims.size();
	
	if(numLeft == 0 || numRight == 0) continue;

	const BV bvLeft  = defaultBVConstructor<T, BV>(lPrims);
	const BV bvRight = defaultBVConstructor<T, BV>(rPrims);

	const T leftArea  = bvLeft.getArea();
	const T rightArea = bvRight.getArea();

	const T C = Ct + (leftArea/curArea)*Ci*numLeft + (rightArea/curArea)*Ci*numRight;

	if(C < minCost){
	  minCost = C;
	  ret     = std::make_pair(lPrims, rPrims);
	}
      }
    }

    return ret;
  };
}

#include <CD_NamespaceFooter.H>

#endif
