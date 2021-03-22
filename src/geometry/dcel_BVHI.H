/*!
  @file   dcel_BVHI.H
  @brief  Implementation of dcel_BVH.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BVHI_H_
#define _DCEL_BVHI_H_

#include "dcel_BVH.H"

template <class T>
using PrimitiveList = std::vector<std::shared_ptr<const faceT<T> > >; 

namespace dcel {

  template <class T, class BV>
  BVH::StopFunctionT<T, faceT<T>, BV> defaultStopFunction = [](const BVH::NodeT<T, faceT<T>, BV>& a_node){
    const auto& primitives = a_node.getPrimitives();
    const int depth        = a_node.getDepth();

    return primitives.size() < 2;
  };
  
  template <class T>
  BVH::PartitionFunctionT<faceT<T> > defaultPartitionFunction = [](const PrimitiveList<T>& a_primitives){
    Vec3T<T> minCorner = a_primitives.front()->getCentroid();
    Vec3T<T> maxCorner = a_primitives.front()->getCentroid();

    for(const auto& p : a_primitives){
      minCorner = min(minCorner, p->getCentroid());
      maxCorner = max(minCorner, p->getCentroid());
    }

    const Vec3T<T> delta = maxCorner - minCorner;
    const int splitDir   = delta.maxDir(true);

    PrimitiveList<T> sortedPrimitives(a_primitives);
    auto comparator = [splitDir](const std::shared_ptr<const faceT<T> >& f1, const std::shared_ptr<const faceT<T> >& f2){
      return f1->getCentroid()[splitDir] < f2->getCentroid()[splitDir];
    };
  
    std::sort(sortedPrimitives.begin(), sortedPrimitives.end(), comparator);

    const int splitIndex = (sortedPrimitives.size()-1)/2;

    PrimitiveList<T> lPrims(sortedPrimitives.begin(), sortedPrimitives.begin() + splitIndex+1);
    PrimitiveList<T> rPrims(sortedPrimitives.begin() + splitIndex + 1, sortedPrimitives.end());
  
    return std::make_pair(lPrims, rPrims);
  };

  template <class T, class BV>
  BVH::PartitionFunctionT<faceT<T> > partitionMinimumOverlap = [](const PrimitiveList<T>& a_primitives){
    const int splitIndex    = (a_primitives.size() - 1)/2;

    T minOverlap = std::numeric_limits<T>::infinity();

    std::pair<PrimitiveList<T>, PrimitiveList<T> > ret;

    constexpr int DIM = 3;
  
    for (int dir = 0; dir < DIM; dir++){

      PrimitiveList<T> sortedPrims(a_primitives);
      std::sort(sortedPrims.begin(), sortedPrims.end(),
		[dir](const std::shared_ptr<const faceT<T> >& f1, const std::shared_ptr<const faceT<T> >& f2){
		  return f1->getCentroid()[dir] < f2->getCentroid()[dir];
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
  BVH::PartitionFunctionT<faceT<T> > partitionSAH = [](const PrimitiveList<T>& a_primitives){
    constexpr int DIM   = 3; 
    constexpr int nBins = 50;
    constexpr T invBins = 1./nBins;
    constexpr T Ct      = 1.0;
    constexpr T Ci      = 1.0;

    const auto curBV   = defaultBVConstructor<T, BV>(a_primitives);
    const auto curArea = curBV.getArea();

    Vec3T<T> lo = a_primitives.front()->getCentroid();
    Vec3T<T> hi = a_primitives.front()->getCentroid();
    for (const auto& p : a_primitives){
      lo = min(lo, p->getCentroid());
      hi = max(hi, p->getCentroid());
    }
  
    const Vec3T<T> delta = (hi-lo)*invBins;

    T minCost = std::numeric_limits<T>::max();

    std::pair<PrimitiveList<T>, PrimitiveList<T> > ret;
  
    for (int dir = 0; dir < DIM; dir++){

      // Sort primitives along dir direction
      PrimitiveList<T> sortedPrimitives = a_primitives;
      std::sort(sortedPrimitives.begin(), sortedPrimitives.end(),
		[dir](const std::shared_ptr<const faceT<T> >& f1, const std::shared_ptr<const faceT<T> >& f2){
		  return f1->getCentroid()[dir] < f2->getCentroid()[dir];
		});

      // Look through spatial bins and determine lowest cost.
      for (int ibin = 0; ibin <= nBins; ibin++){
	const Vec3T<T> pos = lo + T(1.0*ibin)*delta;

	PrimitiveList<T> lPrims;
	PrimitiveList<T> rPrims;

	for (const auto& p : sortedPrimitives){
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

	const auto bvLeft  = defaultBVConstructor<T, BV>(lPrims);
	const auto bvRight = defaultBVConstructor<T, BV>(rPrims);

	const auto leftArea  = bvLeft.getArea();
	const auto rightArea = bvRight.getArea();

	const auto C = Ct + (leftArea/curArea)*Ci*numLeft + (rightArea/curArea)*Ci*numRight;

	if(C < minCost){
	  minCost = C;
	  ret     = std::make_pair(lPrims, rPrims);
	}
      }
    }

    return ret;
  };
  

  template <class T, class BV>
  BVH::BVConstructorT<faceT<T>, BV> defaultBVConstructor = [](const PrimitiveList<T>& a_primitives){
    std::vector<Vec3T<T> > coordinates;

    for (const auto& f : a_primitives){
      const auto faceCoordinates = f->getAllVertexCoordinates();
      coordinates.insert(coordinates.end(), faceCoordinates.begin(), faceCoordinates.end());
    }

    return BV(coordinates);
  };

}

  

#endif
