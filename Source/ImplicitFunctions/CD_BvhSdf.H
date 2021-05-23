/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BvhSdf.H
  @brief  Declaration of an signe distance function class that gets its value function from a DCEL surface Tesselation
  @author Robert Marskar
*/

#ifndef CD_BvhSdf_H
#define CD_BvhSdf_H

// Std includes
#include <memory>
#include <chrono>

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_BVH.H>
#include <CD_DcelFace.H>
#include <CD_NamespaceHeader.H>

using namespace Dcel;
using namespace BVH;

template <class T, class BV>
class BvhSdf : public BaseIF {
public:

  using Vec3  = Vec3T<T>;
  using Face  = FaceT<T>;
  using Node  = NodeT<T, Face, BV>;
  

  //  BvhSdf() = delete;
  BvhSdf();
  BvhSdf(const std::shared_ptr<Node>& a_root, const bool a_flipInside);
  BvhSdf(const BvhSdf& a_object);
  ~BvhSdf();

  Real value(const RealVect& a_point) const override;
  
  BaseIF* newImplicitFunction() const override;

protected:
  std::shared_ptr<Node> m_root;

  mutable long m_numCalled;

  mutable std::chrono::duration<double> m_timespan;
  
  bool m_flipInside;
};

#include <CD_NamespaceFooter.H>

#include <CD_BvhSdfImplem.H>

#endif
