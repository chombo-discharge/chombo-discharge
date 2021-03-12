/*!
  @file   dcel_if.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface tesselation
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_IF_
#define _DCEL_IF_

#include "dcel_BVH.H"
#include "dcel_mesh.H"

#include <memory>

#include <BaseIF.H>

using namespace dcel;

template <class T>
using SortingFunctor = std::function<bool(const faceT<T>&, const faceT<T> &, const int)>;

template <class T>
using CostFunctor    = std::function<std::pair<double, int>(const std::vector<std::shared_ptr<faceT<T> > >&, const int)>;

template <class T>
SortingFunctor<T> defaultSortingFunctor = [](const faceT<T>& f1, const faceT<T>& f2, const int dir){
  return f1.getCentroid()[dir] < f1.getCentroid()[dir];
};

template <class T>
CostFunctor<T> defaultCostFunctor = [](const std::vector<std::shared_ptr<faceT<T> > >& a_faces, const int a_dir){
  return std::make_pair(0,a_dir);
};


/*!
  @brief Class that gets it's value function from a surface tesselation 
*/
template <class T, class BV>
class dcel_if : public BaseIF {
public:

  using Vec3 = dcel::Vec3T<T>;
  using mesh = dcel::meshT<T>;
  using face = dcel::faceT<T>;
  using AABB = dcel::AABBT<T>;
  using Node = dcel::BVH::NodeT<face, BV>;

  dcel_if() = delete;

  dcel_if(const std::shared_ptr<mesh>& a_mesh, const bool a_flipInside);
  dcel_if(const dcel_if& a_object);
  ~dcel_if();

  Real value(const RealVect& a_point) const override;
  
  BaseIF* newImplicitFunction() const override;

  inline
  decltype(auto) buildBVH(SortingFunctor<T> sortFunc = defaultSortingFunctor<T>, CostFunctor<T> costFunc = defaultCostFunctor<T>) noexcept;

protected:

  std::shared_ptr<mesh> m_mesh;

  std::shared_ptr<Node> m_root;
  
  bool m_flipInside;
};

#include "dcel_ifI.H"

#endif
