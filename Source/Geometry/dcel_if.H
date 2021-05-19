/*!
  @file   dcel_if.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface tesselation
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_IF_
#define _DCEL_IF_

#include "dcel_mesh.H"

#include <memory>

#include <BaseIF.H>

#include "CD_NamespaceHeader.H"

using namespace Dcel;

template <class T>
class dcel_if : public BaseIF {
public:

  using Vec3 = Vec3T<T>;
  using mesh = meshT<T>;

  dcel_if() = delete;
  dcel_if(const std::shared_ptr<mesh>& a_mesh, const bool a_flipInside);
  dcel_if(const dcel_if& a_object);
  ~dcel_if();

  Real value(const RealVect& a_point) const override;
  
  BaseIF* newImplicitFunction() const override;

protected:

  std::shared_ptr<mesh> m_mesh;
  
  bool m_flipInside;
};
#include "CD_NamespaceFooter.H"

#include "dcel_ifI.H"

#endif
