/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelSdf.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface tesselation
  @author Robert Marskar
  @date   March 2021
*/

#ifndef CD_DcelSdf_H
#define CD_DcelSdf_H

// Std includes
#include <memory>

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_DcelMesh.H>
#include <CD_NamespaceHeader.H>

using namespace Dcel;

template <class T>
class DcelSdf : public BaseIF {
public:

  using Vec3 = Vec3T<T>;
  using Mesh = MeshT<T>;

  DcelSdf() = delete;
  DcelSdf(const std::shared_ptr<Mesh>& a_mesh, const bool a_flipInside);
  DcelSdf(const DcelSdf& a_object);
  ~DcelSdf();

  Real value(const RealVect& a_point) const override;
  
  BaseIF* newImplicitFunction() const override;

protected:

  std::shared_ptr<Mesh> m_mesh;
  
  bool m_flipInside;
};

#include <CD_NamespaceFooter.H>

#include <CD_DcelSdfImplem.H>

#endif
