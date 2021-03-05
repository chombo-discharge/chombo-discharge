/*!
  @file   dcel_if.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface tesselation
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_IF_
#define _DCEL_IF_

#include "dcel_mesh.H"

#include <memory>

#include <BaseIF.H>

/*!
  @brief Class that gets it's value function from a surface tesselation 
*/
class dcel_if : public BaseIF {
public:

  dcel_if() = delete;

  dcel_if(const std::shared_ptr<dcel::mesh>& a_mesh, const bool a_fluidInside);
  dcel_if(const dcel_if& a_object);
  ~dcel_if();

  Real value(const RealVect& a_point) const override;
  BaseIF* newImplicitFunction() const override;

protected:

  std::shared_ptr<dcel::mesh> m_mesh;
  bool m_fluidInside;
};

#endif
