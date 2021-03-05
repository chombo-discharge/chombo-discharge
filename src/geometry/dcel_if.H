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

  /*!
    @brief Full constructor
  */
  dcel_if(const std::shared_ptr<dcel::dcel_mesh>& a_mesh, const bool a_fluidInside);

  /*!
    @brief Copy constructor
  */
  dcel_if(const dcel_if& a_object);

  /*!
    @brief Destructor
  */
  ~dcel_if();

  /*!
    @brief Value function
  */
  Real value(const RealVect& a_point) const override;

  /*!
    @brief Factory method
  */
  BaseIF* newImplicitFunction() const override;

protected:

  /*!
    @brief Mesh
  */
  std::shared_ptr<dcel::dcel_mesh> m_mesh;

  /*!
    @brief Flip covered region
  */
  bool m_fluidInside;
};

#endif
