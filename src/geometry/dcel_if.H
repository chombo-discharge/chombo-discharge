/*!
  @file   dcel_if.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface tesselation
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_IF_
#define _DCEL_IF_

#include "dcel_mesh.H"

#include <RefCountedPtr.H>
#include <BaseIF.H>

/*!
  @brief Class that gets it's value function from a surface tesselation 
*/
class dcel_if : public BaseIF {
public:

  /*!
    @brief Weak constructor
  */
  dcel_if();

  /*!
    @brief Full constructor
  */
  dcel_if(const RefCountedPtr<dcel_mesh>& a_mesh, const bool a_inside);

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
  virtual Real value(const RealVect& a_point) const;

  /*!
    @brief Factory method
  */
  virtual BaseIF* newImplicitFunction() const;

protected:

  /*!
    @brief Mesh
  */
  RefCountedPtr<dcel_mesh> m_mesh;

  /*!
    @brief Flip covered region
  */
  bool m_inside;
};

#endif
