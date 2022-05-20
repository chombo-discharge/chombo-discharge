/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SignedDistanceDCEL.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface Tesselation
  @author Robert Marskar
  @date   March 2021
*/

#ifndef CD_SignedDistanceDCEL_H
#define CD_SignedDistanceDCEL_H

// Std includes
#include <memory>

// EBGeometry includes
#include <EBGeometry.hpp>

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Signed distance function from a DCEL mesh. 
  @details This does not use a bounding volume hierarchy for performance. The template parameter is the floating point precision. 
*/
template <class T>
class SignedDistanceDCEL : public BaseIF
{
public:
  /*!
    @brief Alias for EBGeometry DCEl mesh with precision T
  */
  using Mesh = EBGeometry::DCEL::MeshT<T>;

  /*!
    @brief Disallowed, use the full constructor
  */
  SignedDistanceDCEL() = delete;

  /*!
    @brief Full constructor
    @param[in] a_mesh DCEL mesh
    @param[in] a_flipInside Hook for turning inside to outside
    @param[in] a_zCoord 2D only. This is the z-coordinate through which we slice the 3D object. 
  */
  SignedDistanceDCEL(const std::shared_ptr<Mesh>& a_mesh, const bool a_flipInside, const Real a_zCoord = 0.0);

  /*!
    @brief Copy constructor
    @param[in] a_object Other distance function
  */
  SignedDistanceDCEL(const SignedDistanceDCEL& a_object);

  /*!
    @brief Destructor (does nothing)
  */
  ~SignedDistanceDCEL();

  /*!
    @brief Value function
    @param[in] a_point 2D/3D point
  */
  Real
  value(const RealVect& a_point) const override;

  /*!
    @brief Factory method
  */
  BaseIF*
  newImplicitFunction() const override;

protected:
  /*!
    @brief DCEL mesh
  */
  std::shared_ptr<Mesh> m_mesh;

  /*!
    @brief Hook for turning inside to outside
  */
  bool m_flipInside;

  /*!
    @brief For 2D only. z-coordinate through which we slice the 3D object
  */
  Real m_zCoord;
};

#include <CD_NamespaceFooter.H>

#include <CD_SignedDistanceDCELImplem.H>

#endif
