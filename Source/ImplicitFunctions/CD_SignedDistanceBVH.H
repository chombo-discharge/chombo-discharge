/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SignedDistanceBVH.H
  @brief  Declaration of an signe distance function class that gets its value function from a DCEL surface Tesselation
  @author Robert Marskar
*/

#ifndef CD_SignedDistanceBVH_H
#define CD_SignedDistanceBVH_H

// Std includes
#include <memory>
#include <chrono>

// EBGeometry includes
#include <EBGeometry.hpp>

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Signed distance function for a DCEL mesh. 
  @details This class uses the EBGeometry submodule for creating signed distance functions from compelx geometries. The template parameters 
  are the precision and the bounding volume type. To use this class, the user must first create the DCEL mesh and then create the BVH. 
*/
template <class T, class BV, int K>
class SignedDistanceBVH : public BaseIF
{
public:
  /*!
    @brief Alias for always-3D vector with template
  */
  using Vec3 = EBGeometry::Vec3T<T>;

  /*!
    @brief DCEL polygon face, which is the primitive type enclosed in the bounding volume hierarchy
  */
  using Face = EBGeometry::DCEL::FaceT<T>;

  /*!
    @brief Alias for BVH node with precision T, enclosing DCEL face with a bounding volume type V
  */
  using Node = EBGeometry::BVH::LinearBVH<T, Face, BV, K>;

  /*!
    @brief Disallowed weak construction
  */
  SignedDistanceBVH() = delete;

  /*!
    @brief Full constructor
    @param[in] a_root       Pointer to root node in the bounding volume hierarchy. 
    @param[in] a_flipInside If true, turns inside to outside
    @param[in] a_zCoord     Special flag for 2D -- indicates the z-position through which we slice the 3D object. Will no longer be signed distance. 
  */
  SignedDistanceBVH(const std::shared_ptr<Node>& a_root, const bool a_flipInside, const Real a_zCoord = 0.0);

  /*!
    @brief Copy constructor
  */
  SignedDistanceBVH(const SignedDistanceBVH& a_object);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~SignedDistanceBVH();

  /*!
    @brief Value function
    @param[in] a_point Spatial point
  */
  Real
  value(const RealVect& a_point) const override;

  /*!
    @brief Factory method. Sends pointers around. 
  */
  BaseIF*
  newImplicitFunction() const override;

protected:
  /*!
    @brief For 2D only. This is the z-coordinate through which we slice the object
  */
  Real m_zCoord;

  /*!
    @brief Pointer to root node in bounding volume hierarchy
  */
  std::shared_ptr<Node> m_root;

  /*!
    @brief Hook for turning inside to outside
  */
  bool m_flipInside;

  /*!
    @brief Number of times the implicit function was called. Use for performance tracking. 
  */
  mutable long m_numCalled;

  /*!
    @brief Total time spent computing the signed distance. Use for performance tracking. 
  */
  mutable std::chrono::duration<double> m_timespan;
};

#include <CD_NamespaceFooter.H>

#include <CD_SignedDistanceBVHImplem.H>

#endif
