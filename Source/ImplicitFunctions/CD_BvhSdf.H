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

/*!
  @brief Signed distance function for a DCEL mesh. 
  @details Template parameters are the precision and the bounding volume type. To use this class, the user must first create the DCEL mesh and then create
  the bounding volume hierarchy. 
*/
template <class T, class BV, int K>
class BvhSdf : public BaseIF {
public:

  /*!
    @brief Alias for always-3D vector with template
  */
  using Vec3  = Vec3T<T>;

  /*!
    @brief DCEL polygon face, which is the primitive type enclosed in the bounding volume hierarchy
  */
  using Face  = FaceT<T>;

  /*!
    @brief Alias for BVH node with precision T, enclosing DCEL face with a bounding volume type V
  */
  using Node  = NodeT<T, Face, BV, K>;

  /*!
    @brief Disallowed weak construction
  */
  BvhSdf() = delete;

  /*!
    @brief Full constructor
    @param[in] a_root       Pointer to root node in bounding volume hierarchy
    @param[in] a_flipInside If true, turns inside to outside
    @param[in] a_zCoord     Special flag for 2D -- indicates the z-position through which we slice the 3D object. 
  */
  BvhSdf(const std::shared_ptr<Node>& a_root, const bool a_flipInside, const Real a_zCoord = 0.0);

  /*!
    @brief Copy constructor
  */
  BvhSdf(const BvhSdf& a_object);

  /*!
    @brief Destructor (does nothing)
  */
  ~BvhSdf();

  /*!
    @brief Value function
    @param[in] a_point Spatial point
  */
  Real value(const RealVect& a_point) const override;

  /*!
    @brief Factory method. Sends pointers around. 
  */
  BaseIF* newImplicitFunction() const override;

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

#include <CD_BvhSdfImplem.H>

#endif
