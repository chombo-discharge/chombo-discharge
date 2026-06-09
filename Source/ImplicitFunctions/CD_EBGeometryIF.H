/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBGeometryIF.H
  @brief  Binding of EBGeometry signed distance functions to Chombo objcts. 
  @author Robert marskar
*/

#ifndef CD_EBGEOMETRYIF_H
#define CD_EBGEOMETRYIF_H

// Std includes
#include <memory>

// Chombo includes
#include <BaseIF.H>
#include <IntersectionIF.H>

// EBGeometry includes
#include <EBGeometry.hpp>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Shallow binding of EBGeometry signed distance functions to Chombo. 
  @details Since EBGeometry is inherently 3D and chombo-discharge is both 2D and 3D, the object will be sliced through a specified z-coordinate doing 2D simulations. 
  @note T is the precision used in EBGeometry.
*/
template <typename T = Real>
class EBGeometryIF : public BaseIF
{
public:
  /**
    @brief Constructor. Sets shape to nullptr. 
  */
  EBGeometryIF();

  /**
    @brief Full constructor
    @param[in] a_flipInside Hook for making inside -> outside.
    @param[in] a_zCoord     z-coordinate through which the object is sliced. 
    @param[in] a_sdf Sdf
  */
  EBGeometryIF(const std::shared_ptr<EBGeometry::ImplicitFunction<T>>& a_sdf,
               const bool                                              a_flipInside,
               const Real                                              a_zCoord = 0.0);

  /**
    @brief Copy constructor
    @param[in] a_inputIF    Other implicit function. 
  */
  EBGeometryIF(const EBGeometryIF& a_inputIF);

  /**
    @brief Destructor. Does nothing
  */
  virtual ~EBGeometryIF();

  /**
    @brief Value function. Returns distance to object. 
    @param[in] a_point Physical position. 
    @return Function value at the given point
  */
  virtual Real
  value(const RealVect& a_point) const override;

  /**
    @brief IF factory method
    @return New copy of this implicit function
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  /**
    @brief Signed distance function
  */
  std::shared_ptr<const EBGeometry::ImplicitFunction<T>> m_sdf;

  /**
    @brief Hook for making outside -> inside
  */
  bool m_flipInside;

  /**
    @brief z-coordinate through which the object is sliced.
  */
  Real m_zCoord;
};

#include <CD_NamespaceFooter.H>

#include <CD_EBGeometryIFImplem.H>

#endif
