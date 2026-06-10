/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RealBox.H
  @brief  Declaration of an axis-aligned box in space
  @author Robert Marskar
*/

#ifndef CD_REALBOX_H
#define CD_REALBOX_H

// Chombo includes
#include <RealVect.H>
#include <Box.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Rectangular box in real space. 
*/
class RealBox
{
public:
  /**
    @brief Default constructor. Must subsequently set the corners. 
  */
  RealBox();

  /**
    @brief Full constructor. Sets corners from input arguments
    @param[in] a_lo Low corner
    @param[in] a_hi High corners
  */
  RealBox(const RealVect& a_lo, const RealVect& a_hi);

  /**
    @brief Full constructor. Sets corners from input box
    @param[in] a_box    IntVect-box
    @param[in] a_probLo Lower-left corner of physical domain
    @param[in] a_dx     Grid resolution.
    @details This sets the low corner as a_probLo + a_box.lowEnd()*a_dx
  */
  RealBox(const Box a_box, const RealVect& a_origin, const Real a_dx);

  /**
    @brief Destructor (does nothing)
  */
  ~RealBox() = default;

  /**
    @brief Get low corner
    @return Lo
  */
  RealVect
  getLo() const;

  /**
    @brief Get high corner
    @return Hi
  */
  RealVect
  getHi() const;

  /**
    @brief Get all the corners of this box (4 in 2D, 8 in 3D);
    @return Corners
  */
  Vector<RealVect>
  getCorners() const;

  /**
    @brief Check if a_box intersects this box. This is done by checking if 
    @param[in] a_box Other box
    @return True if intersection found, false otherwise
  */
  bool
  intersect(const RealBox& a_box) const;

  /**
    @brief Check if a_point is inside this box
    @param[in] a_point Input point
    @return True if point inside, false otherwise
  */
  bool
  isPointInside(const RealVect& a_point) const;

  /**
    @brief Check if a_box is completely inside this box. For this to be true, all the corners of a_box must be 
    completely inside this box
    @param[in] a_box Other box.
    @return True if box inside, false otherwise
  */
  bool
  isBoxInside(const RealBox& a_box) const;

protected:
  /**
    @brief Low corner
  */
  RealVect m_lo;

  /**
    @brief High corner.
  */
  RealVect m_hi;
};

#include <CD_NamespaceFooter.H>

#endif
