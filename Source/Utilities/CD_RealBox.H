/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RealBox.H
  @brief  Declaration of an axis-aligned box in space
  @author Robert Marskar
*/

#ifndef CD_RealBox_H
#define CD_RealBox_H

// Chombo includes
#include <RealVect.H>
#include <Box.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Rectangular box
*/
class RealBox{
public:

  RealBox();
  RealBox(const RealVect a_lo, const RealVect a_hi);
  RealBox(const Box a_box, const RealVect a_origin, const Real a_dx);
  
  ~RealBox();

  /*!
    @brief Get low corner
  */
  RealVect getLo() const;

  /*!
    @brief Get high corner
  */
  RealVect getHi() const;

  /*!
    @brief Get all the corners of this box (4 in 2D, 8 in 3D);
  */
  Vector<RealVect> getCorners() const;

  /*!
    @brief Check if a_box intersects this box. This is done by checking if 
  */
  bool intersect(const RealBox& a_box) const;

  /*!
    @brief Check if a_point is inside this box
  */
  bool isPointInside(const RealVect a_point) const;

  /*!
    @brief Check if a_box is completely inside this box. For this to be true, all the corners of a_box must be 
    completely inside this box
  */
  bool isBoxInside(const RealBox& a_box) const;
  
protected:

  RealVect m_lo;
  RealVect m_hi;

};

#include <CD_NamespaceFooter.H>

#endif
