/*!
  @file   real_box.H
  @brief  Declaration of a rectangular box
  @author Robert Marskar
  @date   June 2018
*/

#ifndef _REAL_BOX_
#define _REAL_BOX_

#include <RealVect.H>
#include <Box.H>

/*!
  @brief Rectangular box
*/
class real_box{
public:

  real_box();
  real_box(const RealVect a_lo, const RealVect a_hi);
  real_box(const Box a_box, const RealVect a_origin, const Real a_dx);
  
  ~real_box();

  /*!
    @brief Get low corner
  */
  RealVect get_lo() const;

  /*!
    @brief Get high corner
  */
  RealVect get_hi() const;

  /*!
    @brief Get all the corners of this box (4 in 2D, 8 in 3D);
  */
  Vector<RealVect> get_corners() const;

  /*!
    @brief Check if a_box intersects this box. This is done by checking if 
  */
  bool intersect(const real_box& a_box) const;

  /*!
    @brief Check if a_point is inside this box
  */
  bool is_point_inside(const RealVect a_point) const;

  /*!
    @brief Check if a_box is completely inside this box. For this to be true, all the corners of a_box must be 
    completely inside this box
  */
  bool is_box_inside(const real_box& a_box) const;
  
protected:

  RealVect m_lo;
  RealVect m_hi;

};

#endif
