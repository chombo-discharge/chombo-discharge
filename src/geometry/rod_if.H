/*!
  @file rod_if.H
  @brief Declaration of a cylinder with spheres at both ends
  @date Nov. 2017
  @author Robert marskar
*/

#ifndef _ROD_IF_
#define _ROD_IF_

#include <BaseIF.H>
#include <IntersectionIF.H>

#include <ParmParse.H>

/*!
  @brief Declaration of a cylinder with rounded caps
*/
class rod_if : public BaseIF{
public:

  /*!
    @brief Constructor
    @param[in] a_center1 Cylinder center (at center of sphere)
    @param[in] a_center2 Cylinder center (at center of sphere)
    @param[in] a_radius  Cylinder and sphere radius
    @param[in] a_inside Domain on inside or outside
  */

  rod_if(const RealVect& a_center1, const RealVect& a_center2, const Real& a_radius, const bool& a_inside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */
  
  rod_if(const rod_if& a_inputIF);

  /*!
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
  */
  virtual Real value(const RealVect& a_point) const; 

  /*!
    @brief IF factory method
  */
  virtual BaseIF* newImplicitFunction() const;

  /*!
    @brief Fast intersection
  */
  virtual bool fastIntersection(const RealVect& a_low, const RealVect& a_high) const{
    ParmParse pp("debug");
    int val = 0;
    pp.query("fast_intersection", val);

    if(val == 0){
      return false;
    }
    else{
      return true;
    }
  }

  /*!
    @brief Check if point is inside or outside
  */
  virtual GeometryService::InOut InsideOutside(const RealVect& a_low, const RealVect& a_high) const {

    RealVect lo = m_center2 - m_radius*RealVect(1.,1, 0);
    RealVect hi = m_center1 + m_radius*RealVect(1.,1, 0);

    std::cout << lo << std::endl;

    Vector<RealVect> corners(8);
    const RealVect DD = a_high - a_low;
    const RealVect DX = RealVect(DD[0], 0.0,   0.0);
    const RealVect DY = RealVect(0.0,   DD[1], 0.0);
    const RealVect DZ = RealVect(0.0,   0.0,   DD[2]);
    corners[0] = a_low;
    corners[1] = a_low + DX;
    corners[2] = a_low + DY;
    corners[3] = a_low + DZ;
    corners[4] = a_low + DX + DY;
    corners[5] = a_low + DX + DZ;
    corners[6] = a_low + DY + DZ;
    corners[7] = a_low + DX + DY + DZ;

    bool outside = true;
    for (int icorner = 0; icorner < corners.size(); icorner++){
      RealVect point = corners[icorner];

      if(point[0] >= lo[0] && point[0] <= hi[0] &&
	 point[1] >= lo[1] && point[1] <= hi[1] &&
	 point[2] >= lo[2] && point[2] <= hi[2]){
	outside = false;
	break;
      }
    }

    if(outside){
      return GeometryService::Regular;
    }
    else{
      pout() << "irr box = " << "\t" << corners[0] << "\t" << corners[7] << endl;
      return GeometryService::Irregular;
    }
  }

protected:

  /*!
    @brief Actual implicit function
  */
  RefCountedPtr<BaseIF> m_baseif;

  RealVect m_center1;
  RealVect m_center2;
  Real     m_radius;
  bool     m_inside;

  
};
#endif
