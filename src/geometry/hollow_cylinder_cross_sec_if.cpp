/*!
  @file hollow_cylinder_cross_sec_if.cpp
  @brief Implementation of hollow_cylinder_cross_sec_if.H
  @date Nov. 2017
  @author Robert Marskar
*/

#include "hollow_cylinder_cross_sec_if.H"

hollow_cylinder_cross_sec_if::hollow_cylinder_cross_sec_if(const RealVect& a_hiCorner,
							   const RealVect& a_loCorner,
							   const Real&     a_radOfCurv,
							   const bool&     a_inside){
 
  CH_assert(a_hiCorner[0] > a_loCorner[0] && a_hiCorner[1] > a_loCorner[1]);
  CH_assert((a_hiCorner[0]-a_loCorner[0] - 2*a_radOfCurv) >= 0);

  m_hiCorner        = a_hiCorner;
  m_loCorner        = a_loCorner;
  m_radOfCurv       = a_radOfCurv;
  m_inside          = a_inside;

}

hollow_cylinder_cross_sec_if::hollow_cylinder_cross_sec_if(const hollow_cylinder_cross_sec_if& a_inputIF){
  m_hiCorner           = a_inputIF.m_hiCorner;
  m_loCorner           = a_inputIF.m_loCorner;
  m_radOfCurv          = a_inputIF.m_radOfCurv;
  m_inside             = a_inputIF.m_inside;
}

hollow_cylinder_cross_sec_if::~hollow_cylinder_cross_sec_if(){
}

Real hollow_cylinder_cross_sec_if::value(const RealVect& a_point) const{
  
  //Declaration of some useful parameters describing the hollow cylindrical cross-section
  Real     retval;
  RealVect lowerRightCorner(D_DECL(m_hiCorner[0], m_loCorner[1], 0));  
  RealVect lowerLeftCorner  = m_loCorner;  
  RealVect upperRightCorner = m_hiCorner;
  RealVect centerRight(D_DECL(lowerRightCorner[0]-m_radOfCurv, lowerRightCorner[1]+m_radOfCurv, 0));  
  RealVect centerLeft(D_DECL(lowerLeftCorner[0]+m_radOfCurv, lowerLeftCorner[1]+m_radOfCurv, 0));

  
  //1. Determine if a_point is inside of the rectangular cross section
  if (a_point[0] >= lowerLeftCorner[0] && a_point[1] >= lowerLeftCorner[1] && a_point[0] <= upperRightCorner[0] && a_point[1] <= upperRightCorner[1]){
    //a_point is inside of the rectangle
    retval = -1;

    
    //2. Determine if a_point is on the inside of the rounded corners.
    //Lower left corner:
    if (a_point[0] >= lowerLeftCorner[0] && a_point[1] >= lowerLeftCorner[1] && a_point[0] <= centerLeft[0] && a_point[1] <= centerLeft[1]){

      //check if a_point is inside of the circular sector (curved corner)

      // The distance squared for center to a_point
      Real distance2;

      // Compute the distance squared
      distance2 = 0.0;
      for (int idir = 0; idir < 2; idir++)
	{
	  Real cur;
	  cur = a_point[idir] - centerLeft[idir];

	  distance2 += cur*cur;
	}

      // Return the difference between the sqaures (zero on the sphere)
      retval = distance2 - m_radOfCurv*m_radOfCurv;
    }
      

    
    //Lower right corner:
    if (a_point[0] <= lowerRightCorner[0] && a_point[0] >= centerRight[0] && a_point[1] >= lowerRightCorner[1] && a_point[1] <= centerRight[1]){

      //check if point is inside of the circular sector (curved corner)

      // The distance squared for center to a_point
      Real distance2;

      // Compute the distance squared
      distance2 = 0.0;
      for (int idir = 0; idir < 2; idir++)
	{
	  Real cur;
	  cur = a_point[idir] - centerRight[idir];

	  distance2 += cur*cur;
	}

      // Return the difference between the sqaures (zero on the sphere)
      retval = distance2 - m_radOfCurv*m_radOfCurv;
    } 
  }//rectangle
   
  //The point is not inside of the rectangle
  else{
    retval = 1;
  }


  if (!m_inside){
    retval = -retval;
  }
  return retval;
  
}

//
BaseIF* hollow_cylinder_cross_sec_if::newImplicitFunction() const{
  hollow_cylinder_cross_sec_if* fresh = new hollow_cylinder_cross_sec_if(m_hiCorner,  
									 m_loCorner,
									 m_radOfCurv,
									 m_inside);

  return static_cast<BaseIF*> (fresh);
}

