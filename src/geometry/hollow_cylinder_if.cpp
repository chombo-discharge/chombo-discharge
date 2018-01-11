/*!
  @file hollow_cylinder_if.cpp
  @brief Implementation of hollow_cylinder_if.H
  @date Jan. 2017
  @author Sigurd Midttun
*/

#include "hollow_cylinder_if.H"
#include "hollow_cylinder_cross_sec_if.H"

hollow_cylinder_if::hollow_cylinder_if(const Real&     a_innerRadius,
				       const Real&     a_outerRadius,
				       const Real&     a_height,
				       const Real&     a_radiusOfCurvature,
				       const RealVect& a_center,
				       const bool&     a_inside,
				       const int&      a_dir){
  
  CH_assert(a_innerRadius > 0 && a_innerRadius < a_outerRadius);
  CH_assert(a_height > 0);
  CH_assert(a_radiusOfCurvature >= 0);
  CH_assert((a_outerRadius-a_innerRadius-2*a_radiusOfCurvature) >= 0);
  
  // Remember the parameters
  m_innerRadius        = a_innerRadius;
  m_outerRadius        = a_outerRadius;
  m_height             = a_height;
  m_radiusOfCurvature  = a_radiusOfCurvature;
  m_center             = a_center;
  m_inside             = a_inside;
  m_dir                = a_dir;

  //Build the cylinder and store it as m_object
  buildHollowCylinder();
}

hollow_cylinder_if::hollow_cylinder_if(const hollow_cylinder_if& a_inputIF){
  
  CH_assert(a_inputIF.m_innerRadius > 0 && a_inputIF.m_innerRadius < a_inputIF.m_outerRadius);
  CH_assert(a_inputIF.m_height > 0);
  CH_assert(a_inputIF.m_radiusOfCurvature >= 0 &&
	    (a_inputIF.m_outerRadius-a_inputIF.m_innerRadius-2*a_inputIF.m_radiusOfCurvature) >= 0);

  // Remember the parameters
  m_innerRadius        = a_inputIF.m_innerRadius;
  m_outerRadius        = a_inputIF.m_outerRadius;
  m_height             = a_inputIF.m_height;
  m_radiusOfCurvature  = a_inputIF.m_radiusOfCurvature;
  m_center             = a_inputIF.m_center;
  m_inside             = a_inputIF.m_inside;
  m_object             = a_inputIF.m_object;
  m_dir                = a_inputIF.m_dir;
}

hollow_cylinder_if::~hollow_cylinder_if(){
  delete m_object;
}

void hollow_cylinder_if::buildHollowCylinder(){
  const RealVect loCorner = RealVect(D_DECL(m_innerRadius, 0,        0));
  const RealVect hiCorner = RealVect(D_DECL(m_outerRadius, m_height, 0));
  
  //Create a hollow cylinder cross section object
  const BaseIF* hollowCylCrossSec = static_cast<BaseIF*> (new hollow_cylinder_cross_sec_if(hiCorner,  
											   loCorner,
											   m_radiusOfCurvature,
											   m_inside));
									         
  //LatheIF makes hollowCylCrossSec into a 3D object by revolving the cross-section around the y-axis in the xy-plane.
  //Afterwards, the z-axis is rotated to the position of the y-axis such that the z-axis points through
  //the center of the cylinder; not the y-axis.
  const LatheIF lathe(*hollowCylCrossSec, m_inside);
  


  //Rotate and translate the 3D object to its specified position. 
  m_object = new TransformIF(*(lathe.newImplicitFunction()));
  const Real pi = 3.14159265358979323846;
  if(m_dir == 0){
    m_object->rotate(0.5*pi, RealVect::Zero, BASISV(1));
  }
  else if(m_dir == 1){
    m_object->rotate(0.5*pi, RealVect::Zero, BASISV(0));
  }
  m_object->translate(m_center);

  //This is no longer needed
  delete hollowCylCrossSec;
}

Real hollow_cylinder_if::value(const RealVect& a_point) const{

  Real val = m_object->value(a_point);
  if(m_inside == 0){
    val = - val;
  }

  return val;
}

BaseIF* hollow_cylinder_if::newImplicitFunction() const{
  
  hollow_cylinder_if* hollowCylinderPtr = new hollow_cylinder_if(m_innerRadius,  
								 m_outerRadius,
								 m_height,
								 m_radiusOfCurvature,
								 m_center,
								 m_inside,
								 m_dir);

  return static_cast<BaseIF*> (hollowCylinderPtr);
}
