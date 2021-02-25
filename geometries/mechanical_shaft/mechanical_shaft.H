/*!
  @file   mechanical_shaft.H
  @brief  Declaration of the mechanical shaft geometry (ABB problem)
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _MECHANICAL_SHAFT_
#define _MECHANICAL_SHAFT_

#include "computational_geometry.H"

/*!
  @brief Mechanical shaft class
*/
class mechanical_shaft : public computational_geometry {
public:

  mechanical_shaft();
  ~mechanical_shaft();

protected:

  void define_electrode();
  void define_dielectric();
  
  RefCountedPtr<BaseIF> get_polygon();
  RefCountedPtr<BaseIF> get_cylinder();
  RefCountedPtr<BaseIF> get_cylinder_profile();
};

#endif
