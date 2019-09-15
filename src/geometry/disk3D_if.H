/*!
  @file disk3D_if.H
*/

#ifndef _DISK3D_IF_
#define _DISK3D_IF_

#include <BaseIF.H>
#include <RealVect.H>

class disk3D_if : public BaseIF {
public:

  disk3D_if(const RealVect a_point, const Real a_R, const Real a_d, const Real a_curv, const bool a_inside);
  virtual ~disk3D_if();

  disk3D_if(const disk3D_if& a_inputIF);
  virtual Real value(const RealVect& a_pos) const;
  virtual BaseIF* newImplicitFunction() const;

protected:

  RefCountedPtr<BaseIF> m_baseif;
};

#endif
