/*!
  @file    mushroom_if.H
  @brief   Declaration of a "mushroom", used in HV vessel
  @author  Robert Marskar
  @date    Nov. 2017
*/

#ifndef _MUSHROOM_IF_
#define _MUSHROOM_IF_

#include <BaseIF.H>
#include <RealVect.H>

class mushroom_if : public BaseIF {
public:  

  mushroom_if(const RealVect a_center,
	      const Real     a_R,
	      const Real     a_r,
	      const Real     a_L,
	      const Real     a_d,
	      const Real     a_curv,
	      const bool     a_inside);

  mushroom_if(const mushroom_if& a_inputIF);
  virtual ~mushroom_if();
  virtual Real value(const RealVect& a_pos) const;
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  RefCountedPtr<BaseIF> m_baseif;

    
};
#endif
