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

namespace ChomboDischarge {

  class mushroom_if : public BaseIF {
  public:  

    mushroom_if(const RealVect a_center,
		const Real     a_R,
		const Real     a_r,
		const Real     a_L,
		const Real     a_d,
		const Real     a_curv,
		const bool     a_fluidInside);

    mushroom_if(const mushroom_if& a_inputIF);
    virtual ~mushroom_if();
    virtual Real value(const RealVect& a_pos) const override;
    virtual BaseIF* newImplicitFunction() const override;
  
  protected:

    bool m_fluidInside;
    RefCountedPtr<BaseIF> m_baseif;

    
  };
}
#endif
