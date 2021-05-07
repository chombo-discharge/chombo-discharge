/*!
  @file coaxial_cable.cpp
  @brief Implementation of coaxial_cable.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "coaxial_cable.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>

#include "cylinder_if.H"
#include "new_sphere_if.H"

namespace ChomboDischarge {

  coaxial_cable::coaxial_cable(){

    ParmParse pp      ("coaxial_cable");
    ParmParse ppOuter ("coaxial_cable.outer");
    ParmParse ppInner ("coaxial_cable.inner");
    ParmParse ppMiddle("coaxial_cable.dielectric");

    // Get centers
    Vector<Real> v(SpaceDim);
    RealVect e1, e2;
    pp.getarr("endpoint1", v, 0, SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    pp.getarr("endpoint2", v, 0, SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));
  
    // Check which ones to use. 
    bool outer, inner, middle;
    ppOuter.get ("on", outer);
    ppInner.get ("on", inner);
    ppMiddle.get("on", middle);

    // Reset
    m_electrodes.resize(0);
    m_dielectrics.resize(0);
    this->set_eps0(1.0);

    if(outer){ // Add outer electrode
      Real rad;
      bool live;

      ppOuter.get("radius", rad);
      ppOuter.get("live",   live);

#if CH_SPACEDIM==2
      RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new new_sphere_if(e1, rad, true));
#elif CH_SPACEDIM==3
      RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new cylinder_if(e1, e2, rad, true));
#endif
      m_electrodes.push_back(electrode(baseif, live));
    }

    if(inner){ // Add inner electrode
      Real rad;
      bool live;

      ppInner.get("radius", rad);
      ppInner.get("live",   live);

#if CH_SPACEDIM==2
      RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new new_sphere_if(e1, rad, false));
#elif CH_SPACEDIM==3
      RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new cylinder_if(e1, e2, rad, false));
#endif
      m_electrodes.push_back(electrode(baseif, live));
    }

    if(middle){ // Add dielectric
      Real rad, eps;

      ppMiddle.get("radius", rad);
      ppMiddle.get("eps",    eps);


#if CH_SPACEDIM==2
      RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new new_sphere_if(e1, rad, false));
#elif CH_SPACEDIM==3
      RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new cylinder_if(e1, e2, rad, false));
#endif
      m_dielectrics.push_back(dielectric(baseif, eps));
    }
  }

  coaxial_cable::~coaxial_cable(){
  
  }
}
