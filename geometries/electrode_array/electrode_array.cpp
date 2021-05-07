/*!
  @file   electrode_array.cpp
  @brief  Implementation of electrode_array.H
  @author Robert Marskar
  @date   Dec. 2018
*/

#include "electrode_array.H"
#include "rod_if.H"
#include <ParmParse.H>

namespace ChomboDischarge {
  
  electrode_array::electrode_array(){
    ParmParse pp("electrode_array");

    bool live;
    Vector<Real> v(SpaceDim);
    Real r, dx, dy;
    RealVect c1, c2;
    int nx, ny;

    pp.get("live",    live);
    pp.get("radius",  r);
    pp.get("delta_x", dx);
    pp.get("delta_x", dy);
    pp.get("num_x",   nx);
    pp.get("num_y",   ny);
  
    pp.getarr("endpoint1", v, 0, SpaceDim);  c1 = RealVect(D_DECL(v[0], v[1], v[2]));
    pp.getarr("endpoint2", v, 0, SpaceDim);  c2 = RealVect(D_DECL(v[0], v[1], v[2]));

#if CH_SPACEDIM==2 // Override in 2D. 
    ny = 1;
    dy = 0.0;
#endif
    for (int ix = 0; ix < nx; ix++){
      for (int iy = 0; iy < ny; iy++){
	const RealVect sx  = ix*dx*BASISREALV(0);
	const RealVect sy  = iy*dy*BASISREALV(1);
      
	const RealVect ic1 = c1 + sx + sy;
	const RealVect ic2 = c2 + sx + sy;

	RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF> (new rod_if(ic1, ic2, r, false));

	m_electrodes.push_back(electrode(rod, live));
      }
    }
  }

  electrode_array::~electrode_array(){

  }
}
