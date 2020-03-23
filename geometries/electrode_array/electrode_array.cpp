/*!
  @file   electrode_array.cpp
  @brief  Implementation of electrode_array.H
  @author Robert Marskar
  @date   Dec. 2018
*/

#include "electrode_array.H"
#include "rod_if.H"
#include <ParmParse.H>

electrode_array::electrode_array(){
  if(SpaceDim != 3){
    MayDay::Abort("electrode_array::electrode_array - this is a pure 3D geometry");
  }
  
  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  ParmParse pp("electrode_array");

  int ilive;
  bool live;
  Vector<Real> a1, a2;
  Real r, dx, dy;
  RealVect c1, c2;
  int nx, ny;

  pp.get("live",    ilive);
  pp.get("radius",  r);
  pp.get("delta_x", dx);
  pp.get("delta_x", dy);
  pp.get("num_x",   nx);
  pp.get("num_y",   ny);
  pp.getarr("center1", a1, 0, SpaceDim);
  pp.getarr("center2", a2, 0, SpaceDim);

  live = (ilive == 0) ? false : true;
  c1 = RealVect(D_DECL(a1[0], a1[1], a1[2]));
  c2 = RealVect(D_DECL(a2[0], a2[1], a2[2]));


  for (int ix = 0; ix < nx; ix++){
    for (int iy = 0; iy < ny; iy++){
      const RealVect sx  = ix*dx*RealVect(BASISV(0));
      const RealVect sy  = iy*dy*RealVect(BASISV(1));
      
      const RealVect ic1 = c1 + sx + sy;
      const RealVect ic2 = c2 + sx + sy;

      RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF> (new rod_if(ic1, ic2, r, false));

      m_electrodes.push_back(electrode(rod, live));
    }
  }
}

electrode_array::~electrode_array(){

}
