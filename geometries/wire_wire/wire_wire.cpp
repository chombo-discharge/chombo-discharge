/*!
  @file   wire_wire.cpp
  @brief  Implementation of wire_wire.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "wire_wire.H"

#include "new_sphere_if.H"
#include "cylinder_if.H"

wire_wire::wire_wire(){
  ParmParse pp1("wire_wire.first");
  ParmParse pp2("wire_wire.second");

  bool useFirst;
  bool useSecond;

  pp1.get("on", useFirst);
  pp2.get("on", useSecond);

  if(useFirst)  this->addWire(pp1);
  if(useSecond) this->addWire(pp2);
}


wire_wire::~wire_wire(){
}

void wire_wire::addWire(ParmParse& a_pp){
    Real r;
    RealVect e1, e2;
    bool live;

    Vector<Real> v;

    a_pp.get("radius", r);
    a_pp.get("live", live);

    a_pp.getarr("endpoint1", v, 0, SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    a_pp.getarr("endpoint2", v, 0, SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));

    Real pot = live ? 0.5 : -0.5;

#if CH_SPACEDIM==2
    RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF> (new new_sphere_if(e1, r, false));
#elif CH_SPACEDIM==3
    RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF> (new cylinder_if(e1, e2, r, false));
#endif

    m_electrodes.push_back(electrode(bif, true, pot));
}
