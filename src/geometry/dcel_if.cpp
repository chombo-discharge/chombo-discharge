/*!
  @file   dcel_if.cpp
  @brief  Implementation of dcel_if.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_if.H"

#include <chrono>

#define enable_dcel_timer 0

dcel_if::dcel_if(const std::shared_ptr<dcel_mesh>& a_mesh, const bool a_fluidInside){
  m_mesh        = a_mesh;
  m_fluidInside = a_fluidInside;
}

dcel_if::dcel_if(const dcel_if& a_object){
  m_mesh   = a_object.m_mesh;
  m_fluidInside = a_object.m_fluidInside;
}

dcel_if::~dcel_if(){

}

Real dcel_if::value(const RealVect& a_point) const{
#if dcel_timer
  auto mesh_start = std::chrono::system_clock::now(); 
#endif
  Real retval = m_mesh->signed_distance(a_point);
#if dcel_timer
  auto mesh_stop = std::chrono::system_clock::now();
  const Real blah = a_point.vectorLength() - 1.0;
  auto imp_stop = std::chrono::system_clock::now();
  std::chrono::duration<double> decl_time = mesh_stop - mesh_start;
  std::chrono::duration<double> impfunc_time = imp_stop-mesh_stop;
  pout() << "Ratio = " << 1.0*decl_time.count()/impfunc_time.count() << endl;
#endif
  if(m_fluidInside){
    retval *= -1;
  }

  return retval;
}

BaseIF* dcel_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new dcel_if(*this));
}
