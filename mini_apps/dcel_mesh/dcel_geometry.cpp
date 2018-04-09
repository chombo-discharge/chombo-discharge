#include "dcel_geometry.H"
#include "dcel_if.H"

dcel_geometry::dcel_geometry(){

}

dcel_geometry::dcel_geometry(const dcel_mesh* const a_mesh){

  m_dielectrics.resize(0);
  m_electrodes.resize(1);

  RefCountedPtr<dcel_mesh> mesh = RefCountedPtr<dcel_mesh> (const_cast<dcel_mesh*> (a_mesh));
  RefCountedPtr<BaseIF> electrode = RefCountedPtr<BaseIF> (new dcel_if(mesh, false)); 

  m_electrodes[0].define(electrode, true);
}

dcel_geometry::~dcel_geometry(){

}
