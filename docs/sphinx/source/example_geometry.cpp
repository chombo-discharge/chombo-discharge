#include "computational_geometry.H"
#include "SphereIF.H" // Chombo class

class my_geometry : public computational_geometry {
public:

  my_geometry(){
    m_electrodes.resize(1);
    m_dielectrics.resize(0);

    const bool live       = true;
    const Real radius     = 1;
    const RealVect center = RealVect::Zero;
    
    RefCountedPtr<BaseIF> my_sphere = RefCountedPtr<BaseIF> (new SphereIF(radius, center, false));

    m_electrodes[0].define(my_sphere, live);
  }

  ~my_geometry(){}
};
