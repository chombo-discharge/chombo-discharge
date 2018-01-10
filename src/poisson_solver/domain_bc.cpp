/*!
  @file domain_bc.cpp
  @brief Impementation of domain_bc.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "domain_bc.H"

domain_bc::domain_bc(Vector<RefCountedPtr<wall_bc> >& a_bc){
  CH_TIME("domain_bc::domain_bc");
  
  this->define(a_bc);

  m_values.resize(2*SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const int which = wall_bc::map_bc(dir, sit());
      m_values[which] = m_bc[which]->get_value();
    }
  }

  
}

domain_bc::~domain_bc(){
  CH_TIME("domain_bc::~domain_bc");
}

void domain_bc::set_potential_ptr(Real (*a_ptr)(const Real a_time)){
  CH_TIME("domain_bc::set_potential_ptr");
  m_ptr = a_ptr;
}

void domain_bc::define(Vector<RefCountedPtr<wall_bc> >& a_bc){
  CH_TIME("domain_bc::define");
  for (int i = 0; i < a_bc.size(); i++){
    CH_assert(!a_bc[i].isNull());
  }

  m_bc = a_bc;
}

Real domain_bc::value(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const{
  CH_TIME("domain_bc::value(1)");

  MayDay::Abort("domain_bc::value - how did you get here, this shouldn't be called anywhere!");
}


Real domain_bc::value(const FaceIndex&      a_face, 
		      const Side::LoHiSide& a_side, 
		      const DataIndex&      a_dit, 
		      const RealVect&       a_point,
		      const RealVect&       a_normal, 
		      const Real&           a_time, 
		      const int&            a_comp) const {
  CH_TIME("domain_bc::value(2)");

  const int dir   = a_face.direction();
  const int which = wall_bc::map_bc(dir, a_side);

  return m_values[which];
}
