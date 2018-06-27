/*!
  @file   conductivitydomainbc_wrapper_factory.cpp
  @brief  Implementation of conductivitydomainbc_wrapper_factory
  @author Robert Marskar
  @date   June 2018
*/


#include "conductivitydomainbc_wrapper_factory.H"

conductivitydomainbc_wrapper_factory::conductivitydomainbc_wrapper_factory(){
  m_hasbc = false;
}

conductivitydomainbc_wrapper_factory::~conductivitydomainbc_wrapper_factory(){

}

void conductivitydomainbc_wrapper_factory::set_wallbc(const Vector<RefCountedPtr<wall_bc> >& a_wallbc){
  m_wallbc = a_wallbc;
  
  m_hasbc = true;
}

void conductivitydomainbc_wrapper_factory::set_potentials(const Vector<RefCountedPtr<BaseBCFuncEval> >& a_potentials){
  m_potentials = a_potentials;
}

conductivitydomainbc_wrapper* conductivitydomainbc_wrapper_factory::create(const ProblemDomain& a_domain,
									   const EBISLayout&    a_ebisl,
									   const RealVect&      a_dx){
  CH_assert(m_hasbc);
    
  conductivitydomainbc_wrapper* fresh = new conductivitydomainbc_wrapper();

  fresh->set_potentials(m_potentials);
  fresh->set_wallbc(m_wallbc);
    
  return fresh;
}
