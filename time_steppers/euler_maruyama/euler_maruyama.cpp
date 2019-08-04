/*!
  @file   euler_maruyama.cpp
  @brief  Implementation of euler_maruyama.H
  @author Robert Marskar
  @date   Aug. 2019
*/

#include "euler_maruyama.H" 

euler_maruyama::euler_maruyama(){

}

euler_maruyama::~euler_maruyama(){

}

Real euler_maruyama::restrict_dt(){
  CH_TIME("euler_maruyama::euler_maruyama");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::euler_maruyama" << endl;
  }

  return 1.E99;
}

Real euler_maruyama::advance(const Real a_dt){
  if(m_verbosity > 5){
    pout() << "euler_maruyama::advance" << endl;
  }

  return a_dt;
}

void euler_maruyama::init_source_terms(){
  CH_TIME("euler_maruyama::init_source_terms");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::init_source_terms" << endl;
  }
}

void euler_maruyama::regrid_internals(){
  CH_TIME("euler_maruyama::regrid_internals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::regrid_internals" << endl;
  }
}

void euler_maruyama::deallocate_internals(){
  CH_TIME("euler_maruyama::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::deallocate_internals" << endl;
  }
}

bool euler_maruyama::need_to_regrid(){
  CH_TIME("euler_maruyama::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::need_to_regrid" << endl;
  }

  return false;
}
