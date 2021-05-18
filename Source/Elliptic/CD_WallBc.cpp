/*!
  @file WallBc.cpp
  @brief Implementation of WallBc.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include <CD_WallBc.H>

#include "CD_NamespaceHeader.H"

WallBc::WallBc(const int a_dir, const Side::LoHiSide a_side, wallbc::which_bc a_which){
  m_dir   = a_dir;
  m_side  = a_side;
  m_which = a_which;
}

WallBc::~WallBc(){
}

void WallBc::set_value(Real a_value){
  m_value = a_value;
}

void WallBc::set_live(bool a_live){
  m_live = a_live;
}

void WallBc::set_function(Real (*a_func)(const RealVect a_pos)){
  m_func = a_func;
}

Real WallBc::get_value(){
  Real value = 0.0;
  
  if(m_which == wallbc::dirichlet){
    value = Real(m_live);
  }
  else if(m_which == wallbc::neumann){
    value = m_value;
  }
  else if(m_which == wallbc::robin){
    value = m_value;
  }

  return value;
}

bool WallBc::is_live(){
  return m_live;
}

wallbc::which_bc WallBc::which_bc(){
  return m_which;
}

int WallBc::map_bc(const int a_dir, const Side::LoHiSide a_side) {
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}

#include "CD_NamespaceFooter.H"
