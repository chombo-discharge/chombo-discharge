/*!
  @file wall_bc.cpp
  @brief Implementation of wall_bc.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "wall_bc.H"

#include "CD_NamespaceHeader.H"

wall_bc::wall_bc(const int a_dir, const Side::LoHiSide a_side, wallbc::which_bc a_which){
  m_dir   = a_dir;
  m_side  = a_side;
  m_which = a_which;
}

wall_bc::~wall_bc(){
}

void wall_bc::set_value(Real a_value){
  m_value = a_value;
}

void wall_bc::set_live(bool a_live){
  m_live = a_live;
}

void wall_bc::set_function(Real (*a_func)(const RealVect a_pos)){
  m_func = a_func;
}

Real wall_bc::get_value(){
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

bool wall_bc::is_live(){
  return m_live;
}

wallbc::which_bc wall_bc::which_bc(){
  return m_which;
}

int wall_bc::map_bc(const int a_dir, const Side::LoHiSide a_side) {
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}

#include "CD_NamespaceFooter.H"
