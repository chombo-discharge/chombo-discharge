/*!
  @file wall_bc.cpp
  @brief Implementation of wall_bc.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "wall_bc.H"

wall_bc::wall_bc(const int a_dir, const Side::LoHiSide a_side, wallbc::which_bc a_which){
  m_dir   = a_dir;
  m_side  = a_side;
  m_which = a_which;
}

wall_bc::~wall_bc(){
}

void wall_bc::set_value(Real a_value){
  CH_assert(m_which == wallbc::neumann);
  m_value = a_value;
}

void wall_bc::set_live(bool a_live){
  CH_assert(m_which == wallbc::dirichlet);
  m_live = a_live;
}

wallbc::which_bc wall_bc::which_bc(){
  return m_which;
}

int wall_bc::map_bc(const int a_dir, const Side::LoHiSide a_side) {
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}
