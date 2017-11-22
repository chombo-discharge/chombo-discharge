/*!
  @file poisson_staircase_gmgs.cpp
  @brief Implementation of poisson_staircase_gmg.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_staircase_gmg.H"

poisson_staircase_gmg::poisson_staircase_gmg(){

}

poisson_staircase_gmg::~poisson_staircase_gmg(){

}

int poisson_staircase_gmg::query_ghost() const {
  return 2; // Need this many cells
}

void poisson_staircase_gmg::solve(){
  if(m_needs_setup){
    this->setup_gmg();
  }
}

void poisson_staircase_gmg::solve(MFAMRCellData& a_state, const MFAMRCellData& a_source){

}

void poisson_staircase_gmg::setup_gmg(){

}
