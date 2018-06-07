/*!
  @file   geo_coarsener.cpp
  @brief  Implementation of geo_coarsener.H
  @author Robert Marskar
  @date   June 2018
*/

#include "geo_coarsener.H"


geo_coarsener::geo_coarsener(){
  m_coarsen_boxes.resize(0);
  m_coarsen_levels.resize(0);
}

geo_coarsener::~geo_coarsener(){
}

Vector<real_box> geo_coarsener::get_coarsen_boxes(){
  return m_coarsen_boxes;
}


Vector<int> geo_coarsener::get_coarsen_levels(){
  return m_coarsen_levels;
}
