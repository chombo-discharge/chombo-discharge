/*!
  @file   lookup_table.cpp
  @brief  Implementation of lookup_table.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "lookup_table.H"

#include "CD_NamespaceHeader.H"
  
lookup_table::lookup_table(){

  m_num_entries = 0;
  m_x.resize(m_num_entries);
  m_y.resize(m_num_entries);
}

lookup_table::~lookup_table(){}

lookup_table::lookup_table(const lookup_table& a_table){
  m_dx = a_table.m_dx;
  m_num_entries = a_table.m_num_entries;
  m_x = a_table.m_x;
  m_y = a_table.m_y;
}
#include "CD_NamespaceFooter.H"
