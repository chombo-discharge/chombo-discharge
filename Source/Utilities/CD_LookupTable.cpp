/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LookupTable.cpp
  @brief  Implementation of CD_LookupTable.H
  @author Robert Marskar
*/

// Our includes
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>
  
LookupTable::LookupTable(){

  m_numEntries = 0;
  m_x.resize(m_numEntries);
  m_y.resize(m_numEntries);
}

LookupTable::~LookupTable(){}

LookupTable::LookupTable(const LookupTable& a_table){
  m_dx = a_table.m_dx;
  m_numEntries = a_table.m_numEntries;
  m_x = a_table.m_x;
  m_y = a_table.m_y;
}

#include <CD_NamespaceFooter.H>
