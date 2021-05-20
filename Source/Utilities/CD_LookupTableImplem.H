/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_LookupTableImplem.H
  @brief  Implementation of CD_LookupTable.H
  @author Robert Marskar
*/

#ifndef CD_LookupTableImplem_H
#define CD_LookupTableImplem_H

// Std includes
#include <iostream>

// Our includes
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>

inline
void LookupTable::addEntry(const Real& a_x, const Real& a_y){
  m_x.push_back(a_x);
  m_y.push_back(a_y);
  m_numEntries++;

  if(m_numEntries == 1){
    m_dx    = 0.;
  }
  else if(m_numEntries >= 2){
    m_dx    = m_x[1] - m_x[0];
  }
}

inline
void LookupTable::addTable(const LookupTable& a_otherTable, const Real a_scale){
  for (int i = 0; i < m_x.size(); i++){
    m_y[i] += a_scale*a_otherTable.getEntry(m_x[i]);
  }
}

inline
void LookupTable::scaleX(const Real& a_scale){
  for (int i = 0; i < m_x.size(); i++){
    m_x[i] *= a_scale;
  }

  m_dx *= a_scale;
}

inline
void LookupTable::scaleY(const Real& a_scale){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] *= a_scale;
  }
}

inline
void LookupTable::dumpTable() const{
  for (int i = 0; i < m_y.size(); i++){
    std::cout << m_x[i] << "\t" << m_y[i] << std::endl;
  }
}

inline
void LookupTable::swapXY(){
  std::vector<Real> tmp = m_x;
  m_x = m_y;
  m_y = tmp;

  if(m_numEntries == 1){
    m_dx    = 0.;
  }
  else if(m_numEntries >= 2){
    m_dx    = m_x[1] - m_x[0];
  }
}

inline
void LookupTable::operator+=(const LookupTable& a_table){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] += a_table.m_y[i];
  }
}

inline
void LookupTable::operator-=(const LookupTable& a_table){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] -= a_table.m_y[i];
  }
}

inline
void LookupTable::operator=(const LookupTable& a_table){
  m_numEntries = a_table.m_numEntries;
  m_dx = a_table.m_dx;
  m_x = a_table.m_x;
  m_y = a_table.m_y;
}

inline
Real LookupTable::getEntry(const Real a_x) const {
  Real value;
  if(m_numEntries == 1){
    value = m_y[0];
  }
  else {

    // Find entry. We will linear interpolate between i and i+1
    // where i: 
    const int i = floor((a_x - m_x[0])/m_dx);

    if(i < 0){
      value = m_y[0];
    }
    else if(i >= m_numEntries - 1){
      value = m_y[m_numEntries - 1];
    }
    else {

      const Real x0  = m_x[i];
      const Real x1  = m_x[i+1];
      const Real y0  = m_y[i];
      const Real y1  = m_y[i+1];

      const Real dydx = (y1-y0)/(x1-x0);
      value = y0 + (dydx)*(a_x - x0);
    }
  }

  return value;
}

inline
Real LookupTable::directLookup(const Real a_x) const {
  Real value;
  if(a_x <= m_x[0]){
    value = m_y[0];
  }
  else if(a_x >= m_x[m_numEntries - 1]){
    value = m_y[m_numEntries - 1];
  }
  else{
    for (int i = 0; i <= m_x.size()-2; i++){
      if(a_x >= m_x[i] && a_x <= m_x[i+1]){
	const Real x0  = m_x[i];
	const Real x1  = m_x[i+1];
	const Real y0  = m_y[i];
	const Real y1  = m_y[i+1];

	const Real dydx  = (y1-y0)/(x1-x0);
	value = y0 + dydx*(a_x-x0);

	break;
      }
    }
  }

  return value;
}

inline
void LookupTable::makeUniform(const int a_num_entries){

  std::vector<Real> new_x(a_num_entries, 0.0);
  std::vector<Real> new_y(a_num_entries, 0.0);

  const Real xmin = m_x[0];
  const Real xmax = m_x[m_numEntries-1];
  const Real dx   = (xmax - xmin)/(a_num_entries - 1);

  // Build x
  for (int i = 0; i < new_x.size(); i++){
    new_x[i] = xmin + i*dx;
  }

  // Build y
  for (int i = 0; i < new_y.size(); i++){
    new_y[i] = directLookup(new_x[i]);
  }


  // Copy to internal data
  m_dx = dx;
  m_x  = new_x;
  m_y  = new_y;
  m_numEntries = a_num_entries;
}

#include <CD_NamespaceFooter.H>

#endif
