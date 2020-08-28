/*!
  @file   lookup_tableI.H
  @brief  Implementation of lookup_table.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _LOOKUP_TABLEI_H_
#define _LOOKUP_TABLEI_H_

#include "lookup_table.H"

#include <iostream>

inline
void lookup_table::add_entry(const Real& a_x, const Real& a_y){
  m_x.push_back(a_x);
  m_y.push_back(a_y);
  m_num_entries++;

  if(m_num_entries == 1){
    m_dx    = 0.;
  }
  else if(m_num_entries >= 2){
    m_dx    = m_x[1] - m_x[0];
  }
}

inline
void lookup_table::add_table(const lookup_table& a_otherTable, const Real a_scale){
  for (int i = 0; i < m_x.size(); i++){
    m_y[i] += a_scale*a_otherTable.get_entry(m_x[i]);
  }
}

inline
void lookup_table::scale_x(const Real& a_scale){
  for (int i = 0; i < m_x.size(); i++){
    m_x[i] *= a_scale;
  }
}

inline
void lookup_table::scale_y(const Real& a_scale){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] *= a_scale;
  }
}

inline
void lookup_table::dump_table() const{
  for (int i = 0; i < m_y.size(); i++){
    std::cout << m_x[i] << "\t" << m_y[i] << std::endl;
  }
}

inline
void lookup_table::swap_xy(){
  std::vector<Real> tmp = m_x;
  m_x = m_y;
  m_y = tmp;

  if(m_num_entries == 1){
    m_dx    = 0.;
  }
  else if(m_num_entries >= 2){
    m_dx    = m_x[1] - m_x[0];
  }
}

inline
void lookup_table::operator+=(const lookup_table& a_table){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] += a_table.m_y[i];
  }
}

inline
void lookup_table::operator-=(const lookup_table& a_table){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] -= a_table.m_y[i];
  }
}

inline
void lookup_table::operator=(const lookup_table& a_table){
  m_num_entries = a_table.m_num_entries;
  m_dx = a_table.m_dx;
  m_x = a_table.m_x;
  m_y = a_table.m_y;
}

inline
Real lookup_table::get_entry(const Real a_x) const {

  Real value;
  if(m_num_entries == 1){
    value = m_y[0];
  }
  else {

    // Find entry. We will linear interpolate between i and i+1
    // where i: 
    const int i = floor((a_x - m_x[0])/m_dx);

    if(i < 0){
      value = m_y[0];
    }
    else if(i >= m_num_entries - 1){
      value = m_y[m_num_entries - 1];
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
Real lookup_table::direct_lookup(const Real a_x) const {
  Real value;
  if(a_x <= m_x[0]){
    value = m_y[0];
  }
  else if(a_x >= m_x[m_num_entries - 1]){
    value = m_y[m_num_entries - 1];
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
void lookup_table::make_uniform(const int a_num_entries){

  std::vector<Real> new_x(a_num_entries, 0.0);
  std::vector<Real> new_y(a_num_entries, 0.0);

  const Real xmin = m_x[0];
  const Real xmax = m_x[m_num_entries-1];
  const Real dx   = (xmax - xmin)/(a_num_entries - 1);

  // Build x
  for (int i = 0; i < new_x.size(); i++){
    new_x[i] = xmin + i*dx;
  }

  // Build y
  for (int i = 0; i < new_y.size(); i++){
    new_y[i] = direct_lookup(new_x[i]);
  }


  // Copy to internal data
  m_dx = dx;
  m_x  = new_x;
  m_y  = new_y;
  m_num_entries = a_num_entries;
}

#endif
