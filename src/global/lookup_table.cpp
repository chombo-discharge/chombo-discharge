/*!
  @file   lookup_table.cpp
  @brief  Implementation of lookup_table.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "lookup_table.H"

lookup_table::lookup_table(){

  m_num_entries = 0;
  m_x.resize(m_num_entries);
  m_y.resize(m_num_entries);
}

lookup_table::~lookup_table(){
}


lookup_table::lookup_table(const lookup_table& a_table){
  m_dx = a_table.m_dx;
  m_num_entries = a_table.m_num_entries;
  m_x = a_table.m_x;
  m_y = a_table.m_y;
}

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

void lookup_table::scale_x(const Real& a_scale){
  for (int i = 0; i < m_x.size(); i++){
    m_x[i] *= a_scale;
  }
}

void lookup_table::scale_y(const Real& a_scale){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] *= a_scale;
  }
}

void lookup_table::dump_table(){
  for (int i = 0; i < m_y.size(); i++){
    pout() << m_x[i] << "\t" << m_y[i] << endl;
  }
}

void lookup_table::operator+=(const lookup_table& a_table){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] += a_table.m_y[i];
  }
}

void lookup_table::operator-=(const lookup_table& a_table){
  for (int i = 0; i < m_y.size(); i++){
    m_y[i] -= a_table.m_y[i];
  }
}

Real lookup_table::get_entry(const Real a_x){

  Real value;
  
  if(m_num_entries == 1){
    value = m_y[0];
  }
  else {

    // Find entry
    int i = floor((a_x - m_x[0])/m_dx);

    if(i < 0){
      i = 0;
      pout() << "lookup_table::get_entry - entry excceds range (low end)" << endl;
    }
    else if(i >= m_num_entries - 1){
      i = m_num_entries - 2;
      pout() << "lookup_table::get_entry - entry excceds range (high end)" << endl;
    }

    const int i_hi = i + 1;
    const Real x0  = m_x[i];
    const Real x1  = m_x[i_hi];
    const Real y0  = m_y[i];
    const Real y1  = m_y[i_hi];

    value = y0 + ((y1-y0)/(x1-x0))*(a_x-x0);

  }

  return value;
}

const Real lookup_table::get_entry(const Real a_x) const {

  Real value;
  
  if(m_num_entries == 1){
    value = m_y[0];
  }
  else {

    // Find entry
    int i = floor((a_x - m_x[0])/m_dx);

    if(i < 0){
      i = 0;
      pout() << "lookup_table::get_entry - entry excceds range (low end)" << endl;
    }
    else if(i >= m_num_entries - 1){
      i = m_num_entries - 2;
      pout() << "lookup_table::get_entry - entry excceds range (high end)" << endl;
    }

    const int i_hi = i + 1;
    const Real x0  = m_x[i];
    const Real x1  = m_x[i_hi];
    const Real y0  = m_y[i];
    const Real y1  = m_y[i_hi];

    value = y0 + ((y1-y0)/(x1-x0))*(a_x-x0);

  }

  return value;
}
