/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
#include <algorithm>
#include <math.h>

// Our includes
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>

template <int N>
inline
LookupTable<N>::LookupTable(){
  m_state     = std::make_tuple(-1, false, false);
}

template <int N>
inline
LookupTable<N>::~LookupTable(){
  
}

template <int N>
inline
LookupTable<N>::LookupTable(const LookupTable& a_table){
  m_data      = a_table.m_data;
  m_state     = a_table.m_state;
  m_delta     = a_table.m_delta;  
}

template <int N>
template <typename... Ts>
inline
void LookupTable<N>::addEntry(const Ts&... x){
  std::array<Real, sizeof...(Ts)> arr = {x...};
  m_data.emplace_back(arr);

  std::get<Sorted> (m_state) = false;
  std::get<Uniform>(m_state) = false;  
}

template <int N>
inline
void LookupTable<N>::sort(const int a_whichColumn){

  auto comparator = [idx = a_whichColumn](const std::array<Real, N>& a, const std::array<Real, N>& b) -> bool {
    return a[idx] < b[idx];
  };
  
  std::sort(m_data.begin(), m_data.end(), comparator);

  std::get<Column>(m_state) = a_whichColumn;
  std::get<Sorted>(m_state) = true;
}

template <int N>
template <int K>
inline
void LookupTable<N>::scale(const Real a_scale){
  for (auto& r : m_data){
    r[K] *= a_scale;
  }

  if(K == std::get<Column>(m_state)){
    m_delta *= a_scale;
  }
}

template <int N>
inline
void LookupTable<N>::swap(const int a_column1, const int a_column2){
  for (auto& r : m_data){
    const Real tmp = r[a_column1];

    r[a_column1] = r[a_column2];
    r[a_column2] = tmp;
  }

  if(std::get<Sorted>(m_state)){
    
    if(a_column1 == std::get<Column>(m_state)){
      std::get<Column> (m_state) = a_column1;
    }
    else if(a_column2 == std::get<Column>(m_state)){
      std::get<Column> (m_state) = a_column2;
    }
  }
}

template <int N>
inline
void LookupTable<N>::clear() {
  m_data.resize(0);
}

template <int N>
inline
void LookupTable<N>::setRange(const Real a_min, const Real a_max, const int a_colRange){

  // Do a backup of the data and reset data holder. 
  const std::vector<std::array<Real, N> > oldData = m_data;
  this->clear();
  
  for (const auto& d : oldData){
    if(d[a_colRange] >= a_min && d[a_colRange] <= a_max){
      m_data.emplace_back(d);
    }
  }

  std::get<Sorted> (m_state) = false;
  std::get<Uniform>(m_state) = false;      
}

template <int N>
inline
void LookupTable<N>::dumpTable(std::ostream& a_outputStream) const{
  for (const auto& r : m_data){
    for (const auto& c : r){
      a_outputStream << c << "\t";
    }
    a_outputStream << "\n";
  }
}

template <int N>
template <int K>
inline
Real LookupTable<N>::getEntry(const Real a_x) const {
  Real ret;
  
  if(m_data.size() == 1){
    ret = m_data.front()[K];
  }
  else{
    const int whichColumn = std::get<Column>(m_state);

    const int idx = std::floor( (a_x - m_data.front()[whichColumn])/m_delta);
    if (idx <= 0){
      ret = m_data.front()[K];
    }
    else if (idx >= m_data.size()-1){
      ret = m_data.back()[K];
    }
    else{
      const Real t = (a_x - m_data[idx][whichColumn])/m_delta;
      ret = m_data[idx][K] + t*(m_data[idx+1][K] - m_data[idx][K]);
    }
  }

  return ret;
}

template <int N>
inline
std::array<Real, N> LookupTable<N>::directLookup(const Real a_x) const {

  std::array<Real, N> ret;
  if(!std::get<Sorted>(m_state)){
    std::cerr << "LookupTable<N>::directLookup -- logic bust. Data must be sorted along the input column!\n";
  }
  else{
    const int column = std::get<Column>(m_state);
      
    if(a_x <= m_data.front()[column]){
      ret = m_data.front();
    }
    else if(a_x >= m_data.back()[column]){
      ret = m_data.back();
    }
    else{ // Do linear interpolation. We look through the rows and find two entries that bracket a_x. We then linearly interpolate between those entries. 
      for (int irow = 0; irow < m_data.size()-1; irow++){
	const std::array<Real, N>& row1 = m_data[irow  ];
	const std::array<Real, N>& row2 = m_data[irow+1];

	const Real x1 = row1[column];
	const Real x2 = row2[column];
	
	if(a_x >= x1 && a_x <= x2){
	  const Real t = (a_x - x1)/(x2-x1);

	  // Now interpolate data in each column. 
	  for (int j = 0; j < N; j++){
	    ret[j] = row1[j] + t * (row2[j] - row1[j]);
	  }
	}
      }
    }
  }

  return ret;
}

template <int N>
inline
void LookupTable<N>::makeUniform(const int a_numRows){

  if(std::get<Sorted>(m_state)){
    const int whichColumn = std::get<Column>(m_state);

    if(a_numRows > 1){

      // Find maximum value and spacing. 
      Real maxVal = -std::numeric_limits<Real>::max();
      Real minVal =  std::numeric_limits<Real>::max();

      for (const auto& r : m_data){
	maxVal = std::max(maxVal, r[whichColumn]);
	minVal = std::min(minVal, r[whichColumn]);
      }

      m_delta = (maxVal - minVal)/(a_numRows - 1);

      // Do direct lookup of the data. 
      std::vector<std::array<Real, N> > newData;

      for (int i = 0; i < a_numRows; i++){
	const Real x = minVal + i*m_delta;

	newData.emplace_back(this->directLookup(x));
      }

      m_data = newData;

      std::get<Uniform>(m_state) = true;
    }
    else{
      std::cerr << "LookupTable<N>::makeUniform -- must have at least two rows!\n";
    }
  }
  else{
    std::cerr << "LookupTable<N>::makeUniform -- must sort the column along one of the coordinates first!\n";
  }

}

template <int N>
inline
int LookupTable<N>::getNumEntries() const {
  return m_data.size();
}

#include <CD_NamespaceFooter.H>

#endif
