/*!
  @file   MultiIndex.cpp
  @brief  Implementation of MultiIndex.H
  @author Robert Marskar
  @date   Aug. 2017
*/

#include "MultiIndex.H"

MultiIndex::MultiIndex(const int a_Q){
  this->define(IntVect::Zero, a_Q);
}

MultiIndex::MultiIndex(const IntVect a_init, const int a_Q){
  this->define(a_init, a_Q);
}


void MultiIndex::define(const IntVect a_index, const int a_Q){
  m_curIndex = a_index;
  m_Q        = a_Q;

  //  Now define the map. 
}
  

IntVect MultiIndex::getIntVect() const {
  return m_curIndex;
}

int MultiIndex::getQ() const {
  return m_Q;
}

int MultiIndex::getLinearIndex(const IntVect a_multiIndex) const {
  return m_mapToLinearIndex.at(a_multiIndex);
}

IntVect MultiIndex::getMultiIndex(const int a_linearIndex) const {
  return m_mapToMultiIndex.at(a_linearIndex);
}

int MultiIndex::operator[](const int a_dir) const {
  return m_curIndex[a_dir];
}

bool MultiIndex::operator==(const MultiIndex& a_index) const {
  return D_TERM(m_curIndex[0] == a_index[0], && m_curIndex[1] == a_index[1], && m_curIndex[2] == a_index[2]);
}

bool MultiIndex::operator==(const int a_Q) const {
  return this->norm() == a_Q;
}

bool MultiIndex::operator>(const MultiIndex& a_index) const {
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if(m_curIndex[dir] > a_index[dir]){
      return true;
    }
    else if(m_curIndex[dir] < a_index[dir]){
      return false;
    }
  }

  return false;
}

bool MultiIndex::operator>(const int a_Q) const {
  return this->norm() > a_Q;
}

bool MultiIndex::operator<(const MultiIndex& a_index) const {
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if(m_curIndex[dir] < a_index[dir]){
      return true;
    }
    else if(m_curIndex[dir] > a_index[dir]){
      return false;
    }
  }

  return false;
}

bool MultiIndex::operator<(const int a_Q) const {
  return this->norm() < a_Q;
}

bool MultiIndex::operator<=(const int a_Q) const {
  if(this->norm() < a_Q || this->norm() == a_Q) {
    return true;
  }
  else{
    return false;
  }
}

bool MultiIndex::ok() const {
  return (this->norm() <= m_Q);
}

int MultiIndex::factorial(const int a_n) const {
  if(a_n > 1){
    return a_n*factorial(a_n - 1);
  }
  else{
    return 1;
  };
}

int MultiIndex::factorial() const {
  int fact = 1;
  for (int dir = 0; dir < SpaceDim; dir++){
    fact *= factorial(m_curIndex[dir]);
  }

  return fact;
}

int MultiIndex::norm() const {
  int retval = 0;
  for (int dir = 0; dir < SpaceDim; dir++){
    retval += m_curIndex[dir];
  }

  return retval;
}

Real MultiIndex::pow(const RealVect& a_vec){

  Real retval = 1.;

  for (int dir = 0; dir < SpaceDim; dir++){
    double base = a_vec[dir];
    double exp  = 1.0*m_curIndex[dir];
    retval *= std::pow(base, exp);
  }

  return retval;
}

void MultiIndex::operator++(){

  if(this->norm() < m_Q) { // Can raise first index.
    m_curIndex[0]++;
  }
  else if(this->norm() == m_Q){ // Can't raise first order, check next index.

    MultiIndex newIndex(IntVect(D_DECL(0, m_curIndex[1]+1, m_curIndex[2])), m_Q);
    if(newIndex.norm() <= m_Q){ // Ok, this is a valid index.
      *this = newIndex;
    }
#if CH_SPACEDIM==3
    else if(newIndex.norm() > m_Q){
      MultiIndex newIndex = MultiIndex(IntVect(D_DECL(0, 0, m_curIndex[2]+1)), m_q);
      if(newIndex.norm() <= m_Q){
	*this = newIndex;
      }
      else{
	m_curIndex[0]++; // Designed to break out of ::ok() loop.
      }
    }
#endif
    else{
      m_curIndex[0]++; // Designed to break out of ::ok() loop.
    }
  }
}

