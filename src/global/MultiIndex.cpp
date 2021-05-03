/*!
  @file   MultiIndex.cpp
  @brief  Implementation of MultiIndex.H
  @author Robert Marskar
  @date   Aug. 2017
*/

#include "MultiIndex.H"

MultiIndex::MultiIndex(){
  m_isDefined = false;
}

MultiIndex::MultiIndex(const IntVect& a_index){
  define(a_index);
}

MultiIndex::MultiIndex(const MultiIndex& a_index){
  define(a_index.getIntVect());
}

MultiIndex::~MultiIndex(){
}

void MultiIndex::define(const IntVect& a_index){
  m_index     = a_index;
  m_isDefined = true;
}

IntVect MultiIndex::getIntVect() const {
  return m_index;
}

int MultiIndex::operator[](const int a_dir) const {
  return m_index[a_dir];
}

MultiIndex& MultiIndex::operator=(const MultiIndex& a_index){
  define(a_index.getIntVect());
}

bool MultiIndex::operator==(const MultiIndex& a_index) const {
  return D_TERM(m_index[0] == a_index[0], && m_index[1] == a_index[1], && m_index[2] == a_index[2]);
}

bool MultiIndex::operator==(const int a_Q) const {
  return abs(*this) == a_Q;
}

bool MultiIndex::operator>(const MultiIndex& a_index) const {
  
  //
  for (int dir = 0; dir < SpaceDim; dir++){
    if(m_index[dir] > a_index[dir]){
      return true;
    }
    else if(m_index[dir] < a_index[dir]){
      return false;
    }
  }

  return false;
}

bool MultiIndex::operator>(const int a_Q) const {
  return abs(*this) > a_Q;
}

bool MultiIndex::operator<(const MultiIndex& a_index) const {
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if(m_index[dir] < a_index[dir]){
      return true;
    }
    else if(m_index[dir] > a_index[dir]){
      return false;
    }
  }

  return false;
}

bool MultiIndex::operator<(const int a_Q) const {
  return abs(*this) < a_Q;
}

bool MultiIndex::operator<=(const int a_Q) const {
  if(abs(*this) < a_Q || abs(*this) == a_Q) {
    return true;
  }
  else{
    return false;
  }
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
  CH_assert(m_isDefined);

  int fact = 1;
  for (int dir = 0; dir < SpaceDim; dir++){
    fact *= factorial(m_index[dir]);
  }

  return fact;
}

int MultiIndex::abs(const MultiIndex& a_index){
  int retval = 0;
  for (int dir = 0; dir < SpaceDim; dir++){
    retval += a_index[dir];
  }

  return retval;
}

Real MultiIndex::pow(const RealVect& a_vec, const MultiIndex& a_index){

  Real retval = 1.;

  for (int dir = 0; dir < SpaceDim; dir++){
    double base = a_vec[dir];
    double exp  = 1.0*a_index[dir];
    retval *= std::pow(base, exp);
  }

  return retval;
}

void MultiIndex::next(const int a_Q){
  if(abs(*this) < a_Q){
    m_index[0] += 1;
  }
  else if(abs(*this) == a_Q){ // Cannot increase first index.

    // Try second index. 
    MultiIndex newIndex(IntVect(D_DECL(0, m_index[1] + 1, m_index[2])));
    if(abs(newIndex) <= a_Q){
      *this = newIndex; // OK, valid index
    }
    
#if CH_SPACEDIM == 3
    // Could not increase second index, try third index. 
    else if(abs(newIndex) > a_Q){
      newIndex = MultiIndex(IntVect(D_DECL(0, 0, m_index[2] + 1)));
      if(abs(newIndex) <= a_Q){
	*this = newIndex;
      }
      else{
	m_index[0] += 1;
      }
    }
#endif
    else{
      m_index[0] += 1;
    }
  }
}

//
bool MultiIndex::hasNext(const int a_Q){
  if(abs(*this) < a_Q){
    return true;
  }
  else if(abs(*this) == a_Q){ // Cannot increase first index. Try second index. 
    MultiIndex newIndex(IntVect(D_DECL(0, m_index[1] + 1, m_index[2])));
    if(abs(newIndex) <= a_Q){
      return true;
    }

    
#if CH_SPACEDIM == 3
    else if(abs(newIndex) > a_Q){
      newIndex = MultiIndex(IntVect(D_DECL(0, 0, m_index[2] + 1)));
      if(abs(newIndex) <= a_Q){
	return true;
      }
    }
#endif
  }

  return false;
}
