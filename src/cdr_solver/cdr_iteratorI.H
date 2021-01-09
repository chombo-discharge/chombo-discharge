/*!
  @file   cdr_iteratorI.H
  @brief  Implementation of cdr_iterator.H
  @author Robert Marskar
  @date   June Dec. 2017
*/

#ifndef _CDR_ITERATORI_H_
#define _CDR_ITERATORI_H_

#include "cdr_iterator.H"

template <class T>
cdr_iterator<T>::cdr_iterator(){
  CH_TIME("cdr_iterator<T>::cdr_iterator");
  m_defined = false;
}

template <class T>
cdr_iterator<T>::cdr_iterator(cdr_layout<T>& a_layout, const species_iteration a_mode){
  CH_TIME("cdr_iterator<T>::cdr_iterator");
    
  m_solvers  = a_layout.get_solvers();
  m_species  = a_layout.get_species();
  m_num      = m_solvers.size();
  m_mode     = a_mode;
  m_defined  = true;

  reset();
}

template <class T>
cdr_iterator<T>::cdr_iterator(cdr_layout<T>& a_layout, Vector<EBAMRCellData*> a_data){
  CH_TIME("cdr_iterator<T>::cdr_iterator");
    
  m_solvers  = a_layout.get_solvers();
  m_species  = a_layout.get_species();
  m_num      = m_solvers.size();
  m_defined  = true;
  m_mode     = species_iteration::all;

  m_celldata = a_data;

  reset();
}

template <class T>
cdr_iterator<T>::~cdr_iterator(){

}

template <class T>
int cdr_iterator<T>::num_solvers(){
  CH_assert(m_defined);
  return m_num;
}

template <class T>
int cdr_iterator<T>::index() const {
  CH_assert(m_defined);
  return m_isolver;
}

template <class T>
void cdr_iterator<T>::reset(){
  CH_assert(m_defined);
  m_isolver = 0;
}

template <class T>
bool cdr_iterator<T>::ok(){
  return (m_isolver < m_num);
}

template <class T>
void cdr_iterator<T>::operator++(){
  CH_assert(m_defined);

  if(m_mode == species_iteration::all){
    m_isolver++;
  }
  else{
    m_isolver++;

    for (int isolver = m_isolver; this->ok(); ++isolver){
      const RefCountedPtr<T>& solver            = m_solvers[isolver];
      const RefCountedPtr<cdr_species>& species = m_species[isolver];
      
      if(m_mode == species_iteration::charged){
	if(species->get_charge() != 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == species_iteration::negative){
	if(species->get_charge() < 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == species_iteration::positive){
	if(species->get_charge() > 0){
	  m_isolver = isolver;
	  break;
	}
      }
    }
  }
}

template <class T>
RefCountedPtr<T>& cdr_iterator<T>::operator() () {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_solvers[m_isolver];
}

template <class T>
RefCountedPtr<cdr_species>& cdr_iterator<T>::get_species() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_species[m_isolver];
}

template <class T>
EBAMRCellData& cdr_iterator<T>::get_data() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return *m_celldata[m_isolver];
}

#endif
