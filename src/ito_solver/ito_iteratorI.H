/*!
  @file   ito_iterator.H
  @brief  Implementation of ito_iterator.H
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _ITO_ITERATORI_H_
#define _ITO_ITERATORI_H_

#include "ito_iterator.H"

#include "CD_NamespaceHeader.H"

template <class T>
ito_iterator<T>::ito_iterator(){
  CH_TIME("ito_iterator<T>::ito_iterator<T>");
  m_defined = false;
}

template <class T>
ito_iterator<T>::ito_iterator(ito_layout<T>& a_layout, const species_iteration a_mode){
  CH_TIME("ito_iterator<T>::ito_iterator<T>");
    
  m_solvers  = a_layout.get_solvers();
  m_species  = a_layout.get_species();
  m_num      = m_solvers.size();
  m_mode     = a_mode;
  m_defined  = true;

  reset();
}

template <class T>
ito_iterator<T>::~ito_iterator<T>(){

}

template <class T>
int ito_iterator<T>::num_solvers(){
  CH_assert(m_defined);
  return m_num;
}

template <class T>
int ito_iterator<T>::index() const {
  CH_assert(m_defined);
  return m_isolver;
}

template <class T>
void ito_iterator<T>::reset(){
  CH_assert(m_defined);
  m_isolver = 0;
}

template <class T>
bool ito_iterator<T>::ok(){
  return (m_isolver < m_num);
}

template <class T>
void ito_iterator<T>::operator++(){
  CH_assert(m_defined);

  if(m_mode == species_iteration::all){
    m_isolver++;
  }
  else{
    m_isolver++;

    for (int isolver = m_isolver; this->ok(); ++isolver){
      const RefCountedPtr<ito_solver>& solver = m_solvers[isolver];
      const RefCountedPtr<ito_species>& species   = m_species[isolver];
      
      if(m_mode == species_iteration::charged){
	if(species->getChargeNumber() != 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == species_iteration::negative){
	if(species->getChargeNumber() < 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == species_iteration::positive){
	if(species->getChargeNumber() > 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == species_iteration::neutral){
	if(species->getChargeNumber() == 0){
	  m_isolver = isolver;
	  break;
	}
      }
    }
  }
}

template <class T>
RefCountedPtr<T>& ito_iterator<T>::operator() () {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_solvers[m_isolver];
}

template <class T>
RefCountedPtr<ito_species>& ito_iterator<T>::get_species() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_species[m_isolver];
}
#include "CD_NamespaceFooter.H"

#endif
