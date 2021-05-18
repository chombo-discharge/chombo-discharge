/*!
  @file   CdrIteratorI.H
  @brief  Implementation of CdrIterator.H
  @author Robert Marskar
  @date   June Dec. 2017
*/

#ifndef _CDR_ITERATORI_H_
#define _CDR_ITERATORI_H_

#include <CD_CdrIterator.H>

#include "CD_NamespaceHeader.H"

template <class T>
CdrIterator<T>::CdrIterator(){
  CH_TIME("CdrIterator<T>::CdrIterator");
  m_defined = false;
}

template <class T>
CdrIterator<T>::CdrIterator(cdr_layout<T>& a_layout, const SpeciesIteration a_mode){
  CH_TIME("CdrIterator<T>::CdrIterator");
    
  m_solvers  = a_layout.get_solvers();
  m_species  = a_layout.get_species();
  m_num      = m_solvers.size();
  m_mode     = a_mode;
  m_defined  = true;

  reset();
}

template <class T>
CdrIterator<T>::CdrIterator(cdr_layout<T>& a_layout, Vector<EBAMRCellData*> a_data){
  CH_TIME("CdrIterator<T>::CdrIterator");
    
  m_solvers  = a_layout.get_solvers();
  m_species  = a_layout.get_species();
  m_num      = m_solvers.size();
  m_defined  = true;
  m_mode     = SpeciesIteration::all;

  m_celldata = a_data;

  reset();
}

template <class T>
CdrIterator<T>::~CdrIterator(){

}

template <class T>
int CdrIterator<T>::num_solvers(){
  CH_assert(m_defined);
  return m_num;
}

template <class T>
int CdrIterator<T>::index() const {
  CH_assert(m_defined);
  return m_isolver;
}

template <class T>
void CdrIterator<T>::reset(){
  CH_assert(m_defined);
  m_isolver = 0;
}

template <class T>
bool CdrIterator<T>::ok(){
  return (m_isolver < m_num);
}

template <class T>
void CdrIterator<T>::operator++(){
  CH_assert(m_defined);

  if(m_mode == SpeciesIteration::all){
    m_isolver++;
  }
  else{
    m_isolver++;

    for (int isolver = m_isolver; this->ok(); ++isolver){
      const RefCountedPtr<T>& solver            = m_solvers[isolver];
      const RefCountedPtr<CdrSpecies>& species = m_species[isolver];
      
      if(m_mode == SpeciesIteration::charged){
	if(species->getChargeNumber() != 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == SpeciesIteration::negative){
	if(species->getChargeNumber() < 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == SpeciesIteration::positive){
	if(species->getChargeNumber() > 0){
	  m_isolver = isolver;
	  break;
	}
      }
    }
  }
}

template <class T>
RefCountedPtr<T>& CdrIterator<T>::operator() () {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_solvers[m_isolver];
}

template <class T>
RefCountedPtr<CdrSpecies>& CdrIterator<T>::get_species() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_species[m_isolver];
}

template <class T>
EBAMRCellData& CdrIterator<T>::get_data() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return *m_celldata[m_isolver];
}
#include "CD_NamespaceFooter.H"

#endif
