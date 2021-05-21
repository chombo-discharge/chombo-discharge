/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ItoIteratorImplem.H
  @brief  Implementation of CD_ItoIterator.H
  @author Robert Marskar
*/

#ifndef CD_ItoIteratorImplem_H
#define CD_ItoIteratorImplem_H

// Our includes
#include <CD_ItoIterator.H>
#include <CD_NamespaceHeader.H>

template <class T>
ItoIterator<T>::ItoIterator(){
  CH_TIME("ItoIterator<T>::ItoIterator<T>");
  m_defined = false;
}

template <class T>
ItoIterator<T>::ItoIterator(ItoLayout<T>& a_layout, const SpeciesIteration a_mode){
  CH_TIME("ItoIterator<T>::ItoIterator<T>");
    
  m_solvers  = a_layout.getSolvers();
  m_species  = a_layout.getSpecies();
  m_num      = m_solvers.size();
  m_mode     = a_mode;
  m_defined  = true;

  reset();
}

template <class T>
ItoIterator<T>::~ItoIterator<T>(){

}

template <class T>
int ItoIterator<T>::getNumberOfSolvers(){
  CH_assert(m_defined);
  return m_num;
}

template <class T>
int ItoIterator<T>::index() const {
  CH_assert(m_defined);
  return m_isolver;
}

template <class T>
void ItoIterator<T>::reset(){
  CH_assert(m_defined);
  m_isolver = 0;
}

template <class T>
bool ItoIterator<T>::ok(){
  return (m_isolver < m_num);
}

template <class T>
void ItoIterator<T>::operator++(){
  CH_assert(m_defined);

  if(m_mode == SpeciesIteration::All){
    m_isolver++;
  }
  else{
    m_isolver++;

    for (int isolver = m_isolver; this->ok(); ++isolver){
      const RefCountedPtr<ItoSolver>& solver = m_solvers[isolver];
      const RefCountedPtr<ItoSpecies>& species   = m_species[isolver];
      
      if(m_mode == SpeciesIteration::Charged){
	if(species->getChargeNumber() != 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == SpeciesIteration::Negative){
	if(species->getChargeNumber() < 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == SpeciesIteration::Positive){
	if(species->getChargeNumber() > 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == SpeciesIteration::Neutral){
	if(species->getChargeNumber() == 0){
	  m_isolver = isolver;
	  break;
	}
      }
    }
  }
}

template <class T>
RefCountedPtr<T>& ItoIterator<T>::operator() () {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_solvers[m_isolver];
}

template <class T>
RefCountedPtr<ItoSpecies>& ItoIterator<T>::getSpecies() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_species[m_isolver];
}

#include <CD_NamespaceFooter.H>

#endif
