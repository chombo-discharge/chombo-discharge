/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtIteratorImplem.H
  @brief  Implementation of RtIterator.H
  @author Robert Marskar
*/

#ifndef CD_RtIteratorImplem_H
#define CD_RtIteratorImplem_H

// Our includes
#include <CD_RtIterator.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Iterator class for RtLayout
*/
template <class T>
RtIterator<T>::RtIterator(RtLayout<T>& a_layout){
  CH_TIME("RtIterator::RtIterator");
    
  m_solvers  = a_layout.getSolvers();
  m_species  = a_layout.getSpecies();
  m_num      = m_solvers.size();

  this->reset();
}

template <class T>
RtIterator<T>::~RtIterator(){
  CH_TIME("RtIterator::~RtIterator");
}

template <class T>
int RtIterator<T>::getNumberOfSolvers(){
  return m_num;
}

template <class T>
int RtIterator<T>::index() const {
  return m_isolver;
}

template <class T>
void RtIterator<T>::reset(){
  m_isolver = 0;
}

template <class T>
bool RtIterator<T>::ok(){
  return (m_isolver < m_num);
}


template <class T>
void RtIterator<T>::operator++(){
  m_isolver++;
}


template <class T>
RefCountedPtr<T>& RtIterator<T>::operator() () {
  CH_assert(m_isolver < m_num);
  return m_solvers[m_isolver];
}
template <class T>
RefCountedPtr<RtSpecies>& RtIterator<T>::getSpecies() {
  CH_assert(m_isolver < m_num);
  return m_species[m_isolver];
}

#include <CD_NamespaceFooter.H>

#endif
