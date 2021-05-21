/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_MFInterfaceFABImplem.H
  @brief  Implementation of CD_MFInterfaceFAB.H
  @author Robert Marskar
*/

#ifndef CD_MFInterfaceFABImplem_H
#define CD_MFInterfaceFABImplem_H

// Our includes
#include <CD_MFInterfaceFAB.H>
#include <CD_NamespaceHeader.H>

template <class T>
MFInterfaceFAB<T>::MFInterfaceFAB(){
}

template <class T>
MFInterfaceFAB<T>::MFInterfaceFAB(const MFLevelGrid& a_mflg, const DataIndex& a_dit){
  this->define(a_mflg, a_dit);
}

template <class T>
MFInterfaceFAB<T>::~MFInterfaceFAB(){
  for (int i = 0; i < m_phase.size(); i++){
    delete m_phase[i];
  }
}

template <class T>
void MFInterfaceFAB<T>::define(const MFLevelGrid& a_mflg, const DataIndex& a_dit){

  const DisjointBoxLayout& dbl    = a_mflg.getGrids();
  const ProblemDomain& domain     = a_mflg.get_domain();
  const RefCountedPtr<MultiFluidIndexSpace>& MultiFluidIndexSpace = a_mflg.getMfIndexSpace();
  const int numPhases            = MultiFluidIndexSpace->numPhases();

  m_ivs = a_mflg.interfaceRegion(dbl.get(a_dit), a_dit);
  m_phase.resize(numPhases, NULL);
    
  for (int i = 0; i < numPhases; i++){
    const EBLevelGrid& eblg = a_mflg.getEBLevelGrid(i);
    const EBISLayout& ebisl = eblg.getEBISL();
    m_phase[i] = new BaseIVFAB<T>(m_ivs, ebisl[a_dit].getEBGraph(), 1);
  }
}

template <class T>
const IntVectSet& MFInterfaceFAB<T>::getIVS() const {
  return m_ivs;
}

template <class T>
BaseIVFAB<T>& MFInterfaceFAB<T>::getIVFAB(const int a_phase){
  return *m_phase[a_phase];
}

template <class T>
const BaseIVFAB<T>& MFInterfaceFAB<T>::getIVFAB(const int a_phase) const {
  return *m_phase[a_phase];
}

template <class T>
BaseIVFAB<T>* MFInterfaceFAB<T>::getPhasePtr(int a_phase){
  return m_phase[a_phase];
}

template <class T>
int MFInterfaceFAB<T>::numPhases() const{
  return m_phase.size();
}

template <class T>
void MFInterfaceFAB<T>::setVal(const T& a_value){
  for (int i = 0; i < m_phase.size(); i++){
    m_phase[i]->setVal(a_value);
  }
}

#include <CD_NamespaceFooter.H>

#endif
