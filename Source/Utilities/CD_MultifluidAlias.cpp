/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_MultifluidAlias.cpp
  @brief  Implementation of CD_MultifluidAlias.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MultifluidAlias.H>
#include <CD_NamespaceHeader.H>

// Public, callable function.
void MultifluidAlias::aliasMF(LevelData<EBCellFAB>& a_alias, const int a_phase, const LevelData<MFCellFAB>& a_input){
  MultifluidAlias::MfCellAliasFactory factory((LevelData<MFCellFAB>*)& a_input, a_phase);
  a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
}

void MultifluidAlias::aliasMF(LevelData<EBFluxFAB>& a_alias, const int a_phase, const LevelData<MFFluxFAB>& a_input){
  MultifluidAlias::MfFluxAliasFactory factory((LevelData<MFFluxFAB>*)& a_input, a_phase);
  a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
}

void MultifluidAlias::aliasMF(LevelData<BaseIVFAB<Real> >& a_alias, const int a_phase, const LevelData<MFBaseIVFAB>& a_input){
  MultifluidAlias::MfIVAliasFactory factory((LevelData<MFBaseIVFAB>*)& a_input, a_phase);
  a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
}

// Private cell stuff below here
MultifluidAlias::MfCellAliasFactory::MfCellAliasFactory(LevelData<MFCellFAB>* a_mf, const int a_phase){
  m_mf    = a_mf;
  m_phase = a_phase;
}

EBCellFAB* MultifluidAlias::MfCellAliasFactory::create(const Box& box, int ncomps, const DataIndex& a_datInd) const{
  return (*m_mf)[a_datInd].getPhasePtr(m_phase);
}

bool MultifluidAlias::MfCellAliasFactory::callDelete() const {
  return false;
}

// Private flux stuff below here
MultifluidAlias::MfFluxAliasFactory::MfFluxAliasFactory(LevelData<MFFluxFAB>* a_mf, const int a_phase){
  m_mf    = a_mf;
  m_phase = a_phase;
}

EBFluxFAB* MultifluidAlias::MfFluxAliasFactory::create(const Box& box, int ncomps, const DataIndex& a_datInd) const{
  return (*m_mf)[a_datInd].getPhasePtr(m_phase);
}

bool MultifluidAlias::MfFluxAliasFactory::callDelete() const {
  return false;
}

// Private iv below here
MultifluidAlias::MfIVAliasFactory::MfIVAliasFactory(LevelData<MFBaseIVFAB>* a_mf, const int a_phase){
  m_mf    = a_mf;
  m_phase = a_phase;
}

BaseIVFAB<Real>* MultifluidAlias::MfIVAliasFactory::create(const Box& box, int ncomps, const DataIndex& a_datInd) const{
  return (*m_mf)[a_datInd].getPhasePtr(m_phase);
}

bool MultifluidAlias::MfIVAliasFactory::callDelete() const {
  return false;
}

#include <CD_NamespaceFooter.H>
