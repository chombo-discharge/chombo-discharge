/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFBaseIVFAB.cpp
  @brief  Implementation of CD_MFBaseIVFAB.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFBaseIVFAB.H>
#include <CD_NamespaceHeader.H>

MFBaseIVFAB::MFBaseIVFAB(){
}

MFBaseIVFAB::MFBaseIVFAB(const Vector<IntVectSet>& a_regions,
			 const Vector<EBGraph>&    a_phase_graphs,
			 const Vector<int>&        a_ncomp){
  CH_TIME("MFBaseIVFAB::MFBaseIVFAB");

  const int numPhases = a_regions.size();
    
  CH_assert(a_phase_graphs.size() == numPhases);
  CH_assert(a_ncomp.size() == numPhases);
    
  m_phase.resize(numPhases, NULL);
  for (int i = 0; i < numPhases; i++){
    m_phase[i] = new BaseIVFAB<Real>(a_regions[i], a_phase_graphs[i], a_ncomp[i]);
  }
}

MFBaseIVFAB::~MFBaseIVFAB(){
  for (int i = 0; i < m_phase.size(); i++){
    delete m_phase[i];
  }
}


BaseIVFAB<Real>& MFBaseIVFAB::getIVFAB(const int a_phase){
  return *m_phase[a_phase];
}

const BaseIVFAB<Real>& MFBaseIVFAB::getIVFAB(const int a_phase) const {
  return *m_phase[a_phase];
}

BaseIVFAB<Real>* MFBaseIVFAB::getPhasePtr(int a_phase){
  return m_phase[a_phase];
}

int MFBaseIVFAB::numPhases(){
  return m_phase.size();
}

void MFBaseIVFAB::setVal(Real a_value){
  for (int i = 0; i < m_phase.size(); i++){
    m_phase[i]->setVal(a_value);
  }
}

void MFBaseIVFAB::copy(const Box& a_from_box,
		       const Interval& a_dst_interv,
		       const Box& a_to_box,
		       const MFBaseIVFAB& a_src,
		       const Interval& a_src_interv){
  for (int i = 0; i < m_phase.size(); i++){
    m_phase[i]->copy(a_from_box, a_dst_interv, a_to_box, *(a_src.m_phase[i]), a_src_interv);
  }
}

int MFBaseIVFAB::preAllocatable(){
  return 1;
}

int MFBaseIVFAB::size(const Box& R, const Interval& comps) const {
  int size = m_phase.size()*sizeof(int);
  for (int i = 0; i < m_phase.size(); ++i){
    size += m_phase[i]->size(R, comps);
  }

  return size;
}

void MFBaseIVFAB::linearOut(void* buf, const Box& R, const Interval& comps) const {
  int* buffer = (int*)buf;
  for (int i=0; i<m_phase.size(); ++i)
    {
      *buffer = m_phase[i]->size(R, comps);
      ++buffer;
    }
  int* size = (int*)buf;
  unsigned char* ebbuffer = (unsigned char*)buffer;
  for (int i=0; i<m_phase.size(); ++i)
    {
      m_phase[i]->linearOut(ebbuffer, R, comps);
      ebbuffer+= size[i];
    }
}

void MFBaseIVFAB::linearIn(void* buf, const Box& R, const Interval& comps){
  int* size = (int*) buf;

  for (int i = 0; i < m_phase.size(); ++i){
    CH_assert(size[i] == m_phase[i]->size(R, comps));
  }

  unsigned char* ebbuffer  = (unsigned char*)(size+m_phase.size());
  for (int i=0; i<m_phase.size(); ++i){
    m_phase[i]->linearIn(ebbuffer, R, comps);
    ebbuffer += size[i];
  }
}

MFBaseIVFABFactory::MFBaseIVFABFactory(Vector<EBISLayout>& a_ebisl, const Vector<int>& a_ncomp){
  CH_TIME("MFBaseIVFABFactory::MFBaseIVFABFactory");
  this->define(a_ebisl, a_ncomp);
}

MFBaseIVFABFactory::~MFBaseIVFABFactory(){
  CH_TIME("MFBaseIVFABFactory::~MFBaseIVFABFactory");
}

void MFBaseIVFABFactory::define(Vector<EBISLayout>& a_ebisl, const Vector<int>& a_ncomp){
  CH_TIME("MFBaseIVFABFactory::define");
  CH_assert(a_ebisl.size() == a_ncomp.size());
  m_ebisl = a_ebisl;
  m_ncomp = a_ncomp;
}

MFBaseIVFAB* MFBaseIVFABFactory::create(const Box& a_box, int a_ignored_argument, const DataIndex& a_dit) const {
  CH_TIME("MFBaseIVFABFactory::create");

  const int numPhases = m_ebisl.size();
    
  Vector<IntVectSet> ivs(numPhases);
  Vector<EBGraph> ebgraph(numPhases);

  for (int i = 0; i < numPhases; ++i){
    ivs[i]     = m_ebisl[i][a_dit].getIrregIVS(a_box);
    ebgraph[i] = m_ebisl[i][a_dit].getEBGraph();
  }

  return new MFBaseIVFAB(ivs, ebgraph, m_ncomp);
}

#include <CD_NamespaceFooter.H>
