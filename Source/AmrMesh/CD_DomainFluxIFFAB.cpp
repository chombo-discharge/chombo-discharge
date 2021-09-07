/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DomainFluxIFFAB.cpp
  @brief  Implementation of CD_DomainFluxIFFAB.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_DomainFluxIFFAB.H>
#include <CD_NamespaceHeader.H>
  
DomainFluxIFFAB::DomainFluxIFFAB(){
  CH_TIME("DomainFluxIFFAB::DomainFluxIFFAB");
  
  this->setDefaultValues();
}

DomainFluxIFFAB::DomainFluxIFFAB(const ProblemDomain& a_domain, const EBISBox& a_ebisbox, const Box& a_box, const int a_nComp){
  CH_TIME("DomainFluxIFFAB::DomainFluxIFFAB(ProblemDomain, EBISBox, Box, int)");
  
  this->setDefaultValues();
  this->define(a_domain, a_ebisbox, a_box, a_nComp);
}

DomainFluxIFFAB::~DomainFluxIFFAB(){
  CH_TIME("DomainFluxIFFAB::~DomainFluxIFFAB");
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if(m_fluxLo[dir] != nullptr) delete m_fluxLo[dir];
    if(m_fluxHi[dir] != nullptr) delete m_fluxHi[dir];
  }
}

const Box& DomainFluxIFFAB::box() const {
  return m_box;
}

int DomainFluxIFFAB::nComp() const {
  return m_nComp;
}

const ProblemDomain& DomainFluxIFFAB::getDomain() const {
  return m_domain;
}

const EBISBox& DomainFluxIFFAB::getEBISBox() const {
  return m_ebisbox;
}

void DomainFluxIFFAB::define(const DomainFluxIFFAB& a_copy){
  CH_TIME("DomainFluxIFFAB::define(DomainFluxIFFAB)");
  
  this->define(a_copy.getDomain(), a_copy.getEBISBox(), a_copy.box(), a_copy.nComp());
}

int DomainFluxIFFAB::size(const Box& R, const Interval& comps) const{
  CH_TIME("DomainFluxIFFAB::size");
  
  int retval = 0;
  for (int dir = 0; dir < SpaceDim; dir++){
    retval += m_fluxLo[dir]->size(R, comps);
    retval += m_fluxHi[dir]->size(R, comps);
  }
  
  return retval;
}

void DomainFluxIFFAB::linearOut(void* buf, const Box& R, const Interval& comps) const {
  CH_TIME("DomainFluxIFFAB::linearOut");
  
  unsigned char* charbuf = (unsigned char*) buf;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    m_fluxLo[dir]->linearOut(charbuf, R, comps);
    int dirsize = m_fluxLo[dir]->size(R, comps);
    charbuf += dirsize;
  }

  for (int dir = 0; dir < SpaceDim; dir++){
    m_fluxHi[dir]->linearOut(charbuf, R, comps);
    int dirsize = m_fluxHi[dir]->size(R, comps);
    charbuf += dirsize;
  }
}

void DomainFluxIFFAB::linearIn(void* buf, const Box& R, const Interval& comps){
  CH_TIME("DomainFluxIFFAB::linearIn");
  
  unsigned char* charbuf = (unsigned char*) buf;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    m_fluxLo[dir]->linearIn(charbuf, R, comps);
    int dirsize = m_fluxLo[dir]->size(R, comps);
    charbuf += dirsize;
  }

  for (int dir = 0; dir < SpaceDim; dir++){
    m_fluxHi[dir]->linearIn(charbuf, R, comps);
    int dirsize = m_fluxHi[dir]->size(R, comps);
    charbuf += dirsize;
  }
}

void DomainFluxIFFAB::define(const ProblemDomain& a_domain, const EBISBox& a_ebisbox, const Box& a_box, const int a_nComp){
  CH_TIME("DomainFluxIFFAB::define");
  
  m_domain    = a_domain;
  m_ebisbox   = a_ebisbox;
  m_nComp     = a_nComp;
  m_isDefined = true;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    Box shrinkbox = a_domain.domainBox();

    Box lobox = adjCellLo(shrinkbox, dir, 1);
    Box hibox = adjCellHi(shrinkbox, dir, 1);

    lobox.shift(dir, 1);
    hibox.shift(dir, -1);

    lobox &= a_box;
    hibox &= a_box;

    lobox &= a_domain;
    hibox &= a_domain;

    if(m_fluxLo[dir] != nullptr) {
      delete m_fluxLo[dir];
      m_fluxLo[dir] = nullptr;
    }
    if(m_fluxHi[dir] != nullptr) {
      delete m_fluxHi[dir];
      m_fluxHi[dir] = nullptr;
    }

    // ivsLo and ivsHi MAY CONTAIN COVERED CELLS!
    const IntVectSet ivsLo(lobox);
    const IntVectSet ivsHi(hibox);

#if 0
    // Get covered cells
    IntVectSet covered_lo = IntVectSet();
    IntVectSet covered_hi = IntVectSet();

    for (IVSIterator iter(ivsLo); iter.ok(); ++iter){
      const IntVect iv = iter();
      if(a_ebisbox.isCovered(iv)) covered_lo |= iv;
    }
    for (IVSIterator iter(ivsHi); iter.ok(); ++iter){
      const IntVect iv = iter();
      if(a_ebisbox.isCovered(iv)) covered_hi |= iv;
    }

    IntVectSet lo = ivsLo - covered_lo;
    IntVectSet hi = ivsHi - covered_hi;
#endif
    
#if 0
    if(!ivsLo.isEmpty()){
      std::cout << ivsLo << std::endl;
    }

    for (IVSIterator iter(ivsLo); iter.ok(); ++iter){
      if(a_ebisbox.isCovered(iter())){
	MayDay::Abort("stop");
      }
    }
#endif

    m_fluxLo[dir] = new BaseIFFAB<Real>(ivsLo, a_ebisbox.getEBGraph(), dir, a_nComp);
    m_fluxHi[dir] = new BaseIFFAB<Real>(ivsHi, a_ebisbox.getEBGraph(), dir, a_nComp);
  }
}

void DomainFluxIFFAB::setDefaultValues(){
  CH_TIME("DomainFluxIFFAB::setDefaultValues");
  
  m_isDefined = false;
  for (int dir = 0; dir < SpaceDim; dir++){
    m_fluxLo[dir] = nullptr;
    m_fluxHi[dir] = nullptr;
  }

  m_box    = Box();
  m_domain = ProblemDomain();
  m_nComp  = 0;
}

void DomainFluxIFFAB::clear(){
  CH_TIME("DomainFluxIFFAB::clear");
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if(m_fluxLo[dir] != nullptr){
      delete m_fluxLo[dir];
      m_fluxLo[dir] = nullptr;
    }
    if(m_fluxHi[dir] != nullptr){
      delete m_fluxHi[dir];
      m_fluxHi[dir] = nullptr;
    }
  }

  this->setDefaultValues();
}

void DomainFluxIFFAB::copy(const Box&             Rfrom,
			   const Interval&        Cdest,
			   const Box&             Rto,
			   const DomainFluxIFFAB& src,
			   const Interval&        Csrc){
  CH_TIME("DomainFluxIFFAB::copy");
  
  for (int dir = 0; dir < SpaceDim; dir++){
    const BaseIFFAB<Real>& src_lo = src(dir, Side::Lo);
    const BaseIFFAB<Real>& src_hi = src(dir, Side::Hi);

    m_fluxLo[dir]->copy(Rfrom, Cdest, Rto, src_lo, Csrc);
    m_fluxHi[dir]->copy(Rfrom, Cdest, Rto, src_hi, Csrc);
  }
}

BaseIFFAB<Real>& DomainFluxIFFAB::operator()(const int a_dir, const Side::LoHiSide a_side){
  CH_assert(m_isDefined);

  BaseIFFAB<Real>* ptr;
  if(a_side == Side::Lo){
    ptr = m_fluxLo[a_dir];
  }
  else if(a_side == Side::Hi){
    ptr = m_fluxHi[a_dir];
  }

  return *ptr;
}

const BaseIFFAB<Real>& DomainFluxIFFAB::operator()(const int a_dir, const Side::LoHiSide a_side) const {
  CH_assert(m_isDefined);

  BaseIFFAB<Real>* ptr;
  if(a_side == Side::Lo){
    ptr = m_fluxLo[a_dir];
  }
  else if(a_side == Side::Hi){
    ptr = m_fluxHi[a_dir];
  }

  return *ptr;
}

#include <CD_NamespaceFooter.H>
