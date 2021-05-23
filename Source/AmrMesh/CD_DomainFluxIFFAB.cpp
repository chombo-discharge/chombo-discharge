/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DomainFluxIFFAB.cpp
  @brief  Implementation of CD_DomainFluxIFFAB.H
  @author Robert Marskar
*/

// Our includes
#include <CD_DomainFluxIFFAB.H>
#include <CD_NamespaceHeader.H>
  
DomainFluxIFFAB::DomainFluxIFFAB(){
  setDefaultValues();
}

DomainFluxIFFAB::DomainFluxIFFAB(const ProblemDomain& a_domain, const EBISBox& a_ebisbox, const Box& a_box, const int a_ncomp){
  setDefaultValues();
  define(a_domain, a_ebisbox, a_box, a_ncomp);
}

DomainFluxIFFAB::~DomainFluxIFFAB(){
  for (int dir = 0; dir < SpaceDim; dir++){
    if(m_flux_lo[dir] != NULL) delete m_flux_lo[dir];
    if(m_flux_hi[dir] != NULL) delete m_flux_hi[dir];
  }
}

int DomainFluxIFFAB::size(const Box& R, const Interval& comps) const{
  int retval = 0;
  for (int dir = 0; dir < SpaceDim; dir++){
    retval += m_flux_lo[dir]->size(R, comps);
    retval += m_flux_hi[dir]->size(R, comps);
  }
  return retval;
}

void DomainFluxIFFAB::linearOut(void* buf, const Box& R, const Interval& comps) const {
  unsigned char* charbuf = (unsigned char*) buf;
  for (int dir = 0; dir < SpaceDim; dir++){
    m_flux_lo[dir]->linearOut(charbuf, R, comps);
    int dirsize = m_flux_lo[dir]->size(R, comps);
    charbuf += dirsize;
  }

  for (int dir = 0; dir < SpaceDim; dir++){
    m_flux_hi[dir]->linearOut(charbuf, R, comps);
    int dirsize = m_flux_hi[dir]->size(R, comps);
    charbuf += dirsize;
  }
}

void DomainFluxIFFAB::linearIn(void* buf, const Box& R, const Interval& comps){
  unsigned char* charbuf = (unsigned char*) buf;
  for (int dir = 0; dir < SpaceDim; dir++){
    m_flux_lo[dir]->linearIn(charbuf, R, comps);
    int dirsize = m_flux_lo[dir]->size(R, comps);
    charbuf += dirsize;
  }

  for (int dir = 0; dir < SpaceDim; dir++){
    m_flux_hi[dir]->linearIn(charbuf, R, comps);
    int dirsize = m_flux_hi[dir]->size(R, comps);
    charbuf += dirsize;
  }
}

void DomainFluxIFFAB::define(const ProblemDomain& a_domain, const EBISBox& a_ebisbox, const Box& a_box, const int a_ncomp){
  m_domain    = a_domain;
  m_ebisbox   = a_ebisbox;
  m_nComp     = a_ncomp;
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

    if(m_flux_lo[dir] != NULL) {
      delete m_flux_lo[dir];
      m_flux_lo[dir] = NULL;
    }
    if(m_flux_hi[dir] != NULL) {
      delete m_flux_hi[dir];
      m_flux_hi[dir] = NULL;
    }

    // ivs_lo and ivs_hi MAY CONTAIN COVERED CELLS!
    const IntVectSet ivs_lo(lobox);
    const IntVectSet ivs_hi(hibox);

#if 0
    // Get covered cells
    IntVectSet covered_lo = IntVectSet();
    IntVectSet covered_hi = IntVectSet();

    for (IVSIterator iter(ivs_lo); iter.ok(); ++iter){
      const IntVect iv = iter();
      if(a_ebisbox.isCovered(iv)) covered_lo |= iv;
    }
    for (IVSIterator iter(ivs_hi); iter.ok(); ++iter){
      const IntVect iv = iter();
      if(a_ebisbox.isCovered(iv)) covered_hi |= iv;
    }

    IntVectSet lo = ivs_lo - covered_lo;
    IntVectSet hi = ivs_hi - covered_hi;
#endif
    
#if 0
    if(!ivs_lo.isEmpty()){
      std::cout << ivs_lo << std::endl;
    }

    for (IVSIterator iter(ivs_lo); iter.ok(); ++iter){
      if(a_ebisbox.isCovered(iter())){
	MayDay::Abort("stop");
      }
    }
#endif

    m_flux_lo[dir] = new BaseIFFAB<Real>(ivs_lo, a_ebisbox.getEBGraph(), dir, a_ncomp);
    m_flux_hi[dir] = new BaseIFFAB<Real>(ivs_hi, a_ebisbox.getEBGraph(), dir, a_ncomp);
  }
}

void DomainFluxIFFAB::setDefaultValues(){
  m_isDefined = false;
  for (int dir = 0; dir < SpaceDim; dir++){
    m_flux_lo[dir] = NULL;
    m_flux_hi[dir] = NULL;
  }

  m_box    = Box();
  m_domain = ProblemDomain();
  m_nComp  = 0;
}

void DomainFluxIFFAB::clear(){
  for (int dir = 0; dir < SpaceDim; dir++){
    if(m_flux_lo[dir] != NULL){
      delete m_flux_lo[dir];
      m_flux_lo[dir] = NULL;
    }
    if(m_flux_hi[dir] != NULL){
      delete m_flux_hi[dir];
      m_flux_hi[dir] = NULL;
    }
  }

  setDefaultValues();
}

void DomainFluxIFFAB::copy(const Box& Rfrom,
			   const Interval& Cdest,
			   const Box& Rto,
			   const DomainFluxIFFAB& src,
			   const Interval& Csrc){
  for (int dir = 0; dir < SpaceDim; dir++){
    const BaseIFFAB<Real>& src_lo = src(dir, Side::Lo);
    const BaseIFFAB<Real>& src_hi = src(dir, Side::Hi);

    m_flux_lo[dir]->copy(Rfrom, Cdest, Rto, src_lo, Csrc);
    m_flux_hi[dir]->copy(Rfrom, Cdest, Rto, src_hi, Csrc);
  }
}

BaseIFFAB<Real>& DomainFluxIFFAB::operator()(const int a_dir, const Side::LoHiSide a_side){
  CH_assert(m_isDefined);

  BaseIFFAB<Real>* ptr;
  if(a_side == Side::Lo){
    ptr = m_flux_lo[a_dir];
  }
  else if(a_side == Side::Hi){
    ptr = m_flux_hi[a_dir];
  }

  return *ptr;
}

const BaseIFFAB<Real>& DomainFluxIFFAB::operator()(const int a_dir, const Side::LoHiSide a_side) const {
  CH_assert(m_isDefined);

  BaseIFFAB<Real>* ptr;
  if(a_side == Side::Lo){
    ptr = m_flux_lo[a_dir];
  }
  else if(a_side == Side::Hi){
    ptr = m_flux_hi[a_dir];
  }

  return *ptr;
}

#include <CD_NamespaceFooter.H>
