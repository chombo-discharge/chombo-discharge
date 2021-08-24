/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzDomainBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzDomainBC::m_comp;
constexpr int EBHelmholtzDomainBC::m_nComp;

EBHelmholtzDomainBC::EBHelmholtzDomainBC(){
  CH_TIME("EBHelmholtzDomainBC::EBHelmholtzDomainBC()");
}

EBHelmholtzDomainBC::~EBHelmholtzDomainBC(){
  CH_TIME("EBHelmholtzDomainBC::~EBHelmholtzDomainBC()");
}

void EBHelmholtzDomainBC::define(const Location::Cell                        a_dataLocation,
				 const EBLevelGrid&                          a_eblg,
				 const RefCountedPtr<LevelData<EBFluxFAB> >& a_Bcoef,
				 const RealVect&                             a_probLo,
				 const Real                                  a_dx){
  CH_TIME("EBHelmholtzDomainBC::define(Location::Cell, EBLevelGrid, RefCountedPtr<LD<EBFluxFAB> >, RealVect, Real)");

  CH_assert(a_dx > 0.0);
  
  m_dataLocation = a_dataLocation;
  m_eblg         = a_eblg;
  m_Bcoef        = a_Bcoef;
  m_probLo       = a_probLo;
  m_dx           = a_dx;
}

#include <CD_NamespaceFooter.H>
