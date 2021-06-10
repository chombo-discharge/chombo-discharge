/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzEBBC::EBHelmholtzEBBC(){

}

EBHelmholtzEBBC::~EBHelmholtzEBBC(){

}

void EBHelmholtzEBBC::applyEBFlux(VoFIterator&       a_vofit,
				  EBCellFAB&         a_Lphi,
				  const EBCellFAB&   a_phi,
				  const RealVect&    a_probLo,
				  const RealVect&    a_dx,
				  const Real&        a_factor,
				  const bool&        a_useHomogeneous) const {
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit){
    this->applyEBFlux(a_Lphi, a_phi, a_vofit(), a_probLo, a_dx, a_factor, a_useHomogeneous);
  }
}

void EBHelmholtzEBBC::define(const EBLevelGrid& a_eblg, const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_Bcoef, const Real a_factor){
  m_Bcoef = a_Bcoef;
  
  this->define(a_eblg, a_factor);
}

const LayoutData<BaseIVFAB<VoFStencil> >& EBHelmholtzEBBC::getFluxStencil() const{
  return m_fluxStencil;
}



#include <CD_NamespaceFooter.H>
