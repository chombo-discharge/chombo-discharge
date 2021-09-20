/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzSaturationChargeBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzSaturationChargeBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzSaturationChargeBCFactory.H>
#include <CD_MFHelmholtzSaturationChargeBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzSaturationChargeBCFactory::MFHelmholtzSaturationChargeBCFactory(const phase::which_phase a_phase){
  CH_TIME("MFHelmholtzSaturationChargeBCFactory::MFHelmholtzSaturationChargeBCFactory");

  m_phase = a_phase;
}

MFHelmholtzSaturationChargeBCFactory::~MFHelmholtzSaturationChargeBCFactory(){
  CH_TIME("MFHelmholtzSaturationChargeBCFactory::~MFHelmholtzSaturationChargeBCFactory");
}

RefCountedPtr<JumpBC> MFHelmholtzSaturationChargeBCFactory::create(const Location::Cell a_dataLocation,
								   const MFLevelGrid&   a_mflg,
								   const BcoefPtr&      a_Bcoef,
								   const Real           a_dx,
								   const int            a_order,
								   const int            a_weight,
								   const int            a_radius,
								   const int            a_ghostCF){

  
  return RefCountedPtr<JumpBC> (new MFHelmholtzSaturationChargeBC(m_phase, a_dataLocation, a_mflg, a_Bcoef, a_dx, a_order, a_weight, a_radius, a_ghostCF));
}

#include <CD_NamespaceFooter.H>
