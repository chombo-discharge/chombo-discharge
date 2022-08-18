/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzJumpBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzJumpBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzJumpBCFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzJumpBCFactory::MFHelmholtzJumpBCFactory()
{
  CH_TIME("MFHelmholtzJumpBCFactory::MFHelmholtzJumpBCFactory");

  this->setDomainDropOrder(-1);
}

MFHelmholtzJumpBCFactory::~MFHelmholtzJumpBCFactory()
{
  CH_TIME("MFHelmholtzJumpBCFactory::~MFHelmholtzJumpBCFactory");
}

void
MFHelmholtzJumpBCFactory::setDomainDropOrder(const int a_domainSize)
{
  m_domainDropOrder = a_domainSize;
}

RefCountedPtr<MFHelmholtzJumpBC>
MFHelmholtzJumpBCFactory::create(const Location::Cell a_dataLocation,
                                 const MFLevelGrid&   a_mflg,
                                 const BcoefPtr&      a_Bcoef,
                                 const Real           a_dx,
                                 const int            a_order,
                                 const int            a_weight,
                                 const int            a_radius,
                                 const int            a_ghostCF)
{
  int order = a_order;

  // Drop order if we must
  for (int dir = 0; dir < SpaceDim; dir++) {
    if(a_mflg.getDomain().size()[dir] <= m_domainDropOrder) {
      order = 1;
    }
  }
  
  return RefCountedPtr<MFHelmholtzJumpBC>(
    new MFHelmholtzJumpBC(a_dataLocation, a_mflg, a_Bcoef, a_dx, order, a_weight, a_radius, a_ghostCF));
}

#include <CD_NamespaceFooter.H>
