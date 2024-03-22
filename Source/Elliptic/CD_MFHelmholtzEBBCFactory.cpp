/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzEBBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzEBBCFactory::MFHelmholtzEBBCFactory()
{
  CH_TIME("MFHelmholtzEBBCFactory::MFHelmholtzEBBCFactory()");
}

MFHelmholtzEBBCFactory::~MFHelmholtzEBBCFactory()
{
  CH_TIME("MFHelmholtzEBBCFactory::~MFHelmholtzEBBCFactory()");
}

#include <CD_NamespaceFooter.H>
