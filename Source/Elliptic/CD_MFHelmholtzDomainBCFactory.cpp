/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  @file   CD_MFHelmholtzDomainBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzDomainBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDomainBCFactory::MFHelmholtzDomainBCFactory()
{
  CH_TIME("MFHelmholtzDomainBCFactory::MFHelmholtzDomainBCFactory()");
}

MFHelmholtzDomainBCFactory::~MFHelmholtzDomainBCFactory()
{
  CH_TIME("MFHelmholtzDomainBCFactory::~MFHelmholtzDomainBCFactory()");
}

#include <CD_NamespaceFooter.H>
