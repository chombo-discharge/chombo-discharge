/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
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
