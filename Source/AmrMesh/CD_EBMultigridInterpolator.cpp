/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBMultigridInterpolator.cpp
  @brief  Implementation of CD_EBMultigridInterpolator.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

EBMultigridInterpolator::EBMultigridInterpolator()
{
  CH_TIME("EBMultigridInterpolator::EBMultigridInterpolator");
}

EBMultigridInterpolator::~EBMultigridInterpolator()
{
  CH_TIME("EBMultigridInterpolator::~EBMultigridInterpolator");
}

#include <CD_NamespaceFooter.H>
