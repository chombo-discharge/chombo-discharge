/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBRedistribution.cpp
 * @brief  Implementation of CD_EBRedistribution.H
 * @author Robert Marskar
 */

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBRedistribution.H>
#include <CD_NamespaceHeader.H>

EBRedistribution::EBRedistribution() noexcept
{
  CH_TIME("EBRedistribution::EBRedistribution");
}

EBRedistribution::~EBRedistribution() noexcept
{
  CH_TIME("EBRedistribution::~EBRedistribution");
}

#include <CD_NamespaceFooter.H>
