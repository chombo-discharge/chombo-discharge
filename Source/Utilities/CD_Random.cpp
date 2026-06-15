/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Random.cpp
  @brief  Implementation of CD_Random.H
  @author Robert Marskar
*/

// Our includes
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

thread_local std::mt19937_64                      Random::s_rng       = std::mt19937_64(0);
thread_local std::uniform_real_distribution<Real> Random::s_uniform01 = std::uniform_real_distribution<Real>(0.0, 1.0);
thread_local std::uniform_real_distribution<Real> Random::s_uniform11 = std::uniform_real_distribution<Real>(-1.0, 1.0);
thread_local std::normal_distribution<Real>       Random::s_normal01  = std::normal_distribution<Real>(0.0, 1.0);

bool Random::s_seeded = false;

//std::once_flag once = std::once_flag();

#include <CD_NamespaceFooter.H>
