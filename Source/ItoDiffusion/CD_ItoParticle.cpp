/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoParticle.cpp
  @brief  Implementation of CD_ItoParticle.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ItoParticle.H>
#include <CD_NamespaceHeader.H>

std::vector<std::string> ItoParticle::s_realVariables = {"weight", "mobility", "diffusion", "energy", "scratch"};
std::vector<std::string> ItoParticle::s_vectVariables = {"old", "v", "scratch"};

#include <CD_NamespaceFooter.H>
