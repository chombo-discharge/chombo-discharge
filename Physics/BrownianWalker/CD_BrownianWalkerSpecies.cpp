/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_BrownianWalkerSpecies.cpp
   @brief  Implementation of CD_BrownianWalkerSpecies.H
   @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_ParticleManagement.H>
#include <CD_BrownianWalkerSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::BrownianWalker;

BrownianWalkerSpecies::BrownianWalkerSpecies()
{
  CH_TIME("BrownianWalkerSpecies::BrownianWalkerSpecies");

  m_name         = "scalar species";
  m_chargeNumber = 0;

  // Parse input options.
  ParmParse    pp("BrownianWalker");
  Vector<Real> v;

  int seed;

  pp.get("seed", seed);
  pp.get("diffusion", m_isDiffusive);
  pp.get("advection", m_isMobile);
  pp.get("num_particles", m_numParticles);
  pp.get("blob_radius", m_blobRadius);

  pp.getarr("blob_center", v, 0, SpaceDim);
  m_blobCenter = RealVect(D_DECL(v[0], v[1], v[2]));

  // Draw initial particles
  this->drawInitParticles();
}

BrownianWalkerSpecies::~BrownianWalkerSpecies()
{
  CH_TIME("BrownianWalkerSpecies::~BrownianWalkerSpecies");
}

void
BrownianWalkerSpecies::drawInitParticles()
{
  CH_TIME("BrownianWalkerSpecies::drawInitParticles");

  // Draw Gaussian particles and set weights to 1
  ParticleManagement::drawGaussianParticles(m_initialParticles, m_numParticles, m_blobCenter, m_blobRadius);

  double* w = m_initialParticles.weightColumn();

  ParticleLoops::loop(m_initialParticles, [&](std::size_t i) {
    w[i] = 1.0;
  });
}

#include <CD_NamespaceFooter.H>
