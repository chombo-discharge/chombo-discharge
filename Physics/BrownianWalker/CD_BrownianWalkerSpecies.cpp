/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BrownianWalkerSpecies.cpp
  @brief  Implementation of CD_BrownianWalkerSpecies.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_Random.H>
#include <CD_ParticleOps.H>
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

BrownianWalkerSpecies::~BrownianWalkerSpecies() { CH_TIME("BrownianWalkerSpecies::~BrownianWalkerSpecies"); }

void
BrownianWalkerSpecies::drawInitParticles()
{
  CH_TIME("BrownianWalkerSpecies::drawInitParticles");

  // TLDR: We draw a bunch of particles from a Gaussian distribution center on m_blobCenter. If we use MPI then
  //       we should ensure that the ranks collectively draw the specified number of starting particles.

  // Nifty little lambda for drawing particle positions from a Gaussian distribution
  std::normal_distribution<Real> gaussianDistribution(0.0, m_blobRadius);
  auto                           randomGaussian = [center = this->m_blobCenter, &gaussianDistribution]() -> RealVect {
    const Real     len = Random::get(gaussianDistribution);
    const RealVect dir = Random::getDirection();

    return center + len * dir;
  };

  m_initialParticles.clear();
  ParticleOps::drawRandomParticles(m_initialParticles, m_numParticles, randomGaussian);

  // Set the particle mass to one.
  for (ListIterator<ItoParticle> lit(m_initialParticles); lit.ok(); ++lit) {
    lit().weight() = 1.0;
  }
}

#include <CD_NamespaceFooter.H>
