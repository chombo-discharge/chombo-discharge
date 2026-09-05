/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMergeStepper.H
 * @brief  Declaration of a Brownian-walker stepper that reseeds and reports on the merge every step
 * @author Robert Marskar
 */

#ifndef CD_PARTICLEMERGESTEPPER_H
#define CD_PARTICLEMERGESTEPPER_H

// Our includes
#include <CD_BrownianWalkerStepper.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Brownian-walker stepper that seeds a fresh uniform population before every advance and can
 * report what the super-particle merge did to it.
 * @details BrownianWalkerStepper::advance already transports the particles and then merges them down to
 * ItoSolver.particles_per_cell, which is the operation under study. This class adds the two things a
 * measurement needs around it: a known population going in, and a report coming out.
 *
 * Before each advance it seeds ParticleMergeStepper.new_particles_per_cell fresh unit-weight particles
 * per cell, uniformly distributed, on top of whatever survived the previous merge. Repeating that is
 * what exposes an error that compounds, which is the failure mode that matters for an operator applied
 * once per step for thousands of steps. The same population is also what initialData() lays down, so
 * the particle time step is bounded from the first step.
 *
 * Every option is read under the ParticleMergeStepper prefix, including the ones inherited from
 * BrownianWalkerStepper, so the stepper has one prefix and it matches the class name.
 *
 * Run it with ParticleMergeStepper.diffusion = false so the only motion is the advective field. Under a
 * constant velocity the population then translates rigidly, which leaves it uniform and keeps the merge
 * the only thing acting on it.
 */
class ParticleMergeStepper : public Physics::BrownianWalker::BrownianWalkerStepper
{
public:
  /**
   * @brief Default constructor. Not used -- the solver has to be supplied.
   */
  ParticleMergeStepper() = delete;

  /**
   * @brief Full constructor.
   * @param[in] a_solver Ito solver to advance and merge.
   */
  ParticleMergeStepper(RefCountedPtr<ItoSolver>& a_solver);

  /**
   * @brief Destructor.
   */
  virtual ~ParticleMergeStepper();

  /**
   * @brief Lay down the initial population.
   * @details The same uniform seeding each advance uses, rather than the base class's blob. Seeding here
   * rather than leaving the container empty is what bounds the first particle time step.
   */
  void
  initialData() override;

  /**
   * @brief Seed a fresh population, advance and merge, then report if asked to.
   * @param[in] a_dt Time step.
   * @return The time step actually taken.
   */
  Real
  advance(const Real a_dt) override;

protected:
  /**
   * @brief Particles seeded per cell before every advance.
   * @details ParticleMergeStepper.new_particles_per_cell.
   */
  int m_newParticlesPerCell;

  /**
   * @brief Whether to print the per-level merge statistics to pout on every step.
   * @details ParticleMergeStepper.print_stats.
   */
  bool m_printStats;

  /**
   * @brief Seed a fresh uniform population of unit-weight particles on every valid cell.
   * @details Cells covered by a finer level are skipped: particles seeded there would be discarded on
   * the next remap, so seeding them would make the population differ from the one the merge sees.
   */
  void
  seedParticles() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
