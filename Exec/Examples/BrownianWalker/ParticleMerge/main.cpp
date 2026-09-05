/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Example which measures what ItoSolver super-particle merging does to a uniform particle population
 * @author Robert Marskar
 */

// Our includes
#include <CD_Driver.H>
#include <CD_ItoSolver.H>
#include <CD_ParticleMergeStepper.H>
#include <CD_ParticleMergeTagger.H>

using namespace ChomboDischarge;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  {
    // No geometry: the default ComputationalGeometry has no electrodes and no dielectrics, so the domain
    // is all fluid and the merge is the only thing acting on the particles.
    auto compgeom = RefCountedPtr<ComputationalGeometry>(new ComputationalGeometry());
    auto amr      = RefCountedPtr<AmrMesh>(new AmrMesh());
    auto solver   = RefCountedPtr<ItoSolver>(new ItoSolver());

    auto tagger      = RefCountedPtr<ParticleMergeTagger>(new ParticleMergeTagger());
    auto timestepper = RefCountedPtr<TimeStepper>(new ParticleMergeStepper(solver));

    tagger->define(amr);

    auto engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

    engine->setupAndRun();
  }

  ChomboDischarge::finalize();
}
