/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Advection radially outwards from a sphere, with refinement boundaries that cross the sphere.
 * @author Robert Marskar
 */

// Std includes
#include <cmath>
#include <string>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_CdrCTU.H>
#include <CD_RoughSphere.H>
#include <CD_RefinementBoxTagger.H>
#include <CD_AdvectionDiffusionStepper.H>

using namespace ChomboDischarge;
using namespace Physics::AdvectionDiffusion;

/**
 * @brief Advection radially outwards from a sphere at the origin, with a static mesh whose refinement boundaries cross
 * the sphere.
 * @details The velocity is RadialOutflow.speed times the unit vector away from the origin. RadialOutflow.initial_data
 * selects the initial state:
 *
 *   - 'blob': a Gaussian exp(-|x - c|^2/r^2) with center c = RadialOutflow.blob_center and radius r =
 *     RadialOutflow.blob_radius,
 *   - 'half_space': 1 for x > 0 and 0 elsewhere,
 *   - 'uniform': 1 everywhere.
 *
 * The sphere is RoughSphere without noise, and the refinement is given by RefinementBoxTagger.
 */
int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  ParmParse pp("RadialOutflow");

  std::string  initialData;
  Real         speed      = 0.0;
  Real         blobRadius = 0.0;
  Vector<Real> v(SpaceDim, 0.0);

  pp.get("initial_data", initialData);
  pp.get("speed", speed);
  pp.get("blob_radius", blobRadius);
  pp.getarr("blob_center", v, 0, SpaceDim);

  const RealVect blobCenter = RealVect(D_DECL(v[0], v[1], v[2]));

  std::function<Real(const RealVect&)> phi0;

  if (initialData == "blob") {
    phi0 = [blobCenter, blobRadius](const RealVect& a_pos) -> Real {
      const RealVect d = a_pos - blobCenter;

      return std::exp(-d.dotProduct(d) / (blobRadius * blobRadius));
    };
  }
  else if (initialData == "half_space") {
    phi0 = [](const RealVect& a_pos) -> Real {
      return (a_pos[0] > 0.0) ? 1.0 : 0.0;
    };
  }
  else if (initialData == "uniform") {
    phi0 = [](const RealVect& /*a_pos*/) -> Real {
      return 1.0;
    };
  }
  else {
    MayDay::Abort("RadialOutflow -- 'initial_data' must be 'blob', 'half_space' or 'uniform'");
  }

  auto velocity = [speed](const RealVect& a_pos) -> RealVect {
    const Real r = a_pos.vectorLength();

    return (r > 0.0) ? speed * a_pos / r : RealVect::Zero;
  };

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RoughSphere());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto solver      = RefCountedPtr<CdrSolver>(new CdrCTU());
  auto timestepper = RefCountedPtr<AdvectionDiffusionStepper>(new AdvectionDiffusionStepper(solver));
  auto tagger      = RefCountedPtr<CellTagger>(new RefinementBoxTagger(amr, Realm::Primal));
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  timestepper->setInitialData(phi0);
  timestepper->setVelocity(velocity);

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
