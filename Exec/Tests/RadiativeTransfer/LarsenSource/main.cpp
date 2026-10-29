/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Test of the Larsen domain boundary condition with a nonzero source, through the stationary EddingtonSP1
 * solver.
 * @details The domain is symmetric about the origin and has no embedded boundaries, the volumetric source is zero, and
 * every domain side has the same Larsen source. The exact solution is then symmetric, so the integrals of phi over the
 * low and the high half of the domain are equal in every direction.
 * @author Robert Marskar
 */

// Std includes
#include <cmath>
#include <limits>
#include <algorithm>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_EddingtonSP1.H>
#include <CD_RegularGeometry.H>
#include <CD_RadiativeTransferStepper.H>
#include <CD_ParallelOps.H>

using namespace ChomboDischarge;
using namespace Physics::RadiativeTransfer;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RegularGeometry());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<RadiativeTransferStepper<EddingtonSP1>>(
    new RadiativeTransferStepper<EddingtonSP1>());
  auto engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  Real maxAsymmetry;

  ParmParse pp("LarsenSource");
  pp.get("max_asymmetry", maxAsymmetry);

  std::string realm;

  ParmParse ppStepper("RadiativeTransferStepper");
  ppStepper.get("realm", realm);

  const EBAMRCellData& phi    = timestepper->getPhi();
  const RealVect       probLo = amr->getProbLo();

  // Integrals of phi over the low (x_dir < 0) and the high (x_dir > 0) half of the domain, over the cells that are not
  // covered by a finer level.
  Real integralLo[SpaceDim];
  Real integralHi[SpaceDim];

  for (int dir = 0; dir < SpaceDim; dir++) {
    integralLo[dir] = 0.0;
    integralHi[dir] = 0.0;
  }

  for (int lvl = 0; lvl <= amr->getFinestLevel(); lvl++) {
    const Real               dx   = amr->getDx()[lvl];
    const Real               dV   = std::pow(dx, SpaceDim);
    const DisjointBoxLayout& dbl  = amr->getGrids(realm)[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : integralLo[ : SpaceDim], integralHi[ : SpaceDim])
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&     din    = dit[mybox];
      const BaseFab<bool>& valid  = (*amr->getValidCells(realm)[lvl])[din];
      const FArrayBox&     phiFAB = (*phi[lvl])[din].getFArrayBox();

      for (BoxIterator bit(dbl[din]); bit.ok(); ++bit) {
        const IntVect iv = bit();

        if (valid(iv)) {
          const RealVect pos = probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * dx;

          for (int dir = 0; dir < SpaceDim; dir++) {
            if (pos[dir] < 0.0) {
              integralLo[dir] += phiFAB(iv, 0) * dV;
            }
            else {
              integralHi[dir] += phiFAB(iv, 0) * dV;
            }
          }
        }
      }
    }
  }

  bool passed = true;

  for (int dir = 0; dir < SpaceDim; dir++) {
    const Real lo = ParallelOps::sum(integralLo[dir]);
    const Real hi = ParallelOps::sum(integralHi[dir]);

    const Real asymmetry = std::abs(hi - lo) / std::max(std::abs(hi + lo), std::numeric_limits<Real>::min());

    pout() << "LarsenSource -- direction " << dir << ": integral of phi over the low half = " << lo
           << ", over the high half = " << hi << ", relative difference = " << asymmetry << endl;

    if (!(asymmetry <= maxAsymmetry)) {
      passed = false;
    }
  }

  if (!passed) {
    MayDay::Abort("LarsenSource -- the solution is not symmetric about the origin");
  }

  pout() << "LarsenSource -- all tests passed" << endl;

  ChomboDischarge::finalize();
}
