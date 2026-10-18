/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Test of FieldSolverAMG on a single level without embedded boundaries.
 * @author Robert Marskar
 */

// Std includes
#include <cmath>
#include <functional>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_FieldSolverGMG.H>
#include <CD_FieldSolverAMG.H>
#include <CD_RegularGeometry.H>
#include <CD_FieldStepper.H>
#include <CD_ParallelOps.H>
#include <CD_Units.H>
#include <CD_DataOps.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

/**
 * @brief Wave number of the manufactured solution.
 */
static constexpr Real s_k = 2.0;

/**
 * @brief Phase shift of the manufactured solution in a coordinate direction.
 * @param[in] a_dir Coordinate direction
 * @return Phase shift.
 */
static Real
phaseShift(const int a_dir)
{
  return 0.5 + 0.25 * a_dir;
}

/**
 * @brief The manufactured solution, phi = prod_d sin(k*x_d + c_d).
 * @param[in] a_pos Position
 * @return phi at a_pos.
 */
static Real
exactPhi(const RealVect& a_pos)
{
  Real phi = 1.0;

  for (int dir = 0; dir < SpaceDim; dir++) {
    phi *= std::sin(s_k * a_pos[dir] + phaseShift(dir));
  }

  return phi;
}

/**
 * @brief Derivative of the manufactured solution in a coordinate direction.
 * @param[in] a_pos Position
 * @param[in] a_dir Coordinate direction
 * @return d(phi)/d(x_dir) at a_pos.
 */
static Real
exactGradient(const RealVect& a_pos, const int a_dir)
{
  Real grad = 1.0;

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (dir == a_dir) {
      grad *= s_k * std::cos(s_k * a_pos[dir] + phaseShift(dir));
    }
    else {
      grad *= std::sin(s_k * a_pos[dir] + phaseShift(dir));
    }
  }

  return grad;
}

/**
 * @brief Set up a field solver on an existing mesh, as FieldStepper does.
 * @param[in] a_amr      AMR mesh
 * @param[in] a_compgeom Computational geometry
 * @param[in] a_realm    Realm
 * @param[in] a_className Class name of the solver, which prefixes its options
 * @return The solver.
 */
template <class T>
static RefCountedPtr<FieldSolver>
makeSolver(const RefCountedPtr<AmrMesh>&               a_amr,
           const RefCountedPtr<ComputationalGeometry>& a_compgeom,
           const std::string&                          a_realm,
           const std::string&                          a_className)
{
  auto solver = RefCountedPtr<FieldSolver>(new T());

  solver->parseOptions();
  solver->setAmr(a_amr);
  solver->setComputationalGeometry(a_compgeom);
  solver->setVoltage([](const Real a_time) -> Real {
    return 1.0;
  });
  solver->setRealm(a_realm);
  solver->setTime(0, 0.0, 0.0);
  solver->allocate();

  // The domain sides are 'dirichlet_custom' or 'neumann_custom' in the input file. Give both the exact solution:
  // the value for Dirichlet, and the inward normal derivative for Neumann.
  ParmParse pp(a_className.c_str());

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const Side::LoHiSide side = sit();
      const std::string    key  = std::string("bc.") + "xyz"[dir] + "." + (side == Side::Lo ? "lo" : "hi");

      std::string bcType;
      pp.get(key.c_str(), bcType);

      if (bcType == "dirichlet_custom") {
        solver->setDomainSideBcFunction(dir, side, [](const RealVect a_pos, const Real a_time) -> Real {
          return exactPhi(a_pos);
        });
      }
      else if (bcType == "neumann_custom") {
        solver->setDomainSideBcFunction(dir, side, [dir, side](const RealVect a_pos, const Real a_time) -> Real {
          return -sign(side) * exactGradient(a_pos, dir);
        });
      }
      else {
        MayDay::Abort("PetscRegular -- domain sides must be 'dirichlet_custom' or 'neumann_custom'");
      }
    }
  }

  // -div(grad(phi)) = rho/eps0, and the Laplacian of phi is -SpaceDim * k^2 * phi.
  solver->setRho([](const RealVect a_pos) -> Real {
    return Units::eps0 * SpaceDim * s_k * s_k * exactPhi(a_pos);
  });

  return solver;
}

/**
 * @brief Get the max-norm difference between two potentials, or between a potential and the exact solution.
 * @param[in] a_amr   AMR mesh
 * @param[in] a_realm Realm
 * @param[in] a_phi   Potential
 * @param[in] a_other Other potential. If null, the exact solution is used.
 * @return Max-norm difference over the gas-phase cells on level 0.
 */
static Real
maxDifference(const AmrMesh&       a_amr,
              const std::string&   a_realm,
              const MFAMRCellData& a_phi,
              const MFAMRCellData* a_other)
{
  constexpr int lvl = 0;

  const DisjointBoxLayout& dbl    = a_amr.getGrids(a_realm)[lvl];
  const Real               dx     = a_amr.getDx()[lvl];
  const RealVect           probLo = a_amr.getProbLo();

  Real maxDiff = 0.0;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const FArrayBox& phi = (*a_phi[lvl])[dit()].getPhase(phase::gas).getFArrayBox();

    for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit) {
      const IntVect  iv    = bit();
      const RealVect pos   = probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * dx;
      const Real     other = (a_other == nullptr) ? exactPhi(pos)
                                                  : (*(*a_other)[lvl])[dit()].getPhase(phase::gas).getFArrayBox()(iv, 0);

      maxDiff = std::max(maxDiff, std::abs(phi(iv, 0) - other));
    }
  }

  return ParallelOps::max(maxDiff);
}

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

#ifdef CH_USE_PETSC
  // Scoped so that the PETSc objects are destroyed before ChomboDischarge::finalize finalizes PETSc and MPI.
  {
    // The driver sets up the mesh. Its FieldStepper uses FieldSolverGMG, so the operators GMG needs are registered.
    auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RegularGeometry());
    auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
    auto tagger      = RefCountedPtr<CellTagger>(nullptr);
    auto timestepper = RefCountedPtr<FieldStepper<FieldSolverGMG>>(new FieldStepper<FieldSolverGMG>());
    auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

    engine->setupAndRun();

    std::string realm;
    Real        maxErrorConstant;
    Real        maxAllowedDifference;

    ParmParse ppStepper("FieldStepper");
    ppStepper.get("realm", realm);

    ParmParse pp("PetscRegular");
    pp.get("max_error_constant", maxErrorConstant);
    pp.get("max_difference", maxAllowedDifference);

    RefCountedPtr<FieldSolver> amg = makeSolver<FieldSolverAMG>(amr, compgeom, realm, "FieldSolverAMG");
    RefCountedPtr<FieldSolver> gmg = makeSolver<FieldSolverGMG>(amr, compgeom, realm, "FieldSolverGMG");

    // There are no embedded boundaries, so the surface charge is zero.
    EBAMRIVData sigma;
    amr->allocate(sigma, realm, phase::gas, 1);
    DataOps::setValue(sigma, 0.0);

    const bool amgConverged = amg->solve(amg->getPotential(), amg->getRho(), sigma, true);
    const bool gmgConverged = gmg->solve(gmg->getPotential(), gmg->getRho(), sigma, true);

    const Real dx       = amr->getDx()[0];
    const Real errorAMG = maxDifference(*amr, realm, amg->getPotential(), nullptr);
    const Real errorGMG = maxDifference(*amr, realm, gmg->getPotential(), nullptr);
    const Real diff     = maxDifference(*amr, realm, amg->getPotential(), &gmg->getPotential());

    pout() << "PetscRegular -- dx = " << dx << ", max error AMG = " << errorAMG << ", max error GMG = " << errorGMG
           << ", max |AMG - GMG| = " << diff << ", error/dx^2 = " << errorAMG / (dx * dx) << endl;

    if (!amgConverged) {
      MayDay::Abort("PetscRegular -- FieldSolverAMG did not converge");
    }
    if (!gmgConverged) {
      MayDay::Abort("PetscRegular -- FieldSolverGMG did not converge");
    }
    if (diff > maxAllowedDifference) {
      MayDay::Abort("PetscRegular -- FieldSolverAMG and FieldSolverGMG differ by more than the solver tolerance");
    }
    if (errorAMG > maxErrorConstant * dx * dx) {
      MayDay::Abort("PetscRegular -- the error is larger than expected for second order");
    }

    pout() << "PetscRegular -- all tests passed" << endl;
  }
#else
  MayDay::Abort("PetscRegular -- this test requires USE_PETSC=TRUE");
#endif

  ChomboDischarge::finalize();
}
