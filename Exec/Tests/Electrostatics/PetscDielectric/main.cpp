/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Test of FieldSolverAMG with an electrode and a dielectric, against FieldSolverGMG, on the RodSphere test
 * setup.
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>
#include <cmath>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_FieldSolverGMG.H>
#include <CD_FieldSolverAMG.H>
#include <CD_RodDielectric.H>
#include <CD_FieldStepper.H>
#include <CD_ParallelOps.H>
#include <CD_DataOps.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

/**
 * @brief Set up a field solver on an existing mesh, as FieldStepper does.
 * @details The electrode is at potential 1, the domain sides are as in the input file, and there is no space charge.
 * @tparam T Field solver class, e.g. FieldSolverAMG or FieldSolverGMG
 * @param[in] a_amr      AMR mesh
 * @param[in] a_compgeom Computational geometry
 * @param[in] a_realm    Realm
 * @return The solver.
 */
template <class T>
static RefCountedPtr<FieldSolver>
makeSolver(const RefCountedPtr<AmrMesh>&               a_amr,
           const RefCountedPtr<ComputationalGeometry>& a_compgeom,
           const std::string&                          a_realm)
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
  solver->setRho([](const RealVect a_pos) -> Real {
    return 0.0;
  });

  return solver;
}

/**
 * @brief Get the max-norm difference between two potentials.
 * @param[in] a_amr   AMR mesh
 * @param[in] a_realm Realm
 * @param[in] a_phi   Potential
 * @param[in] a_other Other potential
 * @return Max-norm difference over the cells of both phases that are neither covered by the EB nor by a finer level, on
 * all levels.
 */
static Real
maxDifference(const AmrMesh&       a_amr,
              const std::string&   a_realm,
              const MFAMRCellData& a_phi,
              const MFAMRCellData& a_other)
{
  Real maxDiff = 0.0;

  for (int lvl = 0; lvl <= a_amr.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl  = a_amr.getGrids(a_realm)[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(max : maxDiff)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&     din   = dit[mybox];
      const BaseFab<bool>& valid = (*a_amr.getValidCells(a_realm)[lvl])[din];

      for (int iphase = 0; iphase < (*a_phi[lvl])[din].numPhases(); iphase++) {
        const EBISBox&   ebisbox = (*a_phi[lvl])[din].getPhase(iphase).getEBISBox();
        const FArrayBox& phi     = (*a_phi[lvl])[din].getPhase(iphase).getFArrayBox();
        const FArrayBox& other   = (*a_other[lvl])[din].getPhase(iphase).getFArrayBox();

        for (BoxIterator bit(dbl[din]); bit.ok(); ++bit) {
          const IntVect iv = bit();

          if (valid(iv) && !(ebisbox.isCovered(iv))) {
            maxDiff = std::max(maxDiff, std::abs(phi(iv, 0) - other(iv, 0)));
          }
        }
      }
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
    auto        compgeom = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
    auto        amr      = RefCountedPtr<AmrMesh>(new AmrMesh());
    std::string realm;

    ParmParse ppStepper("FieldStepper");
    ppStepper.get("realm", realm);

    auto tagger      = RefCountedPtr<CellTagger>(nullptr);
    auto timestepper = RefCountedPtr<FieldStepper<FieldSolverGMG>>(new FieldStepper<FieldSolverGMG>());
    auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

    engine->setupAndRun();

    Real maxAllowedDifference;
    bool writePlotFiles;

    // The 2D inputs allow 1e-6, the 3D inputs 2e-2. In 3D the two solvers agree to the solver tolerance except in solid
    // cut cells with volume fraction 3.7e-4 on the sphere, where they differ by up to 8e-3 also at a relative residual
    // of 1e-15. This is a difference of the discretizations in those cells, not of the convergence.
    ParmParse pp("PetscDielectric");
    pp.get("max_difference", maxAllowedDifference);
    pp.get("write_plot_files", writePlotFiles);

    RefCountedPtr<FieldSolver> amg = makeSolver<FieldSolverAMG>(amr, compgeom, realm);
    RefCountedPtr<FieldSolver> gmg = makeSolver<FieldSolverGMG>(amr, compgeom, realm);

    // A constant surface charge density on the cut cells. Only the interface cells use it.
    Real sigmaValue;
    pp.get("sigma", sigmaValue);

    EBAMRIVData sigma;
    amr->allocate(sigma, realm, phase::gas, 1);
    DataOps::setValue(sigma, sigmaValue);

    const bool amgConverged = amg->solve(amg->getPotential(), amg->getRho(), sigma, true);
    const bool gmgConverged = gmg->solve(gmg->getPotential(), gmg->getRho(), sigma, true);

    const Real diff = maxDifference(*amr, realm, amg->getPotential(), gmg->getPotential());

    // FieldSolverAMG.step0000000.<dim>d.hdf5 and FieldSolverGMG.step0000000.<dim>d.hdf5, with the plot variables of
    // each solver.
    if (writePlotFiles) {
      amg->writePlotFile();
      gmg->writePlotFile();
    }

    pout() << "PetscDielectric -- levels = " << 1 + amr->getFinestLevel() << ", coarsest dx = " << amr->getDx()[0]
           << ", max |AMG - GMG| = " << diff << endl;

    if (!amgConverged) {
      MayDay::Abort("PetscDielectric -- FieldSolverAMG did not converge");
    }
    if (!gmgConverged) {
      MayDay::Abort("PetscDielectric -- FieldSolverGMG did not converge");
    }
    if (diff > maxAllowedDifference) {
      MayDay::Abort("PetscDielectric -- FieldSolverAMG and FieldSolverGMG differ by more than max_difference");
    }

    // Change the permittivities of the set-up solvers, as the semi-implicit schemes of CdrPlasma and ItoKMC do: reset
    // them, double the gas-phase permittivity on the faces and the EB, and hand them back with setSolverPermittivities.
    // This refills the AMG matrix and rebuilds its preconditioner.
    auto changePermittivities = [&](RefCountedPtr<FieldSolver>& a_solver) -> void {
      a_solver->setPermittivities();

      MFAMRCellData& permCell = a_solver->getPermittivityCell();
      MFAMRFluxData& permFace = a_solver->getPermittivityFace();
      MFAMRIVData&   permEB   = a_solver->getPermittivityEB();

      EBAMRFluxData permFaceGas = amr->alias(phase::gas, permFace);
      EBAMRIVData   permEBGas   = amr->alias(phase::gas, permEB);

      DataOps::incr(permFaceGas, permFaceGas, 1.0);
      DataOps::incr(permEBGas, permEBGas, 1.0);

      amr->arithmeticAverage(permFaceGas, realm, phase::gas);
      amr->arithmeticAverage(permEBGas, realm, phase::gas);

      a_solver->setSolverPermittivities(permCell, permFace, permEB);
    };

    changePermittivities(amg);
    changePermittivities(gmg);

    const bool amgConvergedChanged = amg->solve(amg->getPotential(), amg->getRho(), sigma, true);
    const bool gmgConvergedChanged = gmg->solve(gmg->getPotential(), gmg->getRho(), sigma, true);

    const Real diffChanged = maxDifference(*amr, realm, amg->getPotential(), gmg->getPotential());

    pout() << "PetscDielectric -- with the gas permittivity doubled through setSolverPermittivities, max |AMG - GMG| = "
           << diffChanged << endl;

    if (!amgConvergedChanged || !gmgConvergedChanged) {
      MayDay::Abort("PetscDielectric -- a solver did not converge after setSolverPermittivities");
    }
    if (diffChanged > maxAllowedDifference) {
      MayDay::Abort(
        "PetscDielectric -- after setSolverPermittivities, FieldSolverAMG and FieldSolverGMG differ by more "
        "than max_difference");
    }

    pout() << "PetscDielectric -- all tests passed" << endl;
  }
#else
  MayDay::Abort("PetscDielectric -- this test requires USE_PETSC=TRUE");
#endif

  ChomboDischarge::finalize();
}
