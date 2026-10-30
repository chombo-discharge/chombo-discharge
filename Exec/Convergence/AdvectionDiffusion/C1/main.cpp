/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Spatial convergence test for advection-diffusion on a uniform grid.
 * @author Robert Marskar
 */

#include <CD_Driver.H>
#include <CD_CdrCTU.H>
#include <CD_RodDielectric.H>
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionTagger.H>
#include <CD_DischargeIO.H>
#include <CD_ParallelOps.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::AdvectionDiffusion;

// This program runs convergence testing on a uniform grid for an advection-diffusion problem. When
// there are two solutions available (coarse and fine solutions), we can coarsen the fine-grid solution
// and compute the coarse-grid error as
//
//    e = coarsen(phi_fine) - phi_coar
//
// We then compute the various error norms.

int
main(int argc, char* argv[])
{

  ChomboDischarge::initialize(argc, argv);

  // These are the grid resolutions that we run this program for. For the 32^2 grid the "exact solution" is the
  // coarsened solution of the 64^2 grid etc.
  constexpr int refRat = 2;
  constexpr int nComp  = 1;

  std::vector<IntVect> nCells{16 * IntVect::Unit,
                              32 * IntVect::Unit,
                              64 * IntVect::Unit,
                              128 * IntVect::Unit,
                              256 * IntVect::Unit,
                              512 * IntVect::Unit};

  // This stuff is required because the old/new solutions and grids are discarded
  // every time we reinitialize AmrMesh. So we store it here.
  //
  // coarPhi => Solution on coarser grid
  // finePhi => Coarsened fine-grid solution (i.e., coarsen(phi_fine) as explained above)
  // coarDBL => Coarse grid
  // fineDBL => Fine grid
  // coarEBISL => Coarse EB layout
  // fineEBISL => Fine EB layout
  // coarDomain => Coarse grid domain
  EBAMRCellData coarPhi;
  EBAMRCellData finePhi;

  DisjointBoxLayout coarDBL;
  DisjointBoxLayout fineDBL;

  EBISLayout coarEBISL;
  EBISLayout fineEBISL;

  ProblemDomain coarDomain;
  Real          coarDx = 0.0;

  bool notCoarsest = false;

  // Storage for max, L1, and L2 solution error norms.
  std::vector<std::array<Real, 3>> norms;

  // Run the various cases
  for (const auto& cells : nCells) {

    // Set geometry and AMR
    RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
    RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

    // Set up basic AdvectionDiffusion
    auto solver      = RefCountedPtr<CdrSolver>(new CdrCTU());
    auto timestepper = RefCountedPtr<AdvectionDiffusionStepper>(new AdvectionDiffusionStepper(solver));
    auto tagger      = RefCountedPtr<CellTagger>(new AdvectionDiffusionTagger(solver, amr));

    amr->setCoarsestGrid(cells);
    amr->buildDomains();

    // Set up the Driver and run it
    RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));
    engine->setupAndRun();

    // Compute the solution errors if there is a coarser solution available.
    if (notCoarsest) {
      fineDBL   = amr->getGrids("primal")[0];
      fineEBISL = amr->getEBISLayout("primal", phase::gas)[0];

      // Coarsen the fine-grid solution onto the coarse-grid layout
      EBCoarAve
        aveOp(fineDBL, coarDBL, fineEBISL, coarEBISL, coarDomain, refRat, &(*(amr->getEBIndexSpace(phase::gas))));
      aveOp.averageData(*finePhi[0], *(solver->getPhi()[0]), Interval(0, 0), Average::Conservative);

      // Compute the error average(phi_fine) - coar_phi
      DataOps::incr(finePhi, coarPhi, -1.0);

      // Max, L1, and L2 norms of the error on the coarse grid, over all VoFs (cut cells included). L1 and L2 are
      // weighted by the volume kappa*dx^D. The error lives on the coarse layout, so it is not iterated with the
      // AmrMesh, which now holds the fine grids.
      Real Linf = 0.0;
      Real L1   = 0.0;
      Real L2   = 0.0;

      const Real dV = std::pow(coarDx, SpaceDim);

      for (DataIterator dit = coarDBL.dataIterator(); dit.ok(); ++dit) {
        const EBCellFAB& err     = (*finePhi[0])[dit()];
        const EBISBox&   ebisbox = coarEBISL[dit()];

        for (VoFIterator vofit(IntVectSet(coarDBL[dit()]), ebisbox.getEBGraph()); vofit.ok(); ++vofit) {
          const Real e     = std::abs(err(vofit(), 0));
          const Real kappa = ebisbox.volFrac(vofit());

          Linf = std::max(Linf, e);
          L1 += kappa * e * dV;
          L2 += kappa * e * e * dV;
        }
      }

      Linf = ParallelOps::max(Linf);
      L1   = ParallelOps::sum(L1);
      L2   = std::sqrt(ParallelOps::sum(L2));

      norms.emplace_back(std::array<Real, 3>{Linf, L1, L2});
    }

    // Store coarse stuff for next time step.
    coarDBL    = amr->getGrids("primal")[0];
    coarDomain = amr->getDomains()[0];
    coarEBISL  = amr->getEBISLayout("primal", phase::gas)[0];
    coarDx     = amr->getDx()[0];

    // Allocate storage for coarse solution and coarsened fine solution.
    amr->allocate(coarPhi, "primal", phase::gas, 1);
    amr->allocate(finePhi, "primal", phase::gas, 1);

    // Copy the current solution, it becomes the "coarse" solution for the next iteration.
    DataOps::copy(coarPhi, solver->getPhi());

    notCoarsest = true;
  }

  // Print the solution errors
#ifdef CH_MPI
  if (procID() == 0) {
#endif
    // clang-format off
    std::cout << "# cells\t" << "Linf error\t" << "L1 error\t" << "L2 error\n";

    for (int i = 0; i < norms.size(); i++) {
      std::cout << nCells[i][0] << "\t"
		<< std::get<0>(norms[i]) << "\t"
		<< std::get<1>(norms[i]) << "\t"
                << std::get<2>(norms[i]) << "\n";
    }
    // clang-format on
#ifdef CH_MPI
  }
#endif

  ChomboDischarge::finalize();
}
