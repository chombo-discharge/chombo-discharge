#include <CD_Driver.H>
#include <CD_CdrCTU.H>
#include <CD_RodDielectric.H>
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionTagger.H>
#include <CD_DischargeIO.H>
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

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

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

  bool notCoarsest = false;

  // Storage for max, L1, and L2 solution error norms.
  std::vector<std::array<Real, 3>> norms;

  // Run the various cases
  for (const auto& cells : nCells) {

    // Set geometry and AMR
    RefCountedPtr<ComputationalGeometry> compgeom   = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
    RefCountedPtr<AmrMesh>               amr        = RefCountedPtr<AmrMesh>(new AmrMesh());
    RefCountedPtr<GeoCoarsener>          geocoarsen = RefCountedPtr<GeoCoarsener>(new GeoCoarsener());

    // Set up basic AdvectionDiffusion
    auto solver      = RefCountedPtr<CdrSolver>(new CdrCTU());
    auto timestepper = RefCountedPtr<AdvectionDiffusionStepper>(new AdvectionDiffusionStepper(solver));
    auto tagger      = RefCountedPtr<CellTagger>(new AdvectionDiffusionTagger(solver, amr));

    amr->setCoarsestGrid(cells);
    amr->buildDomains();

    // Set up the Driver and run it
    RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
    engine->setupAndRun(input_file);

    // Compute the solution errors if there is a coarser solution available.
    if (notCoarsest) {
      fineDBL   = amr->getGrids("primal")[0];
      fineEBISL = amr->getEBISLayout("primal", phase::gas)[0];

      // Coarsen the fine-grid solution onto the coarse-grid layout
      EBCoarAve aveOp(fineDBL,
                      coarDBL,
                      fineEBISL,
                      coarEBISL,
                      coarDomain,
                      refRat,
                      &(*(amr->getEBIndexSpace(phase::gas))));
      aveOp.averageData(*finePhi[0], *(solver->getPhi()[0]), Interval(0, 0), Average::Conservative);

      // Compute the error average(phi_fine) - coar_phi
      DataOps::incr(finePhi, coarPhi, -1.0);

      // Computes max, L1, and L2 norms.
      const Real Linf = DataOps::norm(*finePhi[0], amr->getDomains()[0], 0, true);
      const Real L1   = DataOps::norm(*finePhi[0], amr->getDomains()[0], 1, true);
      const Real L2   = DataOps::norm(*finePhi[0], amr->getDomains()[0], 2, true);

      norms.emplace_back(std::array<Real, 3>{Linf, L1, L2});
    }

    // Store coarse stuff for next time step.
    coarDBL    = amr->getGrids("primal")[0];
    coarDomain = amr->getDomains()[0];
    coarEBISL  = amr->getEBISLayout("primal", phase::gas)[0];

    // Allocate storage for coarse solution and coarsened fine solution.
    amr->allocate(coarPhi, "primal", phase::gas, 1);
    amr->allocate(finePhi, "primal", phase::gas, 1);

    // Copy the current solution, it becomes the "coarse" solution for the next iteration.
    coarPhi.copy(solver->getPhi());

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
  MPI_Finalize();
#endif
}
