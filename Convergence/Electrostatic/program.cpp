#include <CD_Driver.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_CoaxialCable.H>
#include <CD_FieldStepper.H>
#include <CD_DischargeIO.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  std::vector<IntVect> nCells{32 * IntVect::Unit,
                              64 * IntVect::Unit,
                              128 * IntVect::Unit,
                              256 * IntVect::Unit,
                              512 * IntVect::Unit};

  std::vector<std::array<Real, 3>> norms;

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Get R0, R1 etc from the input script
  ParmParse ppOuter("CoaxialCable.outer");
  ParmParse ppInner("CoaxialCable.inner");
  ParmParse ppDiele("CoaxialCable.dielectric");

  // Parameters from input script
  Real           R0, R1, R2, eps, a;
  constexpr Real phi0 = 1.0;

  ppInner.get("radius", R0);
  ppOuter.get("radius", R2);
  ppDiele.get("radius", R1);
  ppDiele.get("eps", eps);
  a = phi0 / (log(R1 / R2) - log(R1 / R0) / eps);

  // Set up the exact solution.
  auto exactSolution = [R0, R1, R2, eps, phi0, a](const RealVect& x) -> Real {
    const Real r = sqrt(x[0] * x[0] + x[1] * x[1]);

    Real phi = 0.0;

    if (r >= R0 && r <= R1) {
      phi = phi0 + a / eps * log(r / R0);
    }
    else if (r >= R1 && r <= R2) {
      phi = a * log(r / R2);
    }

    return phi;
  };

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom   = RefCountedPtr<ComputationalGeometry>(new CoaxialCable());
  RefCountedPtr<AmrMesh>               amr        = RefCountedPtr<AmrMesh>(new AmrMesh());
  RefCountedPtr<GeoCoarsener>          geocoarsen = RefCountedPtr<GeoCoarsener>(new GeoCoarsener());
  RefCountedPtr<CellTagger>            tagger     = RefCountedPtr<CellTagger>(NULL);

  // Set up basic Poisson, potential = 1
  auto timestepper = RefCountedPtr<FieldStepper<FieldSolverMultigrid>>(new FieldStepper<FieldSolverMultigrid>());

  // Run the various cases
  for (const auto& cells : nCells) {
    amr->setCoarsestGrid(cells);
    amr->buildDomains();

    // Set up the Driver and run our program.
    RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
    engine->setupAndRun(input_file);

    // Compute the error
    MFAMRCellData err;
    amr->allocate(err, "primal", 1);
    DataOps::setValue(err, exactSolution, amr->getProbLo(), amr->getDx(), 0);
    DataOps::incr(err, timestepper->getPotential(), -1.0);

    // Extract error from gas and solid sides
    EBAMRCellData errGas = amr->alias(phase::gas, err);
    //    EBAMRCellData errSol = amr->alias(phase::solid, err);

    // Compute the various norms
    const Real Linf = DataOps::norm(*errGas[0], amr->getDomains()[0], 0, true);
    const Real L1   = DataOps::norm(*errGas[0], amr->getDomains()[0], 1, true);
    const Real L2   = DataOps::norm(*errGas[0], amr->getDomains()[0], 2, true);

    norms.emplace_back(std::array<Real, 3>{Linf, L1, L2});
  }

  // Compute convergence rates
  if (procID() == 0) {
    std::cout << "# cells"
              << "\t"
              << "Linf error"
              << "\t"
              << "L1 error"
              << "\t"
              << "L2 error" << std::endl;
    for (int i = 0; i < norms.size(); i++) {
      std::cout << nCells[i][0] << "\t" << std::get<0>(norms[i]) << "\t" << std::get<1>(norms[i]) << "\t"
                << std::get<2>(norms[i]) << std::endl;
    }
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}
