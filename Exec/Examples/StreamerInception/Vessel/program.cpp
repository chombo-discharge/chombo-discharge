#include <CD_Driver.H>
#include <CD_TracerParticleSolver.H>
#include <CD_Vessel.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_StreamerInceptionStepper.H>
#include <CD_StreamerInceptionTagger.H>
#include <CD_LookupTable.H>
#include <CD_DataParser.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::StreamerInception;

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Read ionization and attachment coefficients and make them into functions.
  constexpr Real N = 2.45E25;

  // clang-format off
  LookupTable<2> ionizationData = DataParser::fractionalFileReadASCII("transport_data.txt", "E/N (Td)	Townsend ioniz. coef. alpha/N (m2)", "");
  LookupTable<2> attachmentData = DataParser::fractionalFileReadASCII("transport_data.txt", "E/N (Td)	Townsend attach. coef. eta/N (m2)", "");
  // clang-format on

  ionizationData.setRange(10, 2000, 0);
  attachmentData.setRange(10, 2000, 0);

  ionizationData.sort(0);
  attachmentData.sort(0);

  ionizationData.setTableSpacing(TableSpacing::Exponential);
  attachmentData.setTableSpacing(TableSpacing::Exponential);

  ionizationData.scale<0>(N * 1.E-21);
  attachmentData.scale<0>(N * 1.E-21);

  ionizationData.scale<1>(N);
  attachmentData.scale<1>(N);

  ionizationData.makeUniform(500);
  attachmentData.makeUniform(500);

  auto alpha = [&](const Real& E) -> Real { return ionizationData.getEntry<1>(E); };

  auto eta = [&](const Real& E) -> Real { return attachmentData.getEntry<1>(E); };

  // Define a background ionization rate.
  auto bgIonization = [N](const Real& E) -> Real { return 2.E6 / (1.17E-4 * exp(2.91E7 / E)); };

  // Define ion mobility and density
  auto ionMobility = [](const Real& E) -> Real { return 2E-4; };

  auto ionDensity = [](const RealVect& x) -> Real { return 1.E10; };

  //

  // Define a lightning impulse voltage curve.
  ParmParse vessel("impulse");
  Real      V0 = 1.0;
  Real      t0 = 0.0;
  Real      t1 = 1.2E-6;
  Real      t2 = 50E-6;

  vessel.get("voltage", V0);
  vessel.get("start", t0);
  vessel.get("t1", t1);
  vessel.get("t2", t2);

  auto voltageCurve = [V0, t0, t1, t2](const Real a_time) -> Real {
    constexpr Real alpha = 1.0 / 50E-6;
    constexpr Real beta  = 1.0 / 1.2E-6;

    return V0 * (exp(-(a_time + t0) / t1) - exp(-(a_time + t0) / t2));
  };

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new Vessel());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

  // Set up time stepper
  auto timestepper = RefCountedPtr<StreamerInceptionStepper<>>(new StreamerInceptionStepper<>());
  auto celltagger =
    RefCountedPtr<StreamerInceptionTagger>(new StreamerInceptionTagger(amr, timestepper->getElectricField()));

  // Set everything.
  timestepper->setAlpha(alpha);
  timestepper->setEta(eta);
  timestepper->setBackgroundRate(bgIonization);
  timestepper->setVoltageCurve(voltageCurve);
  timestepper->setNegativeIonMobility(ionMobility);
  timestepper->setNegativeIonDensity(ionDensity);

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, celltagger));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
