#include <CD_Driver.H>
#include <CD_NoisePlane.H>
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

  // Read the input file into the ParmParse table
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new NoisePlane());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

  // Read ionization and attachment coefficients and make them into functions.
  constexpr Real N = 2.45E25;

  LookupTable1D<2> ionizationData = DataParser::fractionalFileReadASCII("sf6.dat",
                                                                        "E/N (Td)	Townsend ioniz. coef. alpha/N (m2)",
                                                                        "");
  LookupTable1D<2> attachmentData = DataParser::fractionalFileReadASCII("sf6.dat",
                                                                        "E/N (Td)	Townsend attach. coef. eta/N (m2)",
                                                                        "");

  ionizationData.setRange(10, 4000, 0);
  attachmentData.setRange(10, 4000, 0);

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

  // Define transport data
  auto alpha = [&](const Real& E) -> Real {
    return ionizationData.getEntry<1>(E);
  };
  auto eta = [&](const Real& E) -> Real {
    return attachmentData.getEntry<1>(E);
  };
  auto alphaEff = [&](const Real& E) -> Real {
    return alpha(E) - eta(E);
  };
  auto bgRate = [&](const Real& E, const RealVect& x) -> Real {
    return 0.0;
  };
  auto detachRate = [&](const Real& E, const RealVect& x) -> Real {
    return 0.0;
  };
  auto fieldEmission = [&](const Real& E, const RealVect& x) -> Real {
    return 0.0;
  };
  auto ionMobility = [&](const Real& E) -> Real {
    return 0.0;
  };
  auto ionDiffusion = [&](const Real& E) -> Real {
    return 0.0;
  };
  auto ionDensity = [&](const RealVect& x) -> Real {
    return 0.0;
  };
  auto voltageCurve = [&](const Real& time) -> Real {
    return 1.0;
  };

  // Set up time stepper
  auto timestepper = RefCountedPtr<StreamerInceptionStepper<>>(new StreamerInceptionStepper<>());
  auto celltagger  = RefCountedPtr<StreamerInceptionTagger>(
    new StreamerInceptionTagger(amr, timestepper->getElectricField(), alphaEff));

  // Set transport data
  timestepper->setAlpha(alpha);
  timestepper->setEta(eta);
  timestepper->setBackgroundRate(bgRate);
  timestepper->setDetachmentRate(detachRate);
  timestepper->setFieldEmission(fieldEmission);
  timestepper->setIonMobility(ionMobility);
  timestepper->setIonDiffusion(ionDiffusion);
  timestepper->setIonDensity(ionDensity);
  timestepper->setVoltageCurve(voltageCurve);

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, celltagger));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
