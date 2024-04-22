#include <CD_Driver.H>
#include <CD_Vessel.H>
#include <CD_DischargeInceptionStepper.H>
#include <CD_DischargeInceptionTagger.H>
#include <CD_LookupTable.H>
#include <CD_DataParser.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::DischargeInception;

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Read the input file into the ParmParse table
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Read BOLSIG+ data into alpha and eta coefficients
  const Real N  = 2.45E25;
  const Real O2 = 0.2;
  const Real N2 = 0.8;
  const Real T  = 300.0;

  LookupTable1D<> ionizationData = DataParser::fractionalFileReadASCII("transport_data.txt",
                                                                       "E/N (Td)	Townsend ioniz. coef. alpha/N (m2)",
                                                                       "");
  LookupTable1D<> attachmentData = DataParser::fractionalFileReadASCII("transport_data.txt",
                                                                       "E/N (Td)	Townsend attach. coef. eta/N (m2)",
                                                                       "");

  ionizationData.truncate(10, 2000, 0);
  attachmentData.truncate(10, 2000, 0);

  ionizationData.scale<0>(N * 1.E-21);
  attachmentData.scale<0>(N * 1.E-21);

  ionizationData.scale<1>(N);
  attachmentData.scale<1>(N);

  ionizationData.prepareTable(0, 500, LookupTable::Spacing::Exponential);
  attachmentData.prepareTable(0, 500, LookupTable::Spacing::Exponential);

  // Read data for the voltage curve.
  Real peak           = 0.0;
  Real t0             = 0.0;
  Real t1             = 0.0;
  Real t2             = 0.0;
  Real secondTownsend = 0.0;

  ParmParse li("lightning_impulse");
  ParmParse di("DischargeInception");
  li.get("peak", peak);
  li.get("start", t0);
  li.get("tail_time", t1);
  li.get("front_time", t2);

  di.get("second_townsend", secondTownsend);

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new Vessel());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

  // Define transport data
  auto alpha = [&](const Real& E, const RealVect& x) -> Real {
    return ionizationData.interpolate<1>(E);
  };
  auto eta = [&](const Real& E, const RealVect& x) -> Real {
    return attachmentData.interpolate<1>(E);
  };
  auto alphaEff = [&](const Real& E, const RealVect x) -> Real {
    return alpha(E, x) - eta(E, x);
  };
  auto bgRate = [&](const Real& E, const RealVect& x) -> Real {
    return 0.0;
  };
  auto detachRate = [&](const Real& E, const RealVect& x) -> Real {
    const Real Etd = E / (N * 1E-21);
    return 1.24E-11 * 1E-6 * N * exp(-std::pow((179.0 / (8.8 + Etd)), 2));
  };
  auto ionMobility = [&](const Real& E) -> Real {
    return 2E-4;
  };
  auto ionDiffusion = [&](const Real& E) -> Real {
    return ionMobility(E) * Units::kb * T / Units::Qe;
  };
  auto ionDensity = [&](const RealVect& x) -> Real {
    return 4.E6;
  };
  auto voltageCurve = [&](const Real& t) -> Real {
    return peak * (exp(-(t + t0) / t1) - exp(-(t + t0) / t2));
  };
  auto fieldEmission = [&](const Real& E, const RealVect& x) -> Real {
    const Real beta = 1.0; // Field enhancement factor
    const Real phi  = 4.5;
    const Real C1   = 1.54E-6 * std::pow(10, 4.52 / sqrt(phi)) / phi;
    const Real C2   = 2.84E9 * std::pow(phi, 1.5);

    return C1 * (E * E) * exp(-C2 / (beta * E));
  };
  auto secondaryEmission = [&](const Real& E, const RealVect& x) -> Real {
    return secondTownsend;
  };

  // Set up time stepper
  auto timestepper = RefCountedPtr<DischargeInceptionStepper<>>(new DischargeInceptionStepper<>());
  auto celltagger  = RefCountedPtr<DischargeInceptionTagger>(
    new DischargeInceptionTagger(amr, timestepper->getElectricField(), alphaEff));

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
  timestepper->setSecondaryEmission(secondaryEmission);

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, celltagger));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
