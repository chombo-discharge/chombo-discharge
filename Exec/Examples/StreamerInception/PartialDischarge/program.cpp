#include <CD_Driver.H>
#include <CD_Aerosol.H>
#include <CD_StreamerInceptionStepper.H>
#include <CD_StreamerInceptionTagger.H>
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
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new Aerosol());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

  Real gamma  = 0.0;
  Real p      = 1.0;
  Real T      = 300.0;
  Real sigma0 = 0.0;

  pp.get("temperature", T);
  pp.get("pressure", p);
  pp.get("gamma_ions", gamma);
  pp.get("sigma0", sigma0);

  const Real N  = p * Units::atm2pascal / (Units::kb * T);
  const Real O2 = 0.2;
  const Real N2 = 0.8;

  // Read BOLSIG+ data into alpha and eta coefficients
  LookupTable1D<> ionizationData = DataParser::fractionalFileReadASCII("bolsig_air.dat",
                                                                       "E/N (Td)	Townsend ioniz. coef. alpha/N (m2)",
                                                                       "");
  LookupTable1D<> attachmentData = DataParser::fractionalFileReadASCII("bolsig_air.dat",
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

  // Define transport data
  auto alpha = [&](const Real& E, const RealVect& x) -> Real {
    return ionizationData.interpolate<1>(E);
  };
  auto eta = [&](const Real& E, const RealVect& x) -> Real {
    return attachmentData.interpolate<1>(E);
  };
  auto alphaEff = [&](const Real& E, const RealVect& x) -> Real {
    return alpha(E, x) - eta(E, x);
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
  auto secondCoeff = [&](const Real& E, const RealVect& x) -> Real {
    return gamma;
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
  auto sigma = [&](const RealVect& x) -> Real {
    Real xy = 0.0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      if (dir != 1) {
        xy += x[dir] * x[dir];
      }
    }

    xy = sqrt(xy);

    const Real theta = atan2(x[1], xy);

    return sigma0 * std::pow(sin(theta), 3);
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
  timestepper->setSecondaryEmission(secondCoeff);
  timestepper->setIonMobility(ionMobility);
  timestepper->setIonDiffusion(ionDiffusion);
  timestepper->setIonDensity(ionDensity);
  timestepper->setVoltageCurve(voltageCurve);
  timestepper->setSigma(sigma);

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, celltagger));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
