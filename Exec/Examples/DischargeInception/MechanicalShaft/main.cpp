#include <CD_Driver.H>
#include <CD_MechanicalShaft.H>
#include <CD_DischargeInceptionStepper.H>
#include <CD_DischargeInceptionTagger.H>
#include <CD_LookupTable1D.H>
#include <CD_DataParser.H>

using namespace ChomboDischarge;
using namespace Physics::DischargeInception;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  const Real N = 2.45E25;
  const Real T = 300.0;

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
    return 1.0;
  };
  auto fieldEmission = [&](const Real& E, const RealVect& x) -> Real {
    const Real beta = 1.0; // Field enhancement factor
    const Real phi  = 4.5;
    const Real C1   = 1.54E-6 * std::pow(10, 4.52 / sqrt(phi)) / phi;
    const Real C2   = 2.84E9 * std::pow(phi, 1.5);

    return C1 * (E * E) * exp(-C2 / (beta * E));
  };
  auto secondCoeff = [&](const Real& E, const RealVect& x) -> Real {
    return 1E-4;
  };

  // clang-format off
  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new MechanicalShaft());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto timestepper = RefCountedPtr<DischargeInceptionStepper<>>(new DischargeInceptionStepper<>());
  auto celltagger  = RefCountedPtr<DischargeInceptionTagger>(new DischargeInceptionTagger(amr, timestepper->getElectricField(), alphaEff));
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, celltagger));
  // clang-format on

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

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
