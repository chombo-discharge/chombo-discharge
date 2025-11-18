#include <CD_Driver.H>
#include <CD_WireWire.H>
#include <CD_DischargeInceptionStepper.H>
#include <CD_DischargeInceptionTagger.H>
#include <CD_DataParser.H>

using namespace ChomboDischarge;
using namespace Physics::DischargeInception;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  ParmParse pp;

  Real gamma = 0.0;
  Real p     = 1.0;
  Real T     = 300.0;

  pp.get("temperature", T);
  pp.get("pressure", p);
  pp.get("gamma_ions", gamma);

  const Real N = p * Units::atm2pascal / (Units::kb * T);

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

  // clang-format off
  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new WireWire());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto timestepper = RefCountedPtr<DischargeInceptionStepper<>>(new DischargeInceptionStepper<>());
  auto celltagger  = RefCountedPtr<DischargeInceptionTagger>(new DischargeInceptionTagger(amr, timestepper->getElectricField(), alphaEff));
  auto engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, celltagger));
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
