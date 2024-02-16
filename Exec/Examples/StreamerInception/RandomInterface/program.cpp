#include <CD_Driver.H>
#include <CD_RandomInterface.H>
#include <CD_StreamerInceptionStepper.H>
#include <CD_StreamerInceptionTagger.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::StreamerInception;

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Read the input file into the ParmParse table
  const std::string input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());

  // Set geometry and AMR 
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new RandomInterface());
  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());

  // Define transport data
  auto alpha         = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };
  auto eta           = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };
  auto alphaEff      = [&](const Real& E, const RealVect& x) -> Real { return alpha(E,x) - eta(E,x); };
  auto bgRate        = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };
  auto detachRate    = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };
  auto fieldEmission = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };
  auto secondCoeff   = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };
  auto ionMobility   = [&](const Real& E) -> Real { return 0.0; };
  auto ionDiffusion  = [&](const Real& E) -> Real { return 0.0; };
  auto ionDensity    = [&](const RealVect& x) -> Real { return 0.0; };
  auto voltageCurve  = [&](const Real& time) -> Real { return 1.0;};

  // Set up time stepper 
  auto timestepper = RefCountedPtr<StreamerInceptionStepper<>> (new StreamerInceptionStepper<>());
  auto celltagger  = RefCountedPtr<StreamerInceptionTagger> (new StreamerInceptionTagger(amr, timestepper->getElectricField(), alphaEff));

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

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, celltagger));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
