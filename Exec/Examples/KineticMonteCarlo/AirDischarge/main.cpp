#include <CD_Driver.H>
#include <CD_ItoKMCJSON.H>

using namespace ChomboDischarge;
using namespace Physics::ItoKMC;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  Random::seed();

  ParmParse pp;

  std::string basename;

  pp.get("ItoKMCJSON.algorithm", basename);
  setPoutBaseName(basename);

  auto      physics          = RefCountedPtr<ItoKMCPhysics>(new ItoKMCJSON());
  const int numPlasmaSpecies = physics->getNumPlasmaSpecies();
  const int numPhotonSpecies = physics->getNumPhotonSpecies();

  Vector<Real>     particles(numPlasmaSpecies, 0.0);
  Vector<Real>     photons(numPhotonSpecies, 0.0);
  Vector<Real>     phi(numPlasmaSpecies, 0.0);
  Vector<RealVect> gradPhi(numPlasmaSpecies, RealVect::Zero);

  Real initParticles;
  Real physicsDt;
  Real dt;
  Real E;
  Real stopTime;

  pp.get("num_initial_electrons", initParticles);
  pp.get("E", E);
  pp.get("dt", dt);
  pp.get("stop_time", stopTime);

  particles[0] = initParticles;

  Real t = 0.0;

  auto printToFile = [&]() -> void {
    pout() << t << "\t";
    for (int i = 0; i < numPlasmaSpecies; i++) {
      pout() << llround(particles[i]) << "\t";
    }
    pout() << endl;
  };

  physics->defineKMC();
  while (t < stopTime) {
    printToFile();

    physics->advanceKMC(particles, photons, physicsDt, phi, gradPhi, dt, E * RealVect::Unit, RealVect::Zero, 1.0, 1.0);

    t += dt;
  }
  printToFile();

  physics->killKMC();

  ChomboDischarge::finalize();
}
