#include "ParmParse.H"

#include <CD_ItoKMCJSON.H>

using namespace ChomboDischarge;
using namespace Physics::ItoKMC;

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  std::string basename;
  {
    ParmParse pp("ItoKMCJSON");
    pp.get("algorithm", basename);
    setPoutBaseName(basename);
  }

  Random::seed();

  auto      physics          = RefCountedPtr<ItoKMCPhysics>(new ItoKMCJSON());
  const int numPlasmaSpecies = physics->getNumPlasmaSpecies();
  const int numPhotonSpecies = physics->getNumPhotonSpecies();

  Vector<Real>     particles(numPlasmaSpecies, 0.0);
  Vector<Real>     photons(numPhotonSpecies, 0.0);
  Vector<Real>     phi(numPlasmaSpecies, 0.0);
  Vector<RealVect> gradPhi(numPlasmaSpecies, RealVect::Zero);

  Real initParticles;
  Real criticalDt;
  Real nonCriticalDt;
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

    physics->advanceKMC(particles,
                        photons,
                        criticalDt,
                        nonCriticalDt,
                        phi,
                        gradPhi,
                        dt,
                        E * RealVect::Unit,
                        RealVect::Zero,
                        1.0,
                        1.0);

    t += dt;
  }
  printToFile();

  physics->killKMC();

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
