#include <CD_Driver.H>
#include <CD_ItoKMCJSON.H>
#include <CD_ParallelOps.H>

// TLDR: This program runs many independent realizations of the two-group chemistry in chemistry.json,
//
//       e  -> el + el + M+   (slow branching)
//       el -> e              (fast relaxation)
//
//       from a single electron and compares the mean total electron number e + el at stop_time with
//       the exact mean from the moment equations, which close because every reaction is first order.

using namespace ChomboDischarge;
using namespace Physics::ItoKMC;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  Random::seed();

  ParmParse pp;

  int  numRuns;
  Real stopTime;
  Real maxDt;
  Real exactMean;

  pp.get("num_runs", numRuns);
  pp.get("stop_time", stopTime);
  pp.get("max_dt", maxDt);
  pp.get("exact_mean", exactMean);

  auto      physics          = RefCountedPtr<ItoKMCPhysics>(new ItoKMCJSON());
  const int numPlasmaSpecies = physics->getNumPlasmaSpecies();
  const int numPhotonSpecies = physics->getNumPhotonSpecies();

  Vector<Real>     particles(numPlasmaSpecies, 0.0);
  Vector<Real>     photons(numPhotonSpecies, 0.0);
  Vector<Real>     phi(numPlasmaSpecies, 0.0);
  Vector<RealVect> gradPhi(numPlasmaSpecies, RealVect::Zero);

  physics->defineKMC();

  // Every rate in chemistry.json is a constant, so the field is irrelevant and is set to zero.
  Real sum = 0.0;

  for (int run = 0; run < numRuns; run++) {
    for (int i = 0; i < numPlasmaSpecies; i++) {
      particles[i] = 0.0;
    }

    // One electron in the reacting group ("e" is the first plasma species).
    particles[0] = 1.0;

    for (Real t = 0.0; t < stopTime;) {
      const Real dt = std::min(maxDt, stopTime - t);

      Real physicsDt = std::numeric_limits<Real>::max();

      physics->advanceKMC(particles, photons, physicsDt, phi, gradPhi, dt, RealVect::Zero, RealVect::Zero, 1.0, 1.0);

      t += dt;
    }

    sum += particles[0] + particles[1];
  }

  physics->killKMC();

  // Mean over all realizations on all ranks.
  const Real totalSum  = ParallelOps::sum(sum);
  const Real totalRuns = ParallelOps::sum(1.0 * numRuns);
  const Real mean      = totalSum / totalRuns;

  pout() << "realizations               = " << llround(totalRuns) << endl;
  pout() << "mean total electron number = " << mean << endl;
  pout() << "exact mean                 = " << exactMean << endl;
  pout() << "relative deviation         = " << 100.0 * (mean - exactMean) / exactMean << " %" << endl;

  ChomboDischarge::finalize();
}
