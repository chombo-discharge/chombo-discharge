#include <fstream>

#include <nlohmann/json.hpp>

#include <CD_Driver.H>
#include <CD_ItoKMCJSON.H>
#include <CD_ParallelOps.H>

// TLDR: This program runs many independent realizations of the two-group chemistry in chemistry.json,
//
//       e  -> el + el + M+   (slow branching, rate nu_i)
//       el -> e              (fast relaxation, rate nu_r)
//
//       from initial_particles electrons in the e group, and prints the mean populations after every
//       step together with the solution of the moment equations
//
//       d/dt <e>  = -nu_i <e> + nu_r <el>
//       d/dt <el> = 2 nu_i <e> - nu_r <el>
//
//       which are exact because every reaction is first order.

using namespace ChomboDischarge;
using namespace Physics::ItoKMC;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  Random::seed();

  ParmParse pp;

  int  numRuns;
  int  initialParticles;
  Real stopTime;
  Real maxDt;

  pp.get("num_runs", numRuns);
  pp.get("initial_particles", initialParticles);
  pp.get("stop_time", stopTime);
  pp.get("max_dt", maxDt);

  // Rates of the two reactions, read from the chemistry file so that the moment equations always match the
  // chemistry the solver ran.
  std::string chemistryFile;

  ParmParse ppJSON("ItoKMCJSON");

  ppJSON.get("chemistry_file", chemistryFile);

  std::ifstream chemistryStream(chemistryFile);

  nlohmann::json chemistry;

  chemistryStream >> chemistry;

  const Real nuI = chemistry["plasma reactions"][0]["value"].get<Real>();
  const Real nuR = chemistry["plasma reactions"][1]["value"].get<Real>();

  // The step sequence is the same for every realization, so the populations after each step can be
  // accumulated by step index.
  std::vector<Real> times;

  for (Real t = 0.0; t < stopTime;) {
    const Real dt = std::min(maxDt, stopTime - t);

    t += dt;

    times.push_back(t);
  }

  const size_t numSteps = times.size();

  std::vector<Real> sumE(numSteps, 0.0);
  std::vector<Real> sumEl(numSteps, 0.0);
  std::vector<Real> sumTotal2(numSteps, 0.0);

  auto      physics          = RefCountedPtr<ItoKMCPhysics>(new ItoKMCJSON());
  const int numPlasmaSpecies = physics->getNumPlasmaSpecies();
  const int numPhotonSpecies = physics->getNumPhotonSpecies();

  Vector<Real>     particles(numPlasmaSpecies, 0.0);
  Vector<Real>     photons(numPhotonSpecies, 0.0);
  Vector<Real>     phi(numPlasmaSpecies, 0.0);
  Vector<RealVect> gradPhi(numPlasmaSpecies, RealVect::Zero);

  physics->defineKMC();

  // Every rate in chemistry.json is a constant, so the field is irrelevant and is set to zero.
  for (int run = 0; run < numRuns; run++) {
    for (int i = 0; i < numPlasmaSpecies; i++) {
      particles[i] = 0.0;
    }

    // Start in the reacting group ("e" is the first plasma species).
    particles[0] = 1.0 * initialParticles;

    Real t = 0.0;

    for (size_t k = 0; k < numSteps; k++) {
      const Real dt = times[k] - t;

      Real physicsDt = std::numeric_limits<Real>::max();

      physics->advanceKMC(particles, photons, physicsDt, phi, gradPhi, dt, RealVect::Zero, RealVect::Zero, 1.0, 1.0);

      t = times[k];

      const Real total = particles[0] + particles[1];

      sumE[k] += particles[0];
      sumEl[k] += particles[1];
      sumTotal2[k] += total * total;
    }
  }

  physics->killKMC();

  // Integrate the moment equations with RK4 on a fine substep; the fast eigenvalue is nu_i + nu_r.
  std::vector<Real> exactE(numSteps, 0.0);
  std::vector<Real> exactEl(numSteps, 0.0);

  {
    Real e  = 1.0 * initialParticles;
    Real el = 0.0;
    Real t  = 0.0;

    auto rhs = [&](const Real a_e, const Real a_el, Real& a_de, Real& a_del) -> void {
      a_de  = -nuI * a_e + nuR * a_el;
      a_del = 2.0 * nuI * a_e - nuR * a_el;
    };

    for (size_t k = 0; k < numSteps; k++) {
      const int  numSub = 1000;
      const Real h      = (times[k] - t) / numSub;

      for (int s = 0; s < numSub; s++) {
        Real k1e, k1l, k2e, k2l, k3e, k3l, k4e, k4l;

        rhs(e, el, k1e, k1l);
        rhs(e + 0.5 * h * k1e, el + 0.5 * h * k1l, k2e, k2l);
        rhs(e + 0.5 * h * k2e, el + 0.5 * h * k2l, k3e, k3l);
        rhs(e + h * k3e, el + h * k3l, k4e, k4l);

        e += h * (k1e + 2.0 * k2e + 2.0 * k3e + k4e) / 6.0;
        el += h * (k1l + 2.0 * k2l + 2.0 * k3l + k4l) / 6.0;
      }

      t = times[k];

      exactE[k]  = e;
      exactEl[k] = el;
    }
  }

  // Reduce over ranks and print the time series; every rank prints the same reduced data.
  const Real totalRuns = ParallelOps::sum(1.0 * numRuns);

  pout() << "# realizations = " << llround(totalRuns) << ", initial particles = " << initialParticles
         << ", nu_i = " << nuI << ", nu_r = " << nuR << endl;
  pout() << "# time  <e>  <el>  <e+el>  stderr(<e+el>)  e_moments  el_moments  total_moments  deviation(%)" << endl;

  Real finalDeviation = 0.0;

  for (size_t k = 0; k < numSteps; k++) {
    const Real meanE   = ParallelOps::sum(sumE[k]) / totalRuns;
    const Real meanEl  = ParallelOps::sum(sumEl[k]) / totalRuns;
    const Real meanTot = meanE + meanEl;
    const Real meanTo2 = ParallelOps::sum(sumTotal2[k]) / totalRuns;
    const Real stdErr  = sqrt(std::max(meanTo2 - meanTot * meanTot, 0.0) / totalRuns);
    const Real exactT  = exactE[k] + exactEl[k];

    finalDeviation = 100.0 * (meanTot - exactT) / exactT;

    pout() << times[k] << "\t" << meanE << "\t" << meanEl << "\t" << meanTot << "\t" << stdErr << "\t" << exactE[k]
           << "\t" << exactEl[k] << "\t" << exactT << "\t" << finalDeviation << endl;
  }

  pout() << "# final relative deviation = " << finalDeviation << " %" << endl;

  ChomboDischarge::finalize();
}
