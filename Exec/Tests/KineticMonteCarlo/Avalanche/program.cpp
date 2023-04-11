#include <CD_Driver.H>
#include <CD_KMCSolver.H>
#include <CD_KMCSingleStateReaction.H>

using namespace ChomboDischarge;

// FPR is the state representation. Integer and floating point
// types should both work.
using FPR           = long long;
using KMCState      = KMCSingleState<FPR>;
using KMCReaction   = KMCSingleStateReaction<KMCState, FPR>;
using KMCSolverType = KMCSolver<KMCReaction, KMCState, FPR>;

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Read input file.
  ParmParse pp(argc - 2, argv + 2, NULL, argv[1]);

  // Seed the RNG
  Random::setRandomSeed();

  // State that we advance.
  KMCState state(1);

  // Define list of reactions.
  std::vector<std::shared_ptr<const KMCReaction>> reactionList;

  // Define reaction e + null -> e + e + null
  auto ionization = std::make_shared<KMCReaction>(std::list<size_t>{0}, std::list<size_t>{0, 0});

  // Define e -> null
  auto attachment = std::make_shared<KMCReaction>(std::list<size_t>{0}, std::list<size_t>{});

  reactionList.emplace_back(attachment);
  reactionList.emplace_back(ionization);

  // Read initial values.
  Real& ionizationRate = ionization->rate();
  Real& attachmentRate = attachment->rate();

  Real stopTime = 0.0;
  Real SSAlim   = 0.1;
  Real eps      = 0.3;
  int  initVal  = 0;
  int  numSteps = 100;
  int  numCrit  = 5;
  int  numSSA   = 5;

  std::string alg;

  // Read program input variables
  pp.get("nsteps", numSteps);
  pp.get("stop_time", stopTime);
  pp.get("ionization_rate", ionizationRate);
  pp.get("attachment_rate", attachmentRate);
  pp.get("initial_particles", initVal);
  pp.get("algorithm", alg);
  pp.get("num_crit", numCrit);
  pp.get("num_ssa", numSSA);
  pp.get("eps", eps);
  pp.get("ssa_lim", SSAlim);

  state[0] = (long long)initVal;

  // Define the Kinetic Monte Carlo solver and run it until time = 10.
  KMCSolverType kmcSolver(reactionList);

  kmcSolver.setSolverParameters(numCrit, numSSA, eps, SSAlim);

  // Advance problem.
  Real       curDt   = 0.0;
  const Real effRate = ionizationRate - attachmentRate;
  while (curDt < stopTime) {
    pout() << curDt << "\t" << state[0] << "\t" << 1.0 * initVal * exp(effRate * curDt) << endl;

    Real nextDt = std::numeric_limits<Real>::max();
    if (alg == "ssa") {
      nextDt = kmcSolver.getCriticalTimeStep(state);

      kmcSolver.stepSSA(state);
    }
    else if (alg == "tau") {
      nextDt = stopTime / numSteps;

      kmcSolver.advanceTauPlain(state, nextDt);
    }
    else if (alg == "hybrid") {
      nextDt = stopTime / numSteps;

      kmcSolver.advanceHybrid(state, nextDt);
    }
    else {
      const std::string err = "Expected algorithm to be 'ssa', 'tau', or 'hybrid' but got '" + alg + "'";

      MayDay::Error(err.c_str());
    }

    curDt += nextDt;
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
