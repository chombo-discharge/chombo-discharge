#include <CD_Driver.H>
#include <CD_KMCSolver.H>
#include <CD_KMCSingleStateReaction.H>

// TLDR: This program solves for the Schlogl model reactions
//
//       B1 + X + X -> X + X + X (rate c1)
//       X + X + X -> X + X + B1 (rate c2)
//       B2 -> X (rate c3)
//       X -> B2 (rate c4)
//
//       The B1/B2 species are buffered species whose populations don't change much, so we simply modify
//       the c1/c3 reaction rates to include them.

using namespace ChomboDischarge;

constexpr Real c1 = 3E-7;
constexpr Real c2 = 1E-4;
constexpr Real c3 = 1E-3;
constexpr Real c4 = 3.5;

constexpr Real N1 = 1.E5;
constexpr Real N2 = 2.E5;

// FPR is the state representation. Integer and floating point
// types should both work.
using FPR           = Real;
using KMCState      = KMCSingleState<FPR>;
using KMCReaction   = KMCSingleStateReaction<KMCState, FPR>;
using KMCSolverType = KMCSolver<KMCReaction, KMCState, FPR>;

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Read input file
  ParmParse pp(argc - 2, argv + 2, NULL, argv[1]);

  // Seed the RNG
  Random::seed();

  // State that we advance -- there's only one species.
  KMCState state(1);
  state[0] = 250;

  // Define list of reactions.
  std::vector<std::shared_ptr<const KMCReaction>> reactionList;

  auto c1R = std::make_shared<KMCReaction>(std::list<size_t>{0, 0}, std::list<size_t>{0, 0, 0});
  auto c2R = std::make_shared<KMCReaction>(std::list<size_t>{0, 0, 0}, std::list<size_t>{0, 0});
  auto c3R = std::make_shared<KMCReaction>(std::list<size_t>{}, std::list<size_t>{0});
  auto c4R = std::make_shared<KMCReaction>(std::list<size_t>{0}, std::list<size_t>{});

  c1R->rate() = c1 * N1; // Propensity becomes 1/2 * c1 * B1 * X * (X-1)
  c2R->rate() = c2;      // Propensity becomes 1/6 * c2 * X * (X-1) * (X-2)
  c3R->rate() = c3 * N2; // Propensity becomes c3 * B2
  c4R->rate() = c4;      // Propensity becomes c4 * X

  reactionList.emplace_back(c1R);
  reactionList.emplace_back(c2R);
  reactionList.emplace_back(c3R);
  reactionList.emplace_back(c4R);

  // Read input variables
  Real SSAlim   = 0.1;
  Real eps      = 0.3;
  int  numSteps = 100;
  int  numCrit  = 5;
  int  numSSA   = 5;

  std::string alg;

  pp.get("nsteps", numSteps);
  pp.get("algorithm", alg);
  pp.get("num_crit", numCrit);
  pp.get("num_ssa", numSSA);
  pp.get("eps", eps);
  pp.get("ssa_lim", SSAlim);

  // Define the Kinetic Monte Carlo solver and run it until the stop time.
  KMCSolverType kmcSolver(reactionList);

  kmcSolver.setSolverParameters(numCrit, numSSA, eps, SSAlim);

  // Run the SSA algorithm
  Real curDt    = 0.0;
  Real stopTime = 20.0;
  while (curDt < stopTime) {
    pout() << curDt << "\t" << state[0] << endl;

    Real nextDt = std::numeric_limits<Real>::max();
    if (alg == "ssa") {
      nextDt = kmcSolver.getCriticalTimeStep(state);

      kmcSolver.stepSSA(state);
    }
    else if (alg == "tau") {
      nextDt = stopTime / numSteps;

      kmcSolver.advanceTau(state, nextDt);
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
