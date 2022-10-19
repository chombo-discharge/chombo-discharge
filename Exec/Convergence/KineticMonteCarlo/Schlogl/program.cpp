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

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Seed the RNG
  Random::setRandomSeed();

  // State that we advance -- there's only one species.
  KMCSingleState<> state(1);
  state[0] = 250;

  // Define list of reactions.
  std::vector<std::shared_ptr<const KMCSingleStateReaction<>>> reactionList;

  auto c1R = std::make_shared<KMCSingleStateReaction<>>(std::list<size_t>{0, 0}, std::list<size_t>{0, 0, 0});
  auto c2R = std::make_shared<KMCSingleStateReaction<>>(std::list<size_t>{0, 0, 0}, std::list<size_t>{0, 0});
  auto c3R = std::make_shared<KMCSingleStateReaction<>>(std::list<size_t>{}, std::list<size_t>{0});
  auto c4R = std::make_shared<KMCSingleStateReaction<>>(std::list<size_t>{0}, std::list<size_t>{});

  c1R->rate() = c1 * N1; // Propensity becomes 1/2 * c1 * B1 * X * (X-1)
  c2R->rate() = c2;      // Propensity becomes 1/6 * c2 * X * (X-1) * (X-2)
  c3R->rate() = c3 * N2; // Propensity becomes c3 * B2
  c4R->rate() = c4;      // Propensity becomes c4 * X

  reactionList.emplace_back(c1R);
  reactionList.emplace_back(c2R);
  reactionList.emplace_back(c3R);
  reactionList.emplace_back(c4R);

  // Define the Kinetic Monte Carlo solver and run it until time = 10.
  KMCSolver<KMCSingleStateReaction<>, KMCSingleState<>, long long> kmcSolver(reactionList);

  // Run the SSA algorithm
  Real curDt  = 0.0;
  Real stopDt = 20.0;
  while (curDt < stopDt) {
    const Real nextDt = kmcSolver.getCriticalTimeStep(state);

    pout() << curDt << "\t" << nextDt << "\t" << state[0] << endl;

    kmcSolver.stepSSA(state);

    curDt += nextDt;
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
