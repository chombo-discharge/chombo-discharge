#include <CD_Driver.H>
#include <CD_KMCSolver.H>
#include <CD_KMCDualStateReaction.H>

using namespace ChomboDischarge;

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Initial state.
  KMCDualState<> state(1, 0);

  // One initial electron.
  state.getReactiveState()[0] = 1;

  // Define e + null -> e + e + null
  auto ionization = std::make_shared<KMCDualStateReaction<>>(std::list<size_t>{0},
                                                             std::list<size_t>{0, 0},
                                                             std::list<size_t>{});

  // Define e -> null
  auto attachment = std::make_shared<KMCDualStateReaction<>>(std::list<size_t>{0},
                                                             std::list<size_t>{},
                                                             std::list<size_t>{});

  // Set ionization and attachment rates.
  ionization->rate() = 1.0;
  attachment->rate() = 1.0;

  // Define list of reactions.
  std::vector<std::shared_ptr<const KMCDualStateReaction<>>> reactionList;
  reactionList.emplace_back(ionization);
  reactionList.emplace_back(attachment);

  // Define the Kinetic Monte Carlo solver and run it until time = 10.
  KMCSolver<KMCDualStateReaction<>, KMCDualState<>, long long> kmcSolver(reactionList);

  kmcSolver.advanceHybrid(state, 10.0);

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
