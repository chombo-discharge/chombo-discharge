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

  // Read input file.
  ParmParse pp(argc - 2, argv + 2, NULL, argv[1]);

  // Seed the RNG
  Random::setRandomSeed();

  // State that we advance.
  KMCDualState<> state(1, 0);

  // Define list of reactions.
  std::vector<std::shared_ptr<const KMCDualStateReaction<>>> reactionList;

  // Define reaction e + null -> e + e + null
  auto ionization = std::make_shared<KMCDualStateReaction<>>(std::list<size_t>{0},
                                                             std::list<size_t>{0, 0},
                                                             std::list<size_t>{});

  // Define e -> null
  auto attachment = std::make_shared<KMCDualStateReaction<>>(std::list<size_t>{0},
                                                             std::list<size_t>{},
                                                             std::list<size_t>{});

  reactionList.emplace_back(attachment);
  reactionList.emplace_back(ionization);

  // Read initial values.
  Real& ionizationRate = ionization->rate();
  Real& attachmentRate = attachment->rate();

  Real stopDt = 0.0;
  int  initVal;

  // Read program input variables
  pp.get("stop_time", stopDt);
  pp.get("ionization_rate", ionizationRate);
  pp.get("attachment_rate", attachmentRate);
  pp.get("initial_particles", initVal);

  // Set initial number of particles
  state.getReactiveState()[0] = (long long)initVal;

  // Define the Kinetic Monte Carlo solver and run it until time = 10.
  KMCSolver<KMCDualStateReaction<>, KMCDualState<>, long long> kmcSolver(reactionList);

  // Run the SSA algorithm
  Real curDt = 0.0;
  while (curDt < stopDt) {
    pout() << curDt << "\t" << state.getReactiveState()[0] << endl;

    const Real nextDt = kmcSolver.getCriticalTimeStep(state);

    kmcSolver.stepSSA(state);

    curDt += nextDt;
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
