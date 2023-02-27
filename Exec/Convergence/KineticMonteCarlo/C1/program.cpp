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
  int  numRuns  = 1;
  int  initVal  = 0;
  int  numCrit  = 5;
  int  numSSA   = 5;

  Vector<int> numSteps;

  std::string alg;

  // Read program input variables
  pp.get("nruns", numRuns);
  pp.get("stop_time", stopTime);
  pp.get("ionization_rate", ionizationRate);
  pp.get("attachment_rate", attachmentRate);
  pp.get("initial_particles", initVal);
  pp.get("algorithm", alg);
  pp.get("num_crit", numCrit);
  pp.get("num_ssa", numSSA);
  pp.get("eps", eps);
  pp.get("ssa_lim", SSAlim);
  pp.getarr("nsteps", numSteps, 0, pp.countval("nsteps"));

  // Define the Kinetic Monte Carlo solver and run it until time = 10.
  KMCSolverType kmcSolver(reactionList);

  kmcSolver.setSolverParameters(numCrit, numSSA, eps, SSAlim);

  const Real effRate = ionizationRate - attachmentRate;

  Vector<Real> mean(numSteps.size(), 0.0);

  Random::setRandomSeed();

  // Advance problem using different step sizes and population samples
  for (int istep = 0; istep < numSteps.size(); istep++) {
    mean[istep] = 0.0;

    for (int irun = 0; irun < numRuns; irun++) {
      state[0] = (long long)initVal;

      Real curTime = 0.0;

      while (curTime < stopTime) {

        // Print last state.
        if (irun == numRuns - 1 && istep == numSteps.size() - 1) {
          pout() << curTime << "\t" << state[0] << "\t" << 1.0 * initVal * exp(effRate * curTime) << endl;
        }

        Real nextDt = std::numeric_limits<Real>::max();

        if (alg == "ssa") {
          nextDt = std::min(stopTime - curTime, kmcSolver.getCriticalTimeStep(state));

          kmcSolver.stepSSA(state);
        }
        else if (alg == "tau") {
          nextDt = std::min(stopTime - curTime, stopTime / (numSteps[istep] - 1));

          kmcSolver.advanceTau(state, nextDt);
        }
        else if (alg == "heun") {
          nextDt = std::min(stopTime - curTime, stopTime / (numSteps[istep] - 1));

          kmcSolver.advanceHeun(state, nextDt);
        }
        else if (alg == "hybrid") {
          nextDt = std::min(stopTime - curTime, stopTime / (numSteps[istep] - 1));

          kmcSolver.advanceHybrid(state, nextDt);
        }
        else {
          const std::string err = "Expected algorithm to be 'ssa', 'tau', 'heun', or 'hybrid' but got '" + alg + "'";

          MayDay::Error(err.c_str());
        }

        curTime += nextDt;
      }

      if (irun == numRuns - 1 && istep == numSteps.size() - 1) {
        pout() << curTime << "\t" << state[0] << "\t" << 1.0 * initVal * exp(effRate * curTime) << endl;
      }

      mean[istep] += state[0];
    }

#if CH_MPI
    mean[istep] = ParallelOps::sum(mean[istep]);
    mean[istep] /= (numProc() * numRuns);
#else
    mean[istep] /= numRuns;
#endif
  }

  // Print the mean solution error and convergence rate. Oh, and this only makes sense for non-SSA runs
  // since the SSA algorithm does the time steps differently.
  if (procID() == 0 && alg != "ssa" && numSteps.size() > 1) {
    const Real exactSoln = 1.0 * initVal * exp(effRate * stopTime);

    std::cout << std::left << std::setw(15) << "dt" << std::left << std::setw(15) << "Rel. error" << std::left
              << std::setw(15) << "Richardson p" << std::endl;

    for (int istep = 0; istep < numSteps.size(); istep++) {
      const Real curDt  = stopTime / (numSteps[istep] - 1);
      const Real curErr = std::abs(mean[istep] - exactSoln);

      if (istep < numSteps.size() - 1) {
        const Real finerDt  = stopTime / (numSteps[istep + 1] - 1);
        const Real finerErr = std::abs(mean[istep + 1] - exactSoln);
        const Real p        = log(curErr / finerErr) / log(curDt / finerDt);

        std::cout << std::left << std::setw(15) << curDt << std::left << std::setw(15) << curErr / exactSoln
                  << std::left << std::setw(15) << p << "\n";
      }
      else {
        std::cout << std::left << std::setw(15) << curDt << std::left << std::setw(15) << curErr / exactSoln
                  << std::left << std::setw(15) << "*"
                  << "\n";
      }
    }
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
