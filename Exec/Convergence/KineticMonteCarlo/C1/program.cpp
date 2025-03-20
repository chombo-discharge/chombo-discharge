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

  Random::setRandomSeed();

  // First, generate solutions using the SSA.
  Real exactMean = 0.0;
  Real exactVar  = 0.0;

  // Advance problem using non-exact sampling algorithms. The solution with the finest time step is used as a
  // proxy for the "exact" answer
  Vector<Real> algMean(numSteps.size(), 0.0);
  Vector<Real> algVar(numSteps.size(), 0.0);

  for (int istep = numSteps.size() - 1; istep >= 0; istep--) {
    Vector<Real> algSoln;

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
          nextDt = kmcSolver.getCriticalTimeStep(state);

          if (nextDt < stopTime - curTime) {
            kmcSolver.stepSSA(state);
          }
          else {
            nextDt = stopTime - curTime;
          }
        }
        else if (alg == "tau_plain") {
          nextDt = stopTime / numSteps[istep];

          kmcSolver.advanceTauPlain(state, nextDt);
        }
        else if (alg == "tau_midpoint") {
          nextDt = stopTime / numSteps[istep];

          kmcSolver.advanceTauMidpoint(state, nextDt);
        }
        else if (alg == "tau_prc") {
          nextDt = stopTime / numSteps[istep];

          kmcSolver.advanceTauPRC(state, nextDt);
        }
        else if (alg == "hybrid_plain") {
          nextDt = stopTime / numSteps[istep];

          kmcSolver.advanceHybrid(state, nextDt, KMCLeapPropagator::TauPlain);
        }
        else if (alg == "hybrid_midpoint") {
          nextDt = stopTime / numSteps[istep];

          kmcSolver.advanceHybrid(state, nextDt, KMCLeapPropagator::TauMidpoint);
        }
        else if (alg == "hybrid_prc") {
          nextDt = stopTime / numSteps[istep];

          kmcSolver.advanceHybrid(state, nextDt, KMCLeapPropagator::TauPRC);
        }
        else {
          const std::string err = "Don't know the algoritm '" + alg + "'";

          MayDay::Error(err.c_str());
        }

        curTime += nextDt;
      }

      if (irun == numRuns - 1 && istep == numSteps.size() - 1) {
        pout() << curTime << "\t" << state[0] << "\t" << 1.0 * initVal * exp(effRate * curTime) << endl;
      }

      algSoln.push_back(1.0 * state[0]);
    }

    // Compute the mean value.
    algMean[istep] = 0.0;
    for (int irun = 0; irun < numRuns; irun++) {
      algMean[istep] += algSoln[irun];
    }
#if CH_MPI
    algMean[istep] = ParallelOps::sum(algMean[istep]);
    algMean[istep] /= (numProc() * numRuns);
#else
    algMean[istep] /= numRuns;
#endif

    // Compute the variance -- i.e. the deviation from the "exact" solution.
    algVar[istep] = 0.0;
    for (int irun = 0; irun < numRuns; irun++) {
      algVar[istep] += std::pow(algSoln[irun] - algMean.back(), 2);
    }
#if CH_MPI
    algVar[istep] = ParallelOps::sum(algVar[istep]);
    algVar[istep] /= (numProc() * numRuns);
#else
    algVar[istep] /= numRuns;
#endif
    algVar[istep] = sqrt(algVar[istep]);
  }

  exactMean = algMean.back();
  exactVar  = algVar.back();

  // Print the mean solution error and convergence rate. Oh, and this only makes sense for non-SSA runs
  // since the SSA algorithm does the time steps differently.
  if (procID() == 0 && alg != "ssa") {
    // clang-format off
    std::cout << std::left << std::setw(20) << "# dt"
	      << std::left << std::setw(20) << "num steps"
	      << std::left << std::setw(20) << "Mean"      
	      << std::left << std::setw(20) << "Rel. mean error"
	      << std::left << std::setw(20) << "Rel. std error"
	      << std::left << std::setw(20) << "Mean conv. order"
	      << std::left << std::setw(20) << "Var conv. order"
	      << "\n";
    // clang-format on

    for (int istep = 0; istep < numSteps.size(); istep++) {
      const Real curDt   = stopTime / numSteps[istep];
      const Real meanErr = std::abs(algMean[istep] - exactMean);
      const Real varErr  = std::abs(algVar[istep] - exactVar);

      if (istep < numSteps.size() - 1) {

        const Real finerErr = std::abs(algMean[istep + 1] - exactMean);
        const Real finerVar = std::abs(algVar[istep + 1] - exactVar);
        const Real finerDt  = stopTime / numSteps[istep + 1];
        const Real meanConv = log(meanErr / finerErr) / log(curDt / finerDt);
        const Real varConv  = log(varErr / finerVar) / log(curDt / finerDt);

        // clang-format off
	std::cout << std::left << std::setw(20) << curDt
		  << std::left << std::setw(20) << numSteps[istep]
		  << std::left << std::setw(20) << algMean[istep]
		  << std::left << std::setw(20) << meanErr/exactMean
		  << std::left << std::setw(20) << varErr/exactMean
		  << std::left << std::setw(20) << meanConv
		  << std::left << std::setw(20) << varConv
		  << "\n";
        // clang-format on
      }
      else {
        // clang-format off
	std::cout << std::left << std::setw(20) << curDt
		  << std::left << std::setw(20) << numSteps[istep]
		  << std::left << std::setw(20) << algMean[istep]	  
		  << std::left << std::setw(20) << meanErr
		  << std::left << std::setw(20) << varErr
		  << std::left << std::setw(20) << "*"
		  << std::left << std::setw(20) << "*"	  
		  << "\n";
        // clang-format on
      }
    }
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
