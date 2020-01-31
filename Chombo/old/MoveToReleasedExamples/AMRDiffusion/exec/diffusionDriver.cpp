#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MayDay.H"
#include "ParmParse.H"

#include "DiffusionParams.H"
#include "DiffusionSolver.H"
#include "DebugDump.H"
#include "FABView.H"

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
    // Needs to be called with an input file
    if (argc < 2)
    {
      MayDay::Abort("Usage: diffusionDriver input_file");
    }

    // Read command line and input file
    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    // Set and output the parameters for the simulation
    DiffusionParams params;
    params.print();

    // Construct, initialize, and run the solver
    DiffusionSolver solver(params);
    solver.init();
    solver.run();
  }
  // End scoping trick

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif

  exit(0);
}
