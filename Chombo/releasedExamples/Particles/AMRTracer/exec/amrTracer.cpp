#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "AMR.H"
#include "AMRLevelTracer.H"
#include "AdvectionFunctions.H"
#include "DebugDump.H"
#include "memtrack.H"
#include "CH_Attach.H"
#include "FABView.H"
#include "AMRTracerUtils.H"
#include "parstream.H"

#include "UsingNamespace.H"

/***************/
void amrTest(const Real& a_stopTime,
             const int&  a_nstop,
	           const Vector<int>& a_refRat)
{
  // read inputs
  ParmParse pp;

  ProblemDomain prob_domain;
  getProblemDomain(prob_domain);

  AdvectionVelocityFunction velFunc;
  getAdvectionVelocityFunction(velFunc);

  RefCountedPtr<AMRLevelTracerFactory>  amrt_fact;
  getAMRLevelTracerFactory(amrt_fact, velFunc);

  AMR amr;
  defineAMR(amr, amrt_fact, prob_domain, a_refRat);
  
  Vector<Vector<Box> > amrHierarchy;
  setupAMRHierarchy(amrHierarchy);

  setupAMRForFixedHierarchyRun(amr, amrHierarchy);

  // run
  amr.run(a_stopTime, a_nstop);

  // output last pltfile and statistics
  //cleanup
  amr.conclude();
  
}
/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif
  { //scoping trick

    // Check for an input file
    char* inFile = NULL;

    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
    // Parse the command line and the input file (if any)
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
    Real stopTime = 0.0;
    pp.get("max_time",stopTime);

    int nstop = 0;
    pp.get("max_step",nstop);
    int max_level = 0;
    pp.get("max_level",max_level);
    int num_read_levels = Max(max_level,1);
    Vector<int> ref_ratios; // (num_read_levels,1);
    pp.getarr("ref_ratio",ref_ratios,0,num_read_levels+1);
    amrTest(stopTime, nstop, ref_ratios);
  }
#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
