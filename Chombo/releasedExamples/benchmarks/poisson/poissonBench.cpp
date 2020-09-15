#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cerr;

#include "LevelData.H"
#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "Vector.H"
#include "RealVect.H"
#include "CH_Timer.H"
#include "CH_Attach.H"
#include "memusage.H"
#include "memtrack.H"
#include "PoissonF_F.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include <string>
#include <iostream>
#include <sstream>

/******/
int runBench()
{
  //constants

  int nx = 64;
  int maxbox = 32;
  int blocfactor = 8;
  int nghostPhi  = 1;
  int nghostLph  = 0;
  Real domlen = 1.0;
  Real alpha = 0.0;
  Real beta  = 1.0;
  int numApplies = 100;
  //derived 
  Real dx = domlen/nx;
  RealVect dxVect = dx*RealVect::Unit;

  IntVect ivlo = IntVect::Zero;
  IntVect ivhi = (nx-1)*IntVect::Unit;
  ProblemDomain domain(ivlo, ivhi);
  Vector<Box> boxes;
  Vector<int> procs;

  domainSplit(domain, boxes, maxbox, blocfactor);
  
  LoadBalance(procs, boxes);
  DisjointBoxLayout  dbl(boxes, procs, domain);
  LevelData<FArrayBox> phi(dbl, 1, nghostPhi*IntVect::Unit);
  LevelData<FArrayBox> lph(dbl, 1, nghostLph*IntVect::Unit);
  for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
  {
    phi[dit()].setVal(1.0);
    lph[dit()].setVal(1.0);
  }    

  for(int iapply = 0; iapply < numApplies; iapply++)
  {
    CH_TIME("apply_only");
    for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      FORT_SIMPLESTLAPLACIAN(CHF_FRA1(lph[dit()],0),
                             CHF_CONST_FRA1(phi[dit()],0),
                             CHF_CONST_REAL(alpha),
                             CHF_CONST_REAL(beta),
                             CHF_CONST_REALVECT(dxVect),
                             CHF_BOX(dbl[dit()]));
    }
  }
  return 0;
}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;

  // scoping...
  {
    int solverStatus = runBench();
    status += solverStatus;
  }
  //end scoping trick
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}
