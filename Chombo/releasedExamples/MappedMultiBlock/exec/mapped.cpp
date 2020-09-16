#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;

#include "parstream.H"
#include "NewMultiBlockFluxRegister.H"
#include "CubedSphereShellCS.H"
#include "RefCountedPtr.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "FArrayBox.H"
#include "DebugOut.H"
#include "DebugDump.H"
#include "FABView.H"

#include "UsingNamespace.H"



int
anisoCubedSphereShellTest()
{
  int boxsize = 16;
  int nlayers = (SpaceDim == 3) ? 8 : 1;
  RealVect dx = 1.0*RealVect::Unit;
  // Create a cubed sphere coordinate system with six  blocks.
  CubedSphereShellCS* css = new CubedSphereShellCS(boxsize,nlayers);
  // Cubed sphere shell has a specified size of index space 
  IntVect size(D_DECL6(11*boxsize, boxsize, nlayers, 1, 1, 1));
  Box domainBox(IntVect::Zero, size-IntVect::Unit);
  ProblemDomain coarDomain(domainBox);
  ProblemDomain fineDomain = coarDomain;
  fineDomain.refine(IntVect(D_DECL6(2,2,1,1,1,1)));
  css->define(coarDomain, dx);
  RefCountedPtr<MultiBlockCoordSys > coordSys(
      dynamic_cast<MultiBlockCoordSys* >(css));

  //create a coarse dbl over the entire grid
  Vector<int> coarProcs;
  Vector<Box> coarBoxes = coordSys->mappingBlocks();
  LoadBalance(coarProcs, coarBoxes);
  DisjointBoxLayout coarGrids(coarBoxes, coarProcs, coarDomain);

 
  Vector<Box> fineBoxes;
  Vector<int> fineProcs;
  // need to build AMR Hierarchy here  (bvs)

  LoadBalance(fineProcs, fineBoxes);
  DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);
  
  // Create the flux register.
  NewMultiBlockFluxRegister fluxReg(coordSys, fineGrids, coarGrids, 2, 1, true);
  fluxReg.setToZero();
  Real uval = 4;
  Real fval = 7;
  // Create a zero solution U the coarse grids with a single ghost layer.
  LevelData<FArrayBox> Ucoar(coarGrids, 1, IntVect::Unit);
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
    Ucoar[dit()].setVal(uval);
  
  // Create coarse and fine fluxes that straddle the coarse-fine boundary.
  LevelData<FluxBox> coarFlux(coarGrids, 1), fineFlux(fineGrids, 1);
  
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
    coarFlux[dit()].setVal(fval);
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    fineFlux[dit()].setVal(fval);
  
  // Stick the coarse and fine fluxes into the register to compute the defect.
  Interval interv(0, 0);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
            {
              fluxReg.incrementCoarse(coarFlux[dit()][idir], 1.0, dit(), interv, interv, idir, sit());
            }
          for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
            {
              fluxReg.incrementFine(fineFlux[dit()][idir], 1.0, dit(), interv, interv, idir, sit());
            }
        }
    }
  
  // Apply the defect to the solution.
  fluxReg.poutCoarseRegisters();
  fluxReg.poutFineRegisters();
  fluxReg.reflux(Ucoar, interv, 1.0);
  
  
  
  // Git.
  return 0;
}


//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  bool justGetOut = false;
#ifdef CH_USE_FLOAT
  justGetOut = true;
#endif
  if (SpaceDim==1) 
  {
    justGetOut = true;
  }

 int rank=0, number_procs=1;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#endif

  if (rank == 0)
    {
      pout() << " number_procs = " << number_procs << endl;
    }

  // Check for an input file
  char* inFile = NULL;

  if (argc > 1)
    {
      inFile = argv[1];
    }
  else
    {
      pout() << "Usage:  mapped....ex <inputfile>" << endl;
      pout() << "No input file specified" << endl;
      return -1;
    }

  // Parse the command line and the input file (if any)
  ParmParse pp(argc-2,argv+2,NULL,inFile);


  int overallStatus = 0;
  int status=0;

  if (!justGetOut)
  {
 
 
    pout () << "Beginning mapped " << " ..." << endl ;


    status = anisoCubedSphereShellTest();
    if ( status == 0 )
    {
      pout() << ": mapped passed." << std::endl ;
    }
    else
    {
      pout() << 
        ": mapped failed with return code " << status << std::endl ;
      overallStatus += 1;
    }
  }

#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return overallStatus;
}
//-----------------------------------------------------------------------
