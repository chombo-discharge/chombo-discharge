#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include  <iostream>
#include  <cstdio>

#include "LevelOp.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "PoissonBC.H"
#include "ParmParse.H"
#include "AMRSolver.H"
#include "PoissonOp.H"
#include "Vector.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "IntVectSet.H"
#include "LoadBalance.H"
#include "localfunctions.H"
#include "fortranfunctions.H"
#include "DebugOut.H"
//#include "BoxIterator.H"
//#include "LayoutIterator.H"
#include "BRMeshRefine.H"


#include "LocalTimer.H"


LocalTimer Everything("Everything");
LocalTimer InitGrids ("Init Grids",  Everything);
LocalTimer Refine    ("Refine    ",  Everything);
LocalTimer AMRSetup  ("AMR Setup ",  Everything);
LocalTimer AMRSolve  ("AMR Solve ",  Everything);
LocalTimer Output    ("Output    ",  Everything);

int main(int argc, char* argv[])
{
  int SIZE = 64;
  bool PP = false;
  if (argc > 1 )
  {
    SIZE = atoi(argv[1]);
    if (SIZE < 0)
    {
      PP = true;
      SIZE = abs(SIZE);
    }
  }

  int rank, number_procs;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &number_procs);
#else
  rank=0;
  number_procs=1;
#endif
  cout << "rank = " << rank << " number_procs= " << number_procs << endl;
  { // parallel scope

  Everything.start();

  InitGrids.start();
#if CH_SPACEDIM == 2
  const int finestLevel  = 1;
  const int Ncells[2] =
  {
    SIZE, SIZE
  };
#else
  const int finestLevel  = 4;
  const int Ncells[3] =
  {
    8,8,8
  };
#endif

  const int iverbose            = 1;
  const int maxIterations       = 10;
  const int numberVcyclesBottom = 1;
  const int blockFactor         = 8;
  const int maxBoxSize          = 32;
  const int bufferSize          = 1;
  const int baseLevel           = 0;

  const int ibclo[3] =
  {
    0,0,0
  };
  const int ibchi[3] =
  {
    0,0,0
  };

  const Real fillRatio          = 0.80;
  const Real tagThreshold       = 1.0e6;
  //  const Real tagThreshold       = 0.005;
  const Real ProblemLowCord[3]  =
  {
    0.0, 0.0, 0.0
  };
  const Real ProblemHighCord[3] =
  {
    1.0, 1.0, 1.0
  };

  int numLevels;

  bool verbose = (iverbose == 1);
  //const IntVect unitVector(IntVect::Unit);

  Vector<int> refinementRatioVector;
  refinementRatioVector.resize(finestLevel+1);
  refinementRatioVector[0] = 2;
  for (int i=1; i<=finestLevel; i++)
  {
    refinementRatioVector[i] = 2;
  }

  DomainGhostBC domghostbc;
  setDomainBC(domghostbc, ibclo, ibchi, verbose);
  cout << "DomainGhostBC domghostbc built" << endl;

  PoissonOp poissonOperator;
  poissonOperator.setDomainGhostBC(domghostbc);
  cout << "PoissonOp poissonOperator built" << endl;


  Vector<LDFB*>                rhs(finestLevel+1, NULL);
  Vector<LDFB*>                phi(finestLevel+1, NULL);
  Vector<DBL>           DBLVector(finestLevel+1);
  Vector<Box>   boundingBoxVector(finestLevel+1);
  Vector<Real>           DxVector(finestLevel+1);
  Vector<Vector<Box> >   amrBoxes(finestLevel+1);
  Vector<Vector<int> >    procIDs(finestLevel+1);



  // make one big box with Ncells dimensions
  const Box baseDomain = createBox(Ncells, SpaceDim);
  boundingBoxVector[0] = baseDomain;
  cout << " base domain -- boundingBoxVector[0]= " << baseDomain << endl;
  DxVector[0] = (ProblemHighCord[0] - ProblemLowCord[0])/Ncells[0];

  // define base level first
  domainSplit(boundingBoxVector[0], amrBoxes[0], maxBoxSize, blockFactor);

  for (int i=1; i<=finestLevel; i++)
  {
    boundingBoxVector[i] = refine(boundingBoxVector[i-1],
                                  refinementRatioVector[i-1]);
    CH_assert(refinementRatioVector[i-1] > 0);
    DxVector[i] = DxVector[i-1]/refinementRatioVector[i-1];
  }

  procIDs[0].resize(amrBoxes[0].size());
  LoadBalance(procIDs[0], amrBoxes[0]);

  cout << " amrBoxes[0].size() = " << amrBoxes[0].size() << endl;
  //cout << " procIDs[0]= "             << procIDs[0] << endl;
  cout << " procIDs[0].size() = "     << procIDs[0].size() << endl;

  DBLVector[0].define(amrBoxes[0], procIDs[0]);

  rhs[0] = new LDFB(DBLVector[0], 1, unitVector);

  // now initialize data on this level
  LDFB& levelrhs = *rhs[0];
  for (DataIterator dit = levelrhs.dataIterator(); dit.ok(); ++dit)
  {
    FORT_SETRHS(CHF_FRA( levelrhs[dit()]),
                CHF_BOX( levelrhs[dit()].box()),
                CHF_CONST_REAL(DxVector[0]));
  }
  InitGrids.stop();

  Refine.start();
  BRMeshRefine meshRefine(boundingBoxVector[0], refinementRatioVector,
                          fillRatio, blockFactor, bufferSize, maxBoxSize);
  Refine.stop();


  AMRSolver amrSolver(DBLVector,  boundingBoxVector,
                      DxVector,   refinementRatioVector,
                      1, baseLevel, &poissonOperator);
  amrSolver.setVerbose(verbose);
  amrSolver.setMaxIter(maxIterations);
  amrSolver.setNumVCyclesBottom(numberVcyclesBottom);

  if (phi[0] != NULL) delete phi[0];
  phi[0] = new LDFB(DBLVector[0], 1, unitVector);
  amrSolver.solveAMR(phi, rhs);

  cout << " Finished Level 0! " << endl;

  // already have one level defined above, so current level is 1
  int clev = 1;
  bool found_tags=false;

  while (clev <= finestLevel)
  {

     printf("##############################################\n");
     printf("#####  Level %2d                          #####\n", clev);
     printf("##############################################\n");

    Vector<IntVectSet> tagVect(clev+1);
    newTagCells(rhs, phi, clev, DBLVector, boundingBoxVector, DxVector,
                refinementRatioVector, domghostbc, baseLevel, verbose,
                maxIterations, numberVcyclesBottom, tagThreshold, tagVect);

    Vector<Vector<Box> > newamrBoxes(clev+1);

    int newFinestLevel;
    if (procID() == uniqueProc(SerialTask::compute) )
    {
      // Do we need to generate a new level?
      //printf(" clev = %d  empty=%d\n", clev, tagVect[clev-1].isEmpty());
      if (clev >= 0 && !tagVect[clev-1].isEmpty())
      {
        cout << "   found tags, meshRefine... " << endl;
        newFinestLevel = meshRefine.regrid(newamrBoxes, // output
                                           tagVect,  // input, but gets changed!
                                           baseLevel, clev-1,
                                           amrBoxes);
        printf("   after meshRefine newFinestLevel=%d\n", newFinestLevel);

        if (newFinestLevel == clev)
        {
          found_tags = true;
          clev++;
        }
        else
        {
          printf(" regrid returned a newFinestLevel that wasn't kosher -- breaking\n");
          found_tags = false;
        }

      }
      else
      {
        printf(" NO TAGS found, grid resolved, breaking.\n");
        found_tags = false;
      }
    } // end if proc is serial node.

    if (!found_tags) break;

    broadcast(clev,        uniqueProc(SerialTask::compute) );
    broadcast(newamrBoxes, uniqueProc(SerialTask::compute) );

    amrBoxes = newamrBoxes;

    for (int ilev=0; ilev<clev; ilev++)
    {
      procIDs[ilev].resize(amrBoxes[ilev].size());
      LoadBalance(procIDs[ilev], amrBoxes[ilev]);

      //cout << " currentLevel = " << clev << " ilev=" << ilev << endl;
      //cout << "   amrBoxes[ilev].size() = " << amrBoxes[ilev].size() << endl;
      //cout << " procIDs[ilev]= "             << procIDs[ilev] << endl;
      //cout << "   procIDs[ilev].size()  = "  << procIDs[ilev].size() << endl;

      //DBLVector[ilev].define(amrBoxes[ilev], procIDs[ilev]);
      const DBL newDBL(amrBoxes[ilev], procIDs[ilev]);
      DBLVector[ilev] = newDBL;

      if (rhs[ilev] != NULL) delete rhs[ilev];
      rhs[ilev] = new LDFB(DBLVector[ilev], 1, unitVector);

      // now initialize data on this level
      LDFB& levelrhs = *rhs[ilev];
      for (DataIterator dit = levelrhs.dataIterator(); dit.ok(); ++dit)
      {
        FORT_SETRHS(CHF_FRA( levelrhs[dit()]),
                    CHF_BOX( levelrhs[dit()].box()),
                    CHF_CONST_REAL(DxVector[ilev]));
      }
    }

    AMRSolver amrSolver(DBLVector,  boundingBoxVector,
                        DxVector,   refinementRatioVector,
                        clev, baseLevel, &poissonOperator);
    amrSolver.setVerbose(verbose);
    amrSolver.setMaxIter(maxIterations);
    amrSolver.setNumVCyclesBottom(numberVcyclesBottom);

    for (int ilev = 0; ilev < clev; ilev++)
    {
      if (phi[ilev] != NULL) delete phi[ilev];
      phi[ilev] = new LDFB(DBLVector[ilev], 1, unitVector);
    }
    amrSolver.solveAMR(phi, rhs);


  }  // end while on levels...  finestLevel

  cout << " Done!  clev=" << clev << endl;

  numLevels = clev;

  if (found_tags)
  {
    printf(" found_tags = true \n");
  }


  // make an exact HOG (since we know what it is)
  Vector<LDFB*> exact(numLevels, NULL);
  for (int ilev=0; ilev<numLevels; ilev++)
  {
    exact[ilev] = new LDFB(DBLVector[ilev], 1, unitVector);
    LDFB& level_exact = *exact[ilev];
    for (DataIterator dit = level_exact.dataIterator(); dit.ok(); ++dit)
    {
      FORT_SETEXACT(CHF_FRA( level_exact[dit()]),
                    CHF_BOX( level_exact[dit()].box()),
                    CHF_CONST_REAL(DxVector[ilev]));
    }
  }
  outputHDF5(exact, DBLVector, boundingBoxVector,
             refinementRatioVector, numLevels, "exact.hdf5");


  // subtract  phi from exact and put in error
  Vector<LDFB*> error(numLevels, NULL);
  createHOG(error, numLevels, DBLVector);

  // error = exact-phi
  subtractTwoHOGs(exact, phi, numLevels, error);

  outputHDF5(error, DBLVector, boundingBoxVector,
             refinementRatioVector, numLevels, "error.hdf5");

  //outputHDF5(exact, DBLVector, boundingBoxVector,
  //         refinementRatioVector, numLevels, "exact2.hdf5");

  // Compute the Laplacian of the error
  AMRSolver amrSolverError(DBLVector,  boundingBoxVector,
                           DxVector,   refinementRatioVector,
                           numLevels, baseLevel, &poissonOperator);

  Vector<LDFB*> lap_error(numLevels, NULL);
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    lap_error[ilev] = new LDFB(DBLVector[ilev], 1, unitVector);
    amrSolverError.applyAMROperator(error, *lap_error[ilev], ilev);
  }

  multiplyHOGbyConstant(lap_error, numLevels, 3.0);

  outputHDF5(lap_error, DBLVector, boundingBoxVector,
             refinementRatioVector, numLevels, "lap_error.hdf5");


  // Compute the Laplacian of the exact
  AMRSolver amrSolverExact(DBLVector,  boundingBoxVector,
                           DxVector,   refinementRatioVector,
                           numLevels, baseLevel, &poissonOperator);

  Vector<LDFB*> lap_exact(numLevels, NULL);
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    lap_exact[ilev] = new LDFB(DBLVector[ilev], 1, unitVector);
    amrSolverExact.applyAMROperator(exact, *lap_exact[ilev], ilev);
  }

  outputHDF5(lap_exact, DBLVector, boundingBoxVector,
             refinementRatioVector, numLevels, "lap_exact.hdf5");


  outputHDF5(rhs, DBLVector, boundingBoxVector,
             refinementRatioVector, numLevels, "rhs.hdf5");
  outputHDF5(phi, DBLVector, boundingBoxVector,
             refinementRatioVector, numLevels, "phi.hdf5");

  //prettyPrintLDFB(rhs[0], "rhs0");
  //prettyPrintLDFB(rhs[1], "rhs1");

  for (int ilev = 0; ilev < numLevels-1; ilev++)
  {
    delete phi[ilev];
    delete rhs[ilev];
    //delete rhs_lap[ilev];
  }


  Everything.stop();

  LocalTimer::LocalTimerSummary();

  } // end parallel scope
#ifdef CH_MPI
  //dumpmemoryatexit();
  MPI_Finalize();
#endif

  return(0);
}

Box createBox(const int *Ncells, const int dimensions)
{

  IntVect ivlo = IntVect::Zero;
  IntVect ivhi;
  for (int idir = 0; idir < dimensions; idir++)
  {
    ivhi[idir] = Ncells[idir] - 1;
  }
  return Box(ivlo, ivhi);
}


