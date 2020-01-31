#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// This code reads a plotfile and writes another.
// The output file contains a subset of the input file defined by a index
// range on the coarsest level.

#include <iostream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"

// One more function for MPI
void dumpmemoryatexit();

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  // setChomboMPIErrorHandler();
  MPI_Barrier(Chombo_MPI::comm);  // Barrier #1
#endif

  // ------------------------------------------------
  // parse command line and input file
  // ------------------------------------------------
  // infile must be first
  if (argc < 2)
    {
      cerr << "  need inputs file" << endl;
      abort();
    }

  char* in_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, in_file);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  pout() << endl;
  pout() << "Initializing..." << endl;

  // declare variable to store hierarchy
  string                         inFileName;
  Vector<DisjointBoxLayout>      inGrids;
  Vector<LevelData<FArrayBox>* > inData;
  Vector<string>                 inVars;
  Box                            inDomain;
  Real                           inDx;
  Real                           inDt;
  Real                           inTime;
  Vector<int>                    inRefRatio;
  int                            inNumLevels;

  pp.get("infile" ,inFileName);

  string outFileName;
  pp.get("outfile",outFileName);

  Vector<int> loEndVect;
  pp.getarr("lo_end",loEndVect,0,SpaceDim);

  Vector<int> hiEndVect;
  pp.getarr("hi_end",hiEndVect,0,SpaceDim);

  int outNumLevels;
  pp.get("num_levels",outNumLevels);

  IntVect loEnd;
  IntVect hiEnd;

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    loEnd[idir] = loEndVect[idir];
    hiEnd[idir] = hiEndVect[idir];
  }

  Box boxOfInterest(loEnd,hiEnd);

  pout() << "Reading data..." << endl;

  ReadAMRHierarchyHDF5(inFileName,
                       inGrids,
                       inData,
                       inVars,
                       inDomain,
                       inDx,
                       inDt,
                       inTime,
                       inRefRatio,
                       inNumLevels);

  Vector<LevelData<FArrayBox>* > outData(outNumLevels);
  Vector<DisjointBoxLayout>      outGrids(outNumLevels);
  Vector<int>                    outRefRatio(outNumLevels);

  Box outDomain = inDomain & boxOfInterest;

  // allocate output data -- same domain as input
  for (int level = 0; level < outNumLevels; level++)
    {
      const DisjointBoxLayout curDBL = inGrids[level];

      pout() << "Processing level " << level << "..." << endl;

      Vector<Box> outBoxes;
      Vector<int> outProcs;

      pout() << "  Creating data holders..." << endl;

      for (LayoutIterator lit = curDBL.layoutIterator(); lit.ok(); ++lit)
      {
        Box curBox = curDBL.get(lit());

        if (curBox.intersects(boxOfInterest))
        {
          outBoxes.push_back(curBox & boxOfInterest);
        }
      }

      LoadBalance(outProcs,outBoxes);

      outGrids[level].define(outBoxes,outProcs,outDomain);

      outData[level] = new LevelData<FArrayBox>(outGrids[level],
                                                inVars.size(),
                                                inData[level]->ghostVect());

      pout() << "  Copying data..." << endl;

      // copy data for this level
      Interval inInterval(0,inVars.size());
      inData[level]->copyTo(*outData[level]);

      outRefRatio[level] = inRefRatio[level];

      boxOfInterest.refine(inRefRatio[level]);
    }

  pout() << "Writing data..." << endl;

  WriteAMRHierarchyHDF5(outFileName,
                        outGrids,
                        outData,
                        inVars,
                        outDomain,
                        inDx,
                        inDt,
                        inTime,
                        outRefRatio,
                        outNumLevels);

  // clean up memory
  for (int level = 0; level < inNumLevels; level++)
    {
      delete inData[level];
    } 

  for (int level = 0; level < outNumLevels; level++)
    {
      delete outData[level];
    } 

  pout() << "Done..." << endl;
  pout() << endl;

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main
