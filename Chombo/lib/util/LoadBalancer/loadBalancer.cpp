#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// this code is the chombo version of the compare utility
// which takes two plotfiles (one fine "exact" solution, and
// one coarser "computed" solution) and computes L1, L2, and
// Max "errors".  Norms are computed only on valid regions
// of each AMR grid.  assumes that exact solution is a single
// fine grid.

#include <iostream>
#include <fstream>
#include <strstream>
#include <string>
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::string;
#include <cstdio>

#include "AMRIO.H"
#include "ParmParse.H"
#include "DisjointBoxLayout.H"
#include "LoadBalance.H"

void readBoxes(string& a_fileName,
               Vector<Vector<Box> >& a_vectBoxes,
               int& a_numLevels, int a_verbosity);

void writeBoxes(string& a_filename,
                Vector<Vector<Box> >& a_vectBoxes,
                Vector<Vector<int> >& a_vectProcs,
                int a_numLevels, int a_verbosity);

// One more function for MPI
void dumpmemoryatexit();

int
main(int argc, char* argv[])
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
  if (argc< 2)
  {
    cerr << "  need inputs file" << endl;
    abort();
  }

  char* in_file = argv[1];
  ParmParse pp(argc-2, argv+2,NULL,in_file);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  string inFileName, outFileName;
  Vector<Vector<Box> > vectBoxes;
  Vector<Vector<int> > vectProcAssign;
  int numProcs, numLevels;
  int verbosity = 0;

  ParmParse localPP("main");

  localPP.get("inFile", inFileName);
  localPP.get("outFile", outFileName);
  localPP.get("numProcs", numProcs);
  localPP.query("verbosity", verbosity);

  if (verbosity > 0) pout() << "reading Grids file" << endl;
  readBoxes(inFileName, vectBoxes, numLevels, verbosity);
  vectProcAssign.resize(numLevels);

  if (verbosity > 0) pout() << "Calling Load Balance..." << endl;
  // now loop over levels and call LoadBalance function
  for (int lev=0; lev<numLevels; lev++)
    {
      if (verbosity > 1) pout() << "... Level " << lev << endl;
      LoadBalance(vectProcAssign[lev], vectBoxes[lev], numProcs);
    }

  if (verbosity > 0) pout() << "Writing new grids file" << endl;
  writeBoxes(outFileName, vectBoxes, vectProcAssign, numLevels, verbosity);

  if (verbosity > 0) pout() << "exiting..." << endl;

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main

void readBoxes(string& a_fileName,
               Vector<Vector<Box> >& a_vectBoxes,
               int& a_numLevels, int verbosity)
{
  // this is the same as what's in the AMRPoisson example code

  if (verbosity > 2) pout() << "opening " << a_fileName << endl;
  // read grids (and proc assignments, if MPI) from file
  ifstream is(a_fileName.c_str(), ios::in);
  if (is.fail())
    {
      pout() << "Cannot open grids file " << a_fileName << endl;
    }

  // format of file -- number of levels, then for each level starting with
  // level 1, number of grids on level, list of boxes (along with proc
  // assignments if mpi -- since we're going to be calling load balance
  // anyway, we just throw these away.)
  int in_numLevels;
  while (is.get() != '=');
  is >> in_numLevels;
  if (verbosity >3)
    {
      pout() << "NumLevels in grids file = " << in_numLevels << endl;
    }

  a_numLevels = in_numLevels;

  while (is.get() != '\n');
  a_vectBoxes.resize(a_numLevels);
  int ngrid;
  for (int lev=0; lev<a_numLevels; lev++)
    {
      while (is.get() != ':');
      is >> ngrid;
      if (verbosity >3)
        {
          pout() << "level " << lev << " numGrids = " << ngrid << endl;
          pout() << "grids = " ;
        }
      while (is.get() != '\n');
      a_vectBoxes[lev].resize(ngrid);
      Box bx;

      for (int i=0; i<ngrid; i++)
        {
          is >> bx;

          // while there may be processor assigments here
          // in the file, there's no point in even reading them...

          //while (is.get() != '[');
          //is >> this_proc;

          // move on to the next box
          while (char ch = is.get())
            {
              if (ch == '#') break;
              if (ch == '\n') break;
            }

          if (verbosity >4)
          {
            pout() << bx << endl;
          }
          a_vectBoxes[lev][i] = bx;
        } // end loop over boxes on this level
    } // end loop over levels

  if (verbosity > 3) pout() << "done reading file" << endl;
}

void writeBoxes(string& a_filename,
                Vector<Vector<Box> >& a_vectBoxes,
                Vector<Vector<int> >& a_vectProcs,
                int a_numLevels, int verbosity)
{

  if (verbosity > 2) pout() << "opening " << a_filename << endl;

  ofstream os(a_filename.c_str(), ios::out);
  if (os.fail())
    {
      pout() << "cannot open grid output file " << os << endl;
      MayDay::Error();
    }

  os << "NumLevels = " << a_numLevels << endl;
  for (int lev=0; lev < a_numLevels; lev++)
    {
      if (verbosity > 3) pout() << "writing level " << lev << endl;

      // this is a really silly way to do this, but it's
      // definitely the simplest way to get the formatting
      // correct...
      DisjointBoxLayout thisLevelGrids(a_vectBoxes[lev],
                                       a_vectProcs[lev]);

      os << "Number of Grids on Level " << lev << ": "
         << thisLevelGrids.size() << endl;
      // note that we don't need an endl here because
      // BoxLayout::operator<< provides a '\n' at the end
      os << thisLevelGrids;
    }
  os.close();

  if (verbosity > 3) pout() << "done writing file" << endl;

}
