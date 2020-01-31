#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// this code is the chombo version of the diagnostic utility
// which takes two plotfiles (one "exact" solution, and
// one "computed" solution) and writes them out, line by line.

#include <iostream>
#include <strstream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "LevelData.H"
#include "DataIterator.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

void init(string&         a_exactRoot,
          string&         a_computedRoot,
          string&         a_outputRoot,
          Vector<string>& a_errorVars,
          bool&           a_verbose)
{
  ParmParse ppDiagnostic("diagnostic");

  int temp = a_verbose;
  ppDiagnostic.query("verbose", temp);
  a_verbose = (temp == 1);

  ppDiagnostic.get("exact", a_exactRoot);
  ppDiagnostic.get("computed", a_computedRoot);
  ppDiagnostic.get("output", a_outputRoot);

  int nErr = ppDiagnostic.countval("var");
  if (nErr > 0)
  {
    a_errorVars.resize(nErr);
    ppDiagnostic.getarr("var", a_errorVars, 0, nErr);
  }
}


int findComp(const Vector<string>& a_listVars,
             const string& a_thisVar)
{
  bool found = false;
  int retval = -1;
  int listLen = a_listVars.size();
  for (int comp = 0; comp < listLen; ++comp)
    {
      const string& listVarHere = a_listVars[comp];
      if (listVarHere == a_thisVar)
        {
          found = true;
          retval = comp;
        }
    }
  if (!found)
    {
      pout() << "Variable " << a_thisVar
             << " not found in solution!" << endl;
      MayDay::Error();
    }
  return retval;
}

// This function works on two solutions on equivalent grids.
// It subtracts the computed solution from the exact solution
// on valid cells.  Boxes don't need to be the same.
void writeBoth(
               const string&                         a_outputFileRoot,
               const Vector<string>&                 a_errorVars,
               const Vector<LevelData<FArrayBox>* >& a_computedSoln,
               const Vector<string>&                 a_computedVars,
               const Vector<DisjointBoxLayout>&      a_computedGrids,
               const Real                            a_dx, // coarsest
               const Vector<int>&                    a_refRatio,
               const Vector<LevelData<FArrayBox>* >& a_exactSoln,
               const Vector<string>&                 a_exactVars)
{
  int numLevels = a_computedSoln.size();

  CH_assert(a_exactSoln.size() == numLevels);
  CH_assert(a_refRatio.size() >= numLevels - 1);

  Interval zeroInterval(0, 0);
  // loop over variables
  for (int nErr = 0; nErr < a_errorVars.size(); nErr++)
    {
      string thisErrVar = a_errorVars[nErr];

      int exactComp = findComp(a_exactVars, thisErrVar);
      int computedComp = findComp(a_computedVars, thisErrVar);
      Interval exactInterval(exactComp, exactComp);

      ostrstream outputFileStream;
      outputFileStream << a_outputFileRoot
                       << "." << thisErrVar
                       << ".txt"
                       << ends;
      string outputFileName(outputFileStream.str());
      FILE *outputFile;
      outputFile = fopen(outputFileName.c_str(), "w");
      if (outputFile == NULL)
        {
          pout() << "Could not open " << outputFileName
                 << " for write" << endl;
          MayDay::Error();
        }

      // Real dxLevel = a_dx;

      for (int level = 0; level < numLevels; level++)
        {
          LevelData<FArrayBox>& thisLevelComputed = *a_computedSoln[level];
          LevelData<FArrayBox>& thisLevelExact = *a_exactSoln[level];

          const DisjointBoxLayout& levelGrids = thisLevelComputed.getBoxes();
          // const DisjointBoxLayout& exactGrids = thisLevelExact.getBoxes();

          LevelData<FArrayBox> thisLevelExactCopy(levelGrids, 1);

          DataIterator levelDit = levelGrids.dataIterator();
          // copy exact solution -> error
          thisLevelExact.copyTo(exactInterval,
                                thisLevelExactCopy,
                                zeroInterval);
          for (levelDit.reset(); levelDit.ok(); ++levelDit)
            {
              const Box& baseBox = levelGrids[levelDit];
              const FArrayBox& thisComputed = thisLevelComputed[levelDit];
              const FArrayBox& thisExactCopy = thisLevelExactCopy[levelDit];
              BoxIterator bit(baseBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  Real cellComputed = thisComputed(iv, computedComp);
                  Real cellExact = thisExactCopy(iv, 0);
                  fprintf(outputFile, "%.11e %.11e\n",
                          cellComputed, cellExact);
                }
            } // end loop over computed/error grids
        } // end loop over levels
      fclose(outputFile);
    } // end loop over error variables
}


void constructErrorNames(Vector<string>&       a_errorNames,
                         const Vector<string>& a_errorVars)
{
  CH_assert(a_errorNames.size() == a_errorVars.size());

  // for now, don't do anything fancy -- just copy
  for (int i = 0; i < a_errorVars.size(); i++)
  {
    a_errorNames[i] = a_errorVars[i];
  }
}


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

  pout().setf(ios::scientific);
  pout().precision(4);

  // set defaults
  bool verbose = false;

  string exactRoot, computedRoot, outputRoot;

  Vector<string> errorVars;

  init(exactRoot, computedRoot, outputRoot, errorVars, verbose);

  ostrstream exactFile;
  ostrstream computedFile;

  exactFile.fill('0');
  computedFile.fill('0');

  // if not time dependent, file roots are really filenames
  exactFile << exactRoot << ends;
  computedFile << computedRoot << ends;

  //  pout() << "exact Filename = " << exactFile.str() << endl;
  //  pout() << "computed Filename = " << computedFile.str() << endl;

  //  string outputFileName(outputFile.str());

  /*
    Load exact solution.
  */
  // declare memory
  Vector<LevelData<FArrayBox>* > exactSoln;
  Vector<string> exactVars; // exact solution variable names
  Vector<DisjointBoxLayout> exactGrids;
  Box exactDomain;
  Real exactDx, exactDt, exactTime;
  Vector<int> exactRefRatio;
  int exactNumLevels;
  string exactFileName(exactFile.str());

  // get exact solution
  if (verbose)
    {
      pout() << "read exact solution..." << endl;
    }

  ReadAMRHierarchyHDF5(exactFileName,
                       exactGrids,
                       exactSoln,
                       exactVars,
                       exactDomain,
                       exactDx,
                       exactDt,
                       exactTime,
                       exactRefRatio,
                       exactNumLevels);

  if (verbose)
    {
      pout () << "done reading exact soln" << endl;
    }

  /*
    Load computed solution.
  */
  Vector<LevelData<FArrayBox>* > computedSoln;
  Vector<string> computedVars; // computed soln variable names
  Vector<DisjointBoxLayout> computedGrids;
  Box computedDomain;
  Real computedDx, computedDt, computedTime;
  Vector<int> computedRefRatio;
  int computedNumLevels;
  string computedFileName(computedFile.str());

  // now read in computed solution
  if (verbose)
    {
      pout() << "read computed solution..." << endl;
    }

  ReadAMRHierarchyHDF5(computedFileName,
                       computedGrids,
                       computedSoln,
                       computedVars,
                       computedDomain,
                       computedDx,
                       computedDt,
                       computedTime,
                       computedRefRatio,
                       computedNumLevels);

  if (verbose)
    {
      pout() << "done reading computed solution" << endl;
    }

  // reality check
  if (computedDomain != exactDomain)
    {
      MayDay::Error("Incompatible exact and computed domains for sameSize comparison");
    }

  int numExact    = exactVars.size();
  int numComputed = computedVars.size();

  int numError = errorVars.size();
  // If no errorVars were specified
  if (numError == 0)
    {
      // Set errorVars to the intersection of exactVars and computedVars
      // This numVars^2 method should be changed to something more efficient
      for (int iExact = 0; iExact < numExact; iExact++)
        {
          for (int iComp = 0; iComp < numComputed; iComp++)
            {
              if (exactVars[iExact] == computedVars[iComp])
                {
                  errorVars.push_back(exactVars[iExact]);
                  break;
                }
            }
        }

      numError = errorVars.size();
    }
  else
    {
      // if errorVars were specified, then do a quick check that
      // they're present in both exactVars and computedVars
      for (int errVarNo = 0; errVarNo<errorVars.size(); ++errVarNo)
        {
          bool foundComputed = false;
          for (int i=0; i<numComputed; i++)
            {
              if (errorVars[errVarNo] == computedVars[i])
                {
                  foundComputed = true;
                }
            } // end loop over exact variables
          if (!foundComputed)
            {
              pout() << "errorVar " << errorVars[errVarNo]
                     << " not found in computed solution!"
                     << endl;
              MayDay::Error();
            }

          bool foundExact = false;
          for (int i=0; i<numExact; i++)
            {
              if (errorVars[errVarNo] == exactVars[i])
                {
                  foundExact = true;
                }
            } // end loop over exact variables
          if (!foundExact)
            {
              pout() << "errorVar " << errorVars[errVarNo]
                     << " not found in exact solution!"
                     << endl;
              MayDay::Error();
            }

        } // end loop over errorVars
    } // end if errorVars was specified in the inputs file


  Vector<string> errorNames(numError);
  constructErrorNames(errorNames, errorVars);

  // first make sure refRatios are the same
  for (int lev = 0; lev < computedRefRatio.size() - 1; lev++)
    {
      CH_assert(computedRefRatio[lev] == exactRefRatio[lev]);
    }

  CH_assert(exactDx == computedDx);

  if (verbose)
    {
      pout () << "write both..." << endl;
    }

  writeBoth(outputRoot,
            errorVars,
            computedSoln,
            computedVars,
            computedGrids,
            computedDx,
            computedRefRatio,
            exactSoln,
            exactVars);

  if (verbose)
    {
      pout() << "done writing both" << endl;
    }

  // clean up memory
  for (int level = 0; level < exactNumLevels; level++)
    {
      if (exactSoln[level] != NULL)
        {
          delete exactSoln[level];
          exactSoln[level] = NULL;
        }
    }

  for (int level = 0; level < computedNumLevels; level++)
    {
      if (computedSoln[level] != NULL)
        {
          delete computedSoln[level];
          computedSoln[level] = NULL;
        }
    }

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main


#include "NamespaceFooter.H"
