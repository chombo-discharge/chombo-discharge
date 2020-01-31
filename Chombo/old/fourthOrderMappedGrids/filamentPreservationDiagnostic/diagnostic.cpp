#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// This code returns the filament preservation diagnostic.
// For any real tau,
// A(tau, t) = sum_{cell i: scalar[i] > tau} volume[i] .

#include <iostream>
#include <strstream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "LevelData.H"
#include "DataIterator.H"
#include "BoxIterator.H"

#include "computeSum.H"
#include "computeNorm.H"

#include "NamespaceHeader.H"

void init(
          string&         a_fileName,
          Vector<string>& a_errorVars,
          Real&           a_bogusValue,
          bool&           a_verbose,
          Real&           a_tauMin,
          Real&           a_tauInterval,
          int&            a_tauNumber)
{
  ParmParse ppDiagnostic("diagnostic");

  ppDiagnostic.get("file", a_fileName);

  int verboseInt = a_verbose;
  ppDiagnostic.query("verbose", verboseInt);
  a_verbose = (verboseInt == 1);

  ppDiagnostic.query("bogus_value", a_bogusValue);

  int nErr = ppDiagnostic.countval("var");
  if (nErr > 0)
    {
      a_errorVars.resize(nErr);
      ppDiagnostic.getarr("var", a_errorVars, 0, nErr);
    }

  a_tauMin = 0.1;
  ppDiagnostic.query("tau_min", a_tauMin);
  a_tauInterval = 0.05;
  ppDiagnostic.query("tau_interval", a_tauInterval);
  a_tauNumber = 19;
  ppDiagnostic.query("tau_number", a_tauNumber);
}


Real computeDiagnostic(
                       const Vector<LevelData<FArrayBox>* >& a_mySoln,
                       const Vector<LevelData<FArrayBox>* >& a_volume,
                       const Vector<DisjointBoxLayout>&      a_myGrids,
                       Real                                  a_dx, // coarsest
                       const Vector<int>&                    a_refRatio,
                       Real                                  a_tau,
                       Real                                  a_bogus_value) // 0.
{
  int numLevels = a_mySoln.size();

  CH_assert(a_refRatio.size() >= numLevels - 1);

  Real dxLevel = a_dx;

  Real retval = a_bogus_value; // should be 0.
  // outer loop is over levels
  for (int level = 0; level < numLevels; level++)
    {
      LevelData<FArrayBox>& myLevel = *a_mySoln[level];
      LevelData<FArrayBox>& volumeLevel = *a_volume[level];

      const DisjointBoxLayout& levelGrids = myLevel.getBoxes();

      DataIterator levelDit = levelGrids.dataIterator();

      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          const FArrayBox& myFab = myLevel[levelDit];
          const FArrayBox& volumeFab = volumeLevel[levelDit];
          const Box& baseBox = levelGrids[levelDit];
          BoxIterator bit(baseBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              if (myFab(iv, 0) > a_tau)
                {
                  retval += volumeFab(iv, 0);
                }
            }
        }
      // this is a good place to update dx as well
      dxLevel = dxLevel / a_refRatio[level];
    } // end loop over levels

  return retval;
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
  // parse command line
  // ------------------------------------------------
  if (argc < 2)
  {
    cerr << "Sample usage:" << endl;
    cerr << "./diagnostic...ex" << endl;
    cerr << "diagnostic.file=file.hdf5 diagnostic.var=scalar" << endl;
    cerr << "diagnostic.tau_min=0.1 diagnostic.tau_interval=0.05 diagnostic.tau_number=19" << endl;
    abort();
  }

  // char* in_file = argv[1];
  // ParmParse pp(argc-2, argv+2, NULL, in_file);
  ParmParse pp(argc-1, argv+1, NULL);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  pout().setf(ios::scientific);
  pout().precision(4);

  // set defaults
  bool verbose = false;

  // this is the initial value for the error
  Real bogusValue = 0.0;

  string fileName;

  Vector<string> errorVars;

  Real tauMin, tauInterval;
  int tauNumber;

  init(fileName, errorVars, bogusValue, verbose,
       tauMin, tauInterval, tauNumber);

  ostrstream myFile;

  myFile.fill('0');

  // if not time dependent, file roots are really filenames
  myFile << fileName << ends;

  //  pout() << "filename = " << myFile.str() << endl;

  /*
    Load solution.
  */
  // declare memory
  Vector<LevelData<FArrayBox>* > mySoln;
  Vector<string> myVars; // solution variable names
  Vector<DisjointBoxLayout> myGrids;
  Box myDomain;
  Real myDx, myDt, myTime;
  Vector<int> myRefRatio;
  int myNumLevels;
  string myFileName(myFile.str());

  // get solution
  if (verbose)
    {
      pout() << "read solution..." << endl;
    }

  ReadAMRHierarchyHDF5(myFileName,
                       myGrids,
                       mySoln,
                       myVars,
                       myDomain,
                       myDx,
                       myDt,
                       myTime,
                       myRefRatio,
                       myNumLevels);

  if (verbose)
    {
      pout () << "done reading soln" << endl;
    }

  int numVars    = myVars.size();

  int numError = errorVars.size();
  // If no errorVars were specified
  if (numError == 0)
    {
      errorVars = myVars;
      numError = errorVars.size();
    }
  else
    {
      // If errorVars were specified, then do a quick check that
      // they're present in myVars.
      for (int errVarNo = 0; errVarNo<errorVars.size(); ++errVarNo)
        {
          bool found = false;
          for (int i=0; i<numVars; i++)
            {
              if (errorVars[errVarNo] == myVars[i])
                {
                  found = true;
                }
            } // end loop over exact variables
          if (!found)
            {
              pout() << "errorVar " << errorVars[errVarNo]
                     << " not found in input file!"
                     << endl;
              MayDay::Error();
            }
        } // end loop over errorVars
    } // end if errorVars was specified in the inputs file

  /*
    volume is in input file.  Alias it.
  */
  bool foundVolume = false;
  int iVolume;
  for (int icomp = 0; icomp < numVars; icomp++)
    {
      if (myVars[icomp] == "volume")
        {
          iVolume = icomp;
          foundVolume = true;
        }
    }
  if (!foundVolume)
    {
      pout() << "Did not find volume in input file!" << endl;
      MayDay::Error();
    }
  Interval intvlVolume(iVolume, iVolume);

  Vector<LevelData<FArrayBox>* > volume(myNumLevels);
  for (int level = 0; level < myNumLevels; level++)
    {
      volume[level] = new LevelData<FArrayBox>;
      aliasLevelData(*(volume[level]), mySoln[level], intvlVolume);
    }

  /*
    Get adjustedSoln, which is mySoln with zeroes on covered cells.
  */
  Vector<LevelData<FArrayBox>* > adjustedSoln(myNumLevels);
  for (int level = 0; level < myNumLevels; level++)
    {
      const LevelData<FArrayBox>& myLevel = *mySoln[level];
      int ncomp = myLevel.nComp();
      const DisjointBoxLayout& levelGrids = myGrids[level];
      adjustedSoln[level] = new LevelData<FArrayBox>(levelGrids, ncomp);
      LevelData<FArrayBox>& adjustedLevel = *adjustedSoln[level];
      myLevel.copyTo(adjustedLevel);
      if (level < myNumLevels-1)
        { // Copy zero on coarsened myGrids[level+1] to myGrids[level].
          DisjointBoxLayout& finerGrids = myGrids[level+1];
          DisjointBoxLayout finerGridsCoarsened(myGrids[level+1]);
          finerGridsCoarsened.coarsen(myRefRatio[level]);
          LevelData<FArrayBox> zeroData(finerGridsCoarsened, ncomp);
          DataIterator finerDit = finerGrids.dataIterator();
          for (finerDit.begin(); finerDit.ok(); ++finerDit)
            {
              zeroData[finerDit].setVal(0.);
            }
          zeroData.copyTo(adjustedLevel);
        }
    }

  /*
    Now compute the diagnostic.
  */

  if (verbose)
    {
      pout () << "compute diagnostic ..." << endl;
    }

  pout() << "tau";
  for (int iTau = 0; iTau < tauNumber; iTau++)
    {
      Real tau = tauMin + iTau * tauInterval;
      pout() << " " << tau;
    }
  pout() << endl;

  // loop over variables
  for (int nErr = 0; nErr < errorVars.size(); nErr++)
    {
      string thisErrVar = errorVars[nErr];
      bool done = false;
      for (int myComp=0; myComp < myVars.size(); ++myComp)
        {
          string thisMyVar = myVars[myComp];

          // check if this exact variable is "the one"
          if (thisMyVar == thisErrVar)
            {
              // myComp in myVars
              // nErr in errorVars
              pout() << thisMyVar;
              Interval myInterval(myComp, myComp);
              Vector<LevelData<FArrayBox>* > myComponent(myNumLevels);
              for (int level = 0; level < myNumLevels; level++)
                {
                  myComponent[level] = new LevelData<FArrayBox>;
                  aliasLevelData(*(myComponent[level]),
                                 adjustedSoln[level],
                                 myInterval);
                }

              for (int iTau = 0; iTau < tauNumber; iTau++)
                {
                  Real tau = tauMin + iTau * tauInterval;
                  Real diag = computeDiagnostic(myComponent, // comp 0 only
                                                volume,
                                                myGrids,
                                                myDx,
                                                myRefRatio,
                                                tau,
                                                bogusValue);
                  pout() << " " << diag;
                }
              pout() << endl;
              done = true;
              for (int level = 0; level < myNumLevels; level++)
                {
                  delete myComponent[level];
                  myComponent[level] = NULL;
                }
            } // end if this my var is the error var
        } // end loop over my comps
      if (!done)
        {
          pout() << "Variable " << thisErrVar
                 << " not found in input file!" << endl;
          MayDay::Error();
        }
    }

  if (verbose)
    {
      pout() << "done computing diagnostic" << endl;
    }

  // clean up memory
  for (int level = 0; level < myNumLevels; level++)
    {
      if (mySoln[level] != NULL)
        {
          delete mySoln[level];
          delete adjustedSoln[level];
          delete volume[level];
          mySoln[level] = NULL;
          adjustedSoln[level] = NULL;
          volume[level] = NULL;
        }
    }

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main


#include "NamespaceFooter.H"
