#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// This code takes plotfiles of multilevel, cell-averaged data, 
// one finely resolved, one coarsely resolved, and averages 
// finely-resolved data down to coarsely-resolved data.  Averaging 
// is performed only on the valid regions of each AMR grid.  

// At present, the fine data and coarse data must have equivalent 
// ProblemDomains (i.e. they must occupy the same physical space), 
// and corresponding levels in the different data sets must also 
// occupy the same physical space.

// Perhaps someone else will later enable the code to average 
// down cell-centered, point-valued data.

#include <iostream>
#include <strstream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "CoarseAverage.H"

void init(string&         a_aveDownRoot,
          string&         a_crseRoot,
          string&         a_fineRoot,
          int&            a_aveDownRatio,
          int&            a_intFieldSize,
          int&            a_numCrseStart,
          int&            a_numCrseFinish,
          int&            a_crseStep,
          int&            a_crseMult,
          Vector<string>& a_aveDownVars,
          bool&           a_isTimeDep);

void constructAveDownNames(Vector<string>&       a_aveDownNames,
                           const Vector<string>& a_aveDownVars);

void constructPlotFileName(ostrstream&   a_fileName,
                           const string& a_fileRoot,
                           const int     a_intFieldWidth,
                           const int     a_step);

// Function for MPI
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

  // Set defaults and make declarations
  bool isTimeDep = false;
  
  string aveDownRoot, crseRoot, fineRoot;
  
  Vector<string> aveDownVars;
  
  int aveDownRatio;
  int intFieldSize, numCrseStart, numCrseFinish, crseStep, crseMult;
  
  init(aveDownRoot, crseRoot, fineRoot, aveDownRatio, intFieldSize, 
       numCrseStart, numCrseFinish, crseStep, crseMult, aveDownVars, 
       isTimeDep);
  
  int nStep, fineStep;
  
  if (!isTimeDep)
    {
      crseStep = 1;
      numCrseFinish = numCrseStart;
    }
  
  for (nStep = numCrseStart; nStep <= numCrseFinish; nStep += crseStep)
    {
      fineStep = nStep*crseMult;
      
      ostrstream fineFile;
      ostrstream crseFile;
      ostrstream aveDownFile;
      
      fineFile.fill('0');
      crseFile.fill('0');
      aveDownFile.fill('0');

      if (isTimeDep)
        {
          constructPlotFileName(fineFile, fineRoot, intFieldSize, fineStep);
          constructPlotFileName(crseFile, crseRoot, intFieldSize, nStep);
          constructPlotFileName(aveDownFile, aveDownRoot, intFieldSize, nStep);
        }
      else
        {
          // if not time dependent, file roots are really filenames
          fineFile << fineRoot << ends;
          crseFile << crseRoot << ends;
          aveDownFile << aveDownRoot << ends;
        }
      
      //pout() << "fine Filename = " << fineFile.str() << endl;
      //pout() << "coarse Filename = " << crseFile.str() << endl;
      //pout() << "aveDown Filename = " << aveDownFile.str() << endl;

      // declare memory and read in fine data
      Vector<LevelData<FArrayBox>* > fineSoln;
      Vector<string> fineVars; // fine solution variable names
      Vector<DisjointBoxLayout> fineGrids;
      Box fineDomain;
      Real fineDx, fineDt, fineTime;
      Vector<int> fineRefRatio;
      int fineNumLevels;
      IntVect ghostVect = IntVect::Zero;
      string fineFileName(fineFile.str());
      
      ReadAMRHierarchyHDF5(fineFileName,
                           fineGrids,
                           fineSoln,
                           fineVars,
                           fineDomain,
                           fineDx,
                           fineDt,
                           fineTime,
                           fineRefRatio,
                           fineNumLevels);
      
      // declare memory and read in coarse data
      Vector<LevelData<FArrayBox>* > crseSoln;
      Vector<string> crseVars; // crse soln variable names
      Vector<DisjointBoxLayout> crseGrids;
      Box crseDomain;
      Real crseDx, crseDt, crseTime;
      Vector<int> crseRefRatio;
      int crseNumLevels;
      string crseFileName(crseFile.str());
            
      ReadAMRHierarchyHDF5(crseFileName,
                           crseGrids,
                           crseSoln,
                           crseVars,
                           crseDomain,
                           crseDx,
                           crseDt,
                           crseTime,
                           crseRefRatio,
                           crseNumLevels);
      
      int numFine    = fineVars.size();
      int numCrse = crseVars.size();
      
      int numAveDown = aveDownVars.size();
      // If no aveDownVars were specified
      if (numAveDown == 0)
        {
          // Set aveDownVars to the intersection of fineVars and crseVars
          // This numVars^2 method should be changed to something more efficient
          for (int iFine = 0; iFine < numFine; iFine++)
            {
              for (int iComp = 0; iComp < numCrse; iComp++)
                {
                  if (fineVars[iFine] == crseVars[iComp])
                    {
                      aveDownVars.push_back(fineVars[iFine]);
                      break;
                    }
                }
            }
          
          numAveDown = aveDownVars.size();
        }
      else
        {
          // if aveDownVars were specified, then do a quick check that
          // they're present in both fineVars and crseVars
          for (int errVarNo = 0; errVarNo < aveDownVars.size(); ++errVarNo)
            {
              bool foundCrse = false;
              for (int i = 0; i < numCrse; i++)
                {
                  if (aveDownVars[errVarNo] == crseVars[i])
                    {
                      foundCrse = true;
                    }
                } // end loop over fine variables
              if (!foundCrse)
                {
                  pout() << "aveDownVar " << aveDownVars[errVarNo]
                         << " not found in coarse solution!"
                         << endl;
                  MayDay::Error();
                }
              
              bool foundFine = false;
              for (int i = 0; i < numFine; i++)
                {
                  if (aveDownVars[errVarNo] == fineVars[i])
                    {
                      foundFine = true;
                    }
                } // end loop over fine variables
              if (!foundFine)
                {
                  pout() << "aveDownVar " << aveDownVars[errVarNo]
                         << " not found in fine solution!"
                         << endl;
                  MayDay::Error();
                }
              
            } // end loop over aveDownVars
        } // end if aveDownVars was specified in the inputs file
      
      Vector<string> aveDownNames;
      aveDownNames.resize(numAveDown);
      
      constructAveDownNames(aveDownNames, aveDownVars);
      
      Vector<LevelData<FArrayBox>* > aveDown(crseNumLevels);
      
      // loop over levels
      for (int level = 0; level < crseNumLevels; level++)
        {
          // allocate aveDown -- same domain as coarse solution
          aveDown[level] = new LevelData<FArrayBox>(crseGrids[level],
                                                    numAveDown, 
                                                    ghostVect);
          
          // Construct the CoarseAverage object that will allow us to 
          // average down the fine data
          CoarseAverage crseAve(fineGrids[level],
                                crseGrids[level],
                                numAveDown,
                                aveDownRatio,
                                ghostVect);
          
          // Average down the fine data
          crseAve.averageToCoarse(*(aveDown[level]), *(fineSoln[level]));
          
        } // end loop over levels
      
      // Output the plot file
      WriteAMRHierarchyHDF5(aveDownFile.str(),
                            crseGrids,
                            aveDown,
                            aveDownNames,
                            crseDomain,
                            crseDx,
                            crseDt,
                            crseTime,
                            crseRefRatio,
                            crseNumLevels);
      
      // clean up memory
      for (int level = 0; level < fineNumLevels; level++)
        {
          if (fineSoln[level] != NULL)
            {
              delete fineSoln[level];
              fineSoln[level] = NULL;
            }
        }
      
      for (int level = 0; level < crseNumLevels; level++)
        {
          if (crseSoln[level] != NULL)
            {
              delete crseSoln[level];
              crseSoln[level] = NULL;
            }
          
          if (aveDown[level] != NULL)
            {
              delete aveDown[level];
              aveDown[level] = NULL;
            }
        }
      
    } // end loop over time steps
  
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main


void init(string&         a_aveDownRoot,
          string&         a_crseRoot,
          string&         a_fineRoot,
          int&            a_aveDownRatio,
          int&            a_intFieldSize,
          int&            a_numCrseStart,
          int&            a_numCrseFinish,
          int&            a_crseStep,
          int&            a_crseMult,
          Vector<string>& a_aveDownVars,
          bool&           a_isTimeDep)
{
  ParmParse ppAveDown("avedown");
  
  ppAveDown.get("aveDownRoot", a_aveDownRoot);
  ppAveDown.get("crseRoot", a_crseRoot);
  ppAveDown.get("fineRoot", a_fineRoot);
  ppAveDown.get("aveDownRatio", a_aveDownRatio);
  
  int isTimeDepInt = a_isTimeDep;
  ppAveDown.query("isTimeDep", isTimeDepInt);
  a_isTimeDep = (isTimeDepInt == 1);
  
  a_numCrseStart = 0;
  ppAveDown.query("numCrseStart", a_numCrseStart);
  
  if (a_isTimeDep)
    {
      ppAveDown.get("numCrseFinish", a_numCrseFinish);
      a_crseStep = 1;
      ppAveDown.query("crseStep", a_crseStep);
      ppAveDown.get("mult", a_crseMult);
      a_intFieldSize = 4;
      ppAveDown.query("intFieldSize", a_intFieldSize);
    }
  
  int nVars = ppAveDown.countval("aveDownVars");
  if (nVars > 0)
    {
      a_aveDownVars.resize(nVars);
      ppAveDown.getarr("aveDownVars", a_aveDownVars, 0, nVars);
    }
}


void constructAveDownNames(Vector<string>&       a_aveDownNames,
                           const Vector<string>& a_aveDownVars)
{
  CH_assert(a_aveDownNames.size() == a_aveDownVars.size());

  // for now, don't do anything fancy -- just copy
  for (int i = 0; i < a_aveDownVars.size(); i++)
  {
    a_aveDownNames[i] = a_aveDownVars[i];
  }
}


void constructPlotFileName(ostrstream&   a_fileName,
                           const string& a_fileRoot,
                           const int     a_intFieldWidth,
                           const int     a_step)
{
  // this is kinda klugy, but what are ya gonna do?
  a_fileName << a_fileRoot
             << setw(a_intFieldWidth) << a_step
             << "."
             << setw(1) << CH_SPACEDIM
             << "d.hdf5" << ends;
}
