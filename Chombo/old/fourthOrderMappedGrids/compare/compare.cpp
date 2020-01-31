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
#include <strstream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "LayoutIterator.H"

#include "computeSum.H"
#include "computeNorm.H"

#include "NamespaceHeader.H"

void init(string&         a_exactRoot,
          string&         a_computedRoot,
          Vector<string>& a_errorVars,
          Real&           a_bogusValue,
          bool&           a_verbose)
{
  ParmParse ppCompare("compare");

  int temp = a_verbose;
  ppCompare.query("verbose", temp);
  a_verbose = (temp == 1);

  ppCompare.query("bogus_value", a_bogusValue);

  ppCompare.get("exactRoot", a_exactRoot);
  ppCompare.get("computedRoot", a_computedRoot);

  int nErr = ppCompare.countval("error_var");
  if (nErr > 0)
  {
    a_errorVars.resize(nErr);
    ppCompare.getarr("error_var", a_errorVars, 0, nErr);
  }
}


// This function works on two solutions on equivalent grids.
// It subtracts the computed solution from the exact solution
// on valid cells.  Boxes don't need to be the same.
void computeSameSizeError(Vector<LevelData<FArrayBox>* >&       a_error,
                          const Vector<string>&                 a_errorVars,
                          const Vector<LevelData<FArrayBox>* >& a_computedSoln,
                          const Vector<string>&                 a_computedVars,
                          const Vector<DisjointBoxLayout>&      a_computedGrids,
                          const Real                            a_dx, // coarsest
                          const Vector<int>&                    a_refRatio,
                          const Vector<LevelData<FArrayBox>* >& a_exactSoln,
                          const Vector<string>&                 a_exactVars,
                          Real                                  a_bogus_value) // 0.
{
  int numLevels = a_computedSoln.size();

  CH_assert(a_exactSoln.size() == numLevels);
  CH_assert(a_error.size() == numLevels);
  CH_assert(a_refRatio.size() >= numLevels - 1);

  Real dxLevel = a_dx;

  // outer loop is over levels
  for (int level = 0; level < numLevels; level++)
    {
      LevelData<FArrayBox>& thisLevelError = *a_error[level];
      LevelData<FArrayBox>& thisLevelComputed = *a_computedSoln[level];
      LevelData<FArrayBox>& thisLevelExact = *a_exactSoln[level];

      const DisjointBoxLayout levelGrids = thisLevelComputed.getBoxes();
      const DisjointBoxLayout exactGrids = thisLevelExact.getBoxes();

      // Initialize thisLevelError to bogus value (should be 0).
      DataIterator levelDit = levelGrids.dataIterator();
      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          thisLevelError[levelDit].setVal(a_bogus_value);
        }

      // loop over variables
      for (int nErr = 0; nErr < a_errorVars.size(); nErr++)
        {
          string thisErrVar = a_errorVars[nErr];
          bool done = false;

          // first loop over exact variables and copy into error
          for (int exactComp=0; exactComp<a_exactVars.size(); ++exactComp)
            {
              string thisExactVar = a_exactVars[exactComp];

              // check if this exact variable is "the one"
              if (thisExactVar == thisErrVar)
                {
                  // copy exact solution -> error
                  Interval exactInterval(exactComp, exactComp);
                  Interval errInterval(nErr, nErr);
                  thisLevelExact.copyTo(exactInterval,
                                        thisLevelError,
                                        errInterval);
                  done = true;
                } // end if this exact var is the error var
            } // end loop over exact comps

          if (!done)
            {
              pout() << "Variable " << thisErrVar
                     << " not found in exact solution!!!" << endl;
              MayDay::Error();
            }

          done = false;
          int computedComp = 0;
          // now loop over computed variables and subtract computed solution
          while (!done && (computedComp < a_computedVars.size()))
            {
              if (a_computedVars[computedComp] == thisErrVar)
                {
                  for (levelDit.reset(); levelDit.ok(); ++levelDit)
                    {
                      FArrayBox& thisComputed = thisLevelComputed[levelDit()];
                      FArrayBox& thisError = thisLevelError[levelDit()];

                      thisError.minus(thisComputed, computedComp, nErr, 1);
                    } // end loop over computed/error grids

                  done = true;
                } // if a_computedVar is a_errorVar

              computedComp += 1;
            } // end loop over a_computedVars

          if (!done)
            {
              pout() << "Variable " << thisErrVar  << " not found!!!" << endl;
              MayDay::Error();
            }
        } // end loop over error variables

      // now need to set covered regions to 0
      if (level < numLevels - 1)
        {
          // will need to loop over all boxes in finer level, not just
          // those on this processor...
          const BoxLayout& finerGrids = a_computedSoln[level + 1]->boxLayout();
          LayoutIterator fineLit = finerGrids.layoutIterator();

          // outer loop over this level's grids, since there are fewer of them
          DataIterator levelDit = thisLevelError.dataIterator();
          for (levelDit.reset(); levelDit.ok(); ++levelDit)
            {
              const Box& coarseBox = levelGrids[levelDit()];
              FArrayBox& thisError = thisLevelError[levelDit()];
              int numError = thisError.nComp();

              for (fineLit.reset(); fineLit.ok(); ++fineLit)
                {
                  Box fineBox(finerGrids[fineLit()]);
                  // now coarsen box down to this level
                  fineBox.coarsen(a_refRatio[level]);
                  // if coarsened fine box intersects error's box, set
                  // overlap to 0
                  fineBox &= coarseBox;
                  if (!fineBox.isEmpty())
                    {
                      thisError.setVal(0.0, fineBox, 0, numError);
                    }
                } // end loop over finer-level grids
            } // end loop over this-level grids

          // this is a good place to update dx as well
          dxLevel = dxLevel / a_refRatio[level];
        } // end if there is a finer level

      // finally, if we're not doing ghost cells, do an exchange just
      // to "prettify" the output
      thisLevelError.exchange(thisLevelError.interval());
    } // end loop over levels
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


void processData(Vector<LevelData<FArrayBox>* >& a_out,
                 const Vector<LevelData<FArrayBox>* >& a_data,
                 const Vector<LevelData<FArrayBox>* >& a_J,
                 int a_iComp,
                 bool a_abs,
                 bool a_square,
                 bool a_includeJ)
{
  // Assumes numError == 1 !
  int numLevels = a_data.size();
  for (int level = 0; level < numLevels; level++)
    {
      const LevelData<FArrayBox>& thisLevelData = *a_data[level];
      const LevelData<FArrayBox>& thisLevelJ = *a_J[level];
      const DisjointBoxLayout levelGrids = thisLevelData.disjointBoxLayout();
      a_out[level] = new LevelData<FArrayBox>(levelGrids, 1);
      LevelData<FArrayBox>& outLDF = *a_out[level];
      DataIterator levelDit = levelGrids.dataIterator();
      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          FArrayBox& outFab = outLDF[levelDit];
          const FArrayBox& dataFab = thisLevelData[levelDit];
          outFab.copy(dataFab, levelGrids[levelDit],
                      a_iComp, levelGrids[levelDit], 0, 1);
          if (a_abs) outFab.abs();
          if (a_square) outFab.mult(dataFab, levelGrids[levelDit],
                                    a_iComp, 0, 1);
          if (a_includeJ)
            {
              const FArrayBox& JFab = thisLevelJ[levelDit];
              outFab.mult(JFab, levelGrids[levelDit],
                          0, 0, 1);
            }
        }
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

  // this is the initial value for the error
  Real bogusValue = 0.0;

  string exactRoot, computedRoot;

  Vector<string> errorVars;

  init(exactRoot, computedRoot,
       errorVars, bogusValue, verbose);

  ostrstream exactFile;
  ostrstream computedFile;

  exactFile.fill('0');
  computedFile.fill('0');

  // if not time dependent, file roots are really filenames
  exactFile << exactRoot << ends;
  computedFile << computedRoot << ends;

  //  pout() << "exact Filename = " << exactFile.str() << endl;
  //  pout() << "computed Filename = " << computedFile.str() << endl;

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

  /*
    Find error = computedSoln - exactSoln.
  */
  Vector<LevelData<FArrayBox>* > error(computedNumLevels);
  // allocate error -- same domain as computed solution
  for (int level = 0; level < computedNumLevels; level++)
    {
      error[level] = new LevelData<FArrayBox>(computedGrids[level],
                                              numError);
    }

  // first make sure refRatios are the same
  for (int lev = 0; lev < computedRefRatio.size() - 1; lev++)
    {
      CH_assert(computedRefRatio[lev] == exactRefRatio[lev]);
    }

  CH_assert(exactDx == computedDx);

  if (verbose)
    {
      pout () << "compute sameSize error..." << endl;
    }

  computeSameSizeError(error,
                       errorVars,
                       computedSoln,
                       computedVars,
                       computedGrids,
                       computedDx,
                       computedRefRatio,
                       exactSoln,
                       exactVars,
                       bogusValue);

  if (verbose)
    {
      pout() << "done computing sameSize error" << endl;
    }

  /*
    J is in exact solution.  Alias it.
  */
  bool foundJ = false;
  int iJ;
  for (int iExact = 0; iExact < numExact; iExact++)
    {
      if (exactVars[iExact] == "J")
        {
          iJ = iExact;
          foundJ = true;
        }
    }
  if (!foundJ)
    {
      pout() << "Did not find J in exact solution!" << endl;
      MayDay::Error();
    }
  Interval intvlJ(iJ, iJ);

  Vector<LevelData<FArrayBox>* > J(exactNumLevels);
  for (int level = 0; level < exactNumLevels; level++)
    {
      J[level] = new LevelData<FArrayBox>;
      aliasLevelData(*(J[level]), exactSoln[level], intvlJ);
    }
  // Jerr is the same as J but lives on layout of error.
  Vector<LevelData<FArrayBox>* > Jerr(exactNumLevels);
  for (int level = 0; level < exactNumLevels; level++)
    {
      Jerr[level] = new LevelData<FArrayBox>(computedGrids[level], 1);
      (*J[level]).copyTo(*Jerr[level]);
    }

  /*
    Get component indices of errorVars in exactVars and computedVars.
   */
  // loop over error variables:  assume there is only one.
  int iCompInExact;
  int iCompInComputed;
  CH_assert(numError == 1);
  for (int nErr = 0; nErr < numError; nErr++)
    {
      string thisErrVar = errorVars[nErr];
      // Look for thisErrVar in exact variables.
      for (int exactComp = 0; exactComp < exactVars.size(); exactComp++)
        {
          string thisExactVar = exactVars[exactComp];
          // check if this exact variable is "the one"
          if (thisExactVar == thisErrVar)
            {
              iCompInExact = exactComp;
            }
        }
      // Look for thisErrVar in exact variables.
      for (int computedComp = 0; computedComp < numComputed; computedComp++)
        {
          string thisComputedVar = computedVars[computedComp];
          // check if this computed variable is "the one"
          if (thisComputedVar == thisErrVar)
            {
              iCompInComputed = computedComp;
            }
        }
    }
  Interval intvlCompInExact(iCompInExact, iCompInExact);
  Interval intvlCompInComputed(iCompInComputed, iCompInComputed);
  Interval intvlCompInError(0, 0);

  Vector<LevelData<FArrayBox>* > exactSolnAbs(exactNumLevels);
  processData(exactSolnAbs, exactSoln, J, iCompInExact, true, false, false);
  Vector<LevelData<FArrayBox>* > errorAbs(exactNumLevels);
  processData(errorAbs, error, J, 0, true, false, false);

  // We'll take integrals of these.
  Vector<LevelData<FArrayBox>* > exactSolnAbsJ(exactNumLevels);
  processData(exactSolnAbsJ, exactSoln, J, iCompInExact, true, false, true);
  Vector<LevelData<FArrayBox>* > errorAbsJ(exactNumLevels);
  processData(errorAbsJ, error, Jerr, 0, true, false, true);

  Vector<LevelData<FArrayBox>* > exactSoln2J(exactNumLevels);
  processData(exactSoln2J, exactSoln, J, iCompInExact, false, true, true);
  Vector<LevelData<FArrayBox>* > error2J(exactNumLevels);
  processData(error2J, error, Jerr, 0, false, true, true);

  /*
    Find max and min of exact and computed solution of error variables.
   */

  Vector<Real> exactMin(numError);
  Vector<Real> exactMax(numError);
  Vector<Real> exactAbsMax(numError);
  Vector<Real> computedMin(numError);
  Vector<Real> computedMax(numError);
  Vector<Real> errorMin(numError);
  Vector<Real> errorMax(numError);
  Vector<Real> errorAbsMax(numError);
  int lbase = 0;
  Interval intvl0(0, 0);
  // loop over error variables
  for (int nErr = 0; nErr < numError; nErr++)
    {
      errorMin[nErr] =
        computeMin(error, computedRefRatio, intvl0, lbase);
      errorMax[nErr] =
        computeMax(error, computedRefRatio, intvl0, lbase);
      errorMax[nErr] =
        computeMax(error, computedRefRatio, intvl0, lbase);
      errorAbsMax[nErr] =
        computeMax(errorAbs, computedRefRatio, intvl0, lbase);
      exactAbsMax[nErr] =
        computeMax(exactSolnAbs, exactRefRatio, intvl0, lbase);
      string thisErrVar = errorVars[nErr];
      // Look for thisErrVar in exact variables.
      for (int exactComp = 0; exactComp < exactVars.size(); exactComp++)
        {
          string thisExactVar = exactVars[exactComp];
          // check if this exact variable is "the one"
          if (thisExactVar == thisErrVar)
            {
              Interval intvlComp(exactComp, exactComp);
              exactMin[nErr] =
                computeMin(exactSoln, exactRefRatio, intvlComp, lbase);
              exactMax[nErr] =
                computeMax(exactSoln, exactRefRatio, intvlComp, lbase);
            }
        }
      // Look for thisErrVar in computed variables.
      for (int computedComp = 0; computedComp < numComputed; computedComp++)
        {
          string thisComputedVar = computedVars[computedComp];
          // check if this computed variable is "the one"
          if (thisComputedVar == thisErrVar)
            {
              Interval intvlComp(computedComp, computedComp);
              computedMin[nErr] =
                computeMin(computedSoln, computedRefRatio, intvlComp, lbase);
              computedMax[nErr] =
                computeMax(computedSoln, computedRefRatio, intvlComp, lbase);
            }
        }
    }

  // now compute norms

  // If you really want the integral, multiply computeSum by scaling.
  // But we take ratios of integrals, so it doesn't matter.
  Real scaling = 1. / (4. * M_PI);

  Real exactSolnAbsIntegral =
    computeSum(exactSolnAbsJ, exactRefRatio, exactDx, intvl0, lbase);

  Real exactSoln2Integral =
    computeSum(exactSoln2J, exactRefRatio, exactDx, intvl0, lbase);

  Real errorAbsIntegral =
    computeSum(errorAbsJ, exactRefRatio, exactDx, intvl0, lbase);

  Real error2Integral =
    computeSum(error2J, exactRefRatio, exactDx, intvl0, lbase);

  Real L1 = errorAbsIntegral / exactSolnAbsIntegral;
  Real L2 = sqrt(error2Integral / exactSoln2Integral);
  Real Linf = errorMax[0] / exactAbsMax[0];
  Real deltaPhi0 = exactMax[0] - exactMin[0];
  Real phiMax = (computedMax[0] - exactMax[0]) / deltaPhi0;
  Real phiMin = (computedMin[0] - exactMin[0]) / deltaPhi0;

  pout() << L1 << " "
         << L2 << " "
         << Linf << " "
         << phiMax << " "
         << phiMin << endl;

  Real exactSolnIntegral = scaling *
    computeSum(exactSoln, exactRefRatio, exactDx, intvlCompInExact, lbase);
  Real computedSolnIntegral = scaling *
    computeSum(computedSoln, exactRefRatio, exactDx, intvlCompInComputed, lbase);
  Real diffSolnIntegral = computedSolnIntegral - exactSolnIntegral;
  if (false)
    {
      pout() << "sums "
             << exactSolnIntegral << " "
             << computedSolnIntegral << " "
             << diffSolnIntegral << endl;
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
      if (error[level] != NULL)
        {
          delete error[level];
          error[level] = NULL;
          delete exactSolnAbs[level];
          delete errorAbs[level];
          delete exactSolnAbsJ[level];
          delete errorAbsJ[level];
          delete exactSoln2J[level];
          delete error2J[level];
          delete J[level];
          delete Jerr[level];
        }
    }

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main


#include "NamespaceFooter.H"
