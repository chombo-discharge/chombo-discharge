#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "localfunctions.H"
#include "parstream.H"
//#include "tdglF_F.H"

//#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "BRMeshRefine.H"
#include "FineInterp.H"


int outputHDF5(const Vector<LevelData<FArrayBox>* >& vectData,
               const Vector<DisjointBoxLayout>& vectGrids,
               const Vector<Box>& vectDomain,
               const Vector<int>& vectRatio,
               const int numlevels,
               const string& filename)
{
  string dataName("blue_agave_level");
  Vector<string> vectName(1);
  vectName[0] = dataName;

  Box domain = vectDomain[0];
  Real dx = 1.;
  Real dt = 1.;
  Real time = 1.;

#ifdef CH_USE_HDF5
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectData,
                        vectName,
                        domain,
                        dx, dt, time,
                        vectRatio,
                        numlevels);
#endif

  return 0;
}


int outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
               const Vector<LevelData<FArrayBox>* >& vectRhs,
               const Vector<DisjointBoxLayout>& vectGrids,
               const Vector<Box>& vectDomain,
               const Vector<int>& vectRatio,
               const int numlevels,
               const bool verbose)
{
  string phiName("phi");
  string rhsName("rhs");
  Vector<string> vectName(2);
  vectName[0] = phiName;
  vectName[1] = rhsName;
  Box domain = vectDomain[0];
  Real dx = 1.;
  Real dt = 1.;
  Real time = 1.;
  Vector<LevelData<FArrayBox>* > vectPhiAndRHS(numlevels, NULL);
  for (int ilev = 0; ilev < numlevels; ilev++)
  {
    vectPhiAndRHS[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], 2);
    Interval phiInterval(0,0);
    Interval rhsInterval(1,1);
    vectPhi[ilev]->copyTo(vectPhi[ilev]->interval(),
                          *vectPhiAndRHS[ilev],
                          phiInterval);
    vectRhs[ilev]->copyTo(vectRhs[ilev]->interval(),
                          *vectPhiAndRHS[ilev],
                          rhsInterval);
  }
#ifdef CH_USE_HDF5
  string filename("poissonOut.hdf5");
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhiAndRHS,
                        vectName,
                        domain,
                        dx, dt, time,
                        vectRatio,
                        numlevels);
#endif

  for (int ilev = 0; ilev < numlevels; ilev++)
  {
    delete vectPhiAndRHS[ilev];
  }

  return 0;
}

