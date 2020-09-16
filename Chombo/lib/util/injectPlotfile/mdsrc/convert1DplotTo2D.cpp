#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "convert1DplotTo2D.H"

#define SAVE_SPACEDIM CH_SPACEDIM
#include "LevelData.H.multidim"
#include "FArrayBox.H.multidim"
#include "Vector.H"
#include "AMRIO.H.multidim"

// At this point, CH_SPACEDIM=0; the .H.multidim files do that.

#include "Slicing.H.transdim"
#include "Injection.H.transdim"

//#define CH_SPACEDIM SAVE_SPACEDIM
//#undef SAVE_SPACEDIM
//#include "UsingNamespace.H"

using namespace Chombo;
using namespace CH_MultiDim;
using std::endl;

int convert1DPlotTo2D(const string& fname1d,
                      const string& fname2d,
                      int num_width)
{
  
  int status =0;

  // first, define 1D dataholders
  D1::Vector<D1::DisjointBoxLayout> vectGrids1d;
  D1::Vector<D1::LevelData<D1::FArrayBox>* >  vectData1d;
  Vector<string> vectNames;
  D1::Box domain1d;
  Real dx;
  Real dt;
  Real time;
  Vector<int> refRatio;
  int numLevels;

  status = D1::ReadAMRHierarchyHDF5(fname1d,
                                    vectGrids1d,
                                    vectData1d,
                                    vectNames,
                                    domain1d,
                                    dx,
                                    dt,
                                    time,
                                    refRatio,
                                    numLevels);


  // now, define 2D dataholders
  D2::Vector<D2::DisjointBoxLayout> vectGrids2d(vectGrids1d.size());
  D2::Vector<D2::LevelData<D2::FArrayBox>* >  vectData2d(vectData1d.size(), NULL);
;
  D2::Box domain2d;
  
  // now do the magic injections from 1d->2d
  D2::SliceSpec slice(1,0);
  injectBox(domain2d, domain1d, slice);
  for (int lev=0; lev<vectGrids1d.size(); lev++)
    {
      injectDisjointBoxLayout(vectGrids2d[lev],
                              vectGrids1d[lev],
                              slice);
      vectData2d[lev] = new D2::LevelData<D2::FArrayBox>;
      injectLevelData(*vectData2d[lev], *vectData1d[lev], slice);
    }

  D2::WriteAMRHierarchyHDF5(fname2d,
                            vectGrids2d,
                            vectData2d,
                            vectNames,
                            domain2d,
                            dx,
                            dt,
                            time,
                            refRatio,
                            numLevels);
  


  return status;
}
