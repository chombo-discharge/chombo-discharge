#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "slice4Dplotfile.H"
#include "Vector.H"


#define SAVE_SPACEDIM CH_SPACEDIM
#include "LevelData.H.multidim"
#include "FArrayBox.H.multidim"
#include "AMRIO.H.multidim"

// At this point, CH_SPACEDIM=0; the .H.multidim files do that.

#include "Slicing.H.transdim"
#include "Injection.H.transdim"
#include "RealVect.H.multidim"

//#define CH_SPACEDIM SAVE_SPACEDIM
//#undef SAVE_SPACEDIM
//#include "UsingNamespace.H"

using namespace Chombo;
using namespace CH_MultiDim;
using std::endl;

int slice4Dplotfile(const string& fname4d,
                    const string& slicefname,
                    const Vector<int>& slice_dir,
                    const Vector<int>& slice_comp)
{
  int status =0;

  // first, define 4D dataholders
  D4::Vector<D4::DisjointBoxLayout> vectGrids4d;
  D4::Vector<D4::LevelData<D4::FArrayBox>* >  vectData4d;
  Vector<string> vectNames;
  D4::Box domain4d;
  D4::RealVect dx4d;
  Real dt;
  Real time;
  Vector<D4::IntVect> refRatio4d;
  int numLevels;

  status = D4::ReadAnisotropicAMRHierarchyHDF5(fname4d,
                                               vectGrids4d,
                                               vectData4d,
                                               vectNames,
                                               domain4d,
                                               dx4d,
                                               dt,
                                               time,
                                               refRatio4d,
                                               numLevels);


  // now, define 3D dataholders
  D3::Vector<D3::DisjointBoxLayout> vectGrids3d(vectGrids4d.size());
  D3::Vector<D3::LevelData<D3::FArrayBox>* >  vectData3d(vectData4d.size(), NULL);


  // now do the magic slicing from 4d->3d
  D3::Box domain3d;
  D4::SliceSpec firstSlice(slice_dir[0],slice_comp[0]);
  bool outOfBounds = false;
  sliceBox(domain3d, domain4d, firstSlice, &outOfBounds);
  
  if (outOfBounds) 
    {
      MayDay::Abort("slice out of bounds in domain");
    }

  Vector<D3::IntVect> refRatio3d(refRatio4d.size());  
  D3::RealVect dx3d;
  // SliceRealVect doesn't exist, so do it ourselves
  //sliceRealVect(dx3d, dx4d, firstSlice);
  
  D4::IntVect iv4D;;
  for (int dir=0; dir<4; dir++) iv4D[dir] = dir;
  D3::IntVect iv3D;
  sliceIntVect(iv3D, iv4D, firstSlice);
  // use sliced IV to inform slicing of RealVect
  for (int dir=0; dir<3; dir++) dx3d[dir] = dx4d[iv3D[dir]];
  
  
  for (int lev=0; lev<refRatio4d.size(); lev++)
    {
      sliceIntVect(refRatio3d[lev], refRatio4d[lev], firstSlice);
    }

  // for now, assume no AMR (revisit this when we need to)
  CH_assert(numLevels == 1);

  for (int lev=0; lev<vectGrids4d.size(); lev++)
    {
      bool maintainProcAssign = true;
      sliceDisjointBoxLayout(vectGrids3d[lev],
                             vectGrids4d[lev],
                             firstSlice, 
                             maintainProcAssign);
      vectData3d[lev] = new D3::LevelData<D3::FArrayBox>;
      sliceLevelDataLocalOnly(*vectData3d[lev], *vectData4d[lev], 
                              firstSlice);
    }

  
  // if we're only slicing to 3d, then stop here and write plotfile
  if (slice_dir.size() == 1)
    {
      D3::WriteAnisotropicAMRHierarchyHDF5(slicefname,
                                           vectGrids3d,
                                           vectData3d,
                                           vectNames,
                                           domain3d,
                                           dx3d,
                                           dt,
                                           time,
                                           refRatio3d,
                                           numLevels);

    }
  else
    {
      // continue on to slice down to 2D

      // define 2D dataholdersssssss
      D2::Vector<D2::DisjointBoxLayout> vectGrids2d(vectGrids3d.size());
      D2::Vector<D2::LevelData<D2::FArrayBox>* >  vectData2d(vectData3d.size(),
                                                             NULL);
      D2::Box domain2d;
      D3::SliceSpec secondSlice(slice_dir[1],slice_comp[1]);
      sliceBox(domain2d, domain3d, secondSlice);

      Vector<D2::IntVect> refRatio2d(refRatio3d.size());  
      D2::RealVect dx2d;
      //sliceRealVect(dx2d, dx3d, secondSlice);

      for (int dir=0; dir<3; dir++) iv3D[dir] = dir;
      D2::IntVect iv2D;
      sliceIntVect(iv2D, iv3D, secondSlice);
      // use sliced IV to inform slicing of RealVect
      for (int dir=0; dir<2; dir++) dx2d[dir] = dx3d[iv2D[dir]];

      for (int lev=0; lev<refRatio3d.size(); lev++)
        {
          sliceIntVect(refRatio2d[lev], refRatio3d[lev], secondSlice);
        }

      
      
      for (int lev=0; lev<vectGrids3d.size(); lev++)
        {
          sliceDisjointBoxLayout(vectGrids2d[lev],
                                 vectGrids3d[lev],
                                 secondSlice);
          vectData2d[lev] = new D2::LevelData<D2::FArrayBox>;
          sliceLevelDataLocalOnly(*vectData2d[lev], *vectData3d[lev], 
                                  secondSlice);
        }

      
      D2::WriteAnisotropicAMRHierarchyHDF5(slicefname,
                                           vectGrids2d,
                                           vectData2d,
                                           vectNames,
                                           domain2d,
                                           dx2d,
                                           dt,
                                           time,
                                           refRatio2d,
                                           numLevels);

    }
  
  
  return status;
}
