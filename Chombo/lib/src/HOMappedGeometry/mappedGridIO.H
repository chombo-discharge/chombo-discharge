#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MAPPEDGRIDIO_H_
#define _MAPPEDGRIDIO_H_

#include "FourthOrderCoordSys.H"
#include "AMRIO.H"
#include "NodeAMRIO.H"
#include "REAL.H"
#include <string>
#include "NamespaceHeader.H"
using std::string;
using std::istream;


#ifdef CH_USE_HDF5
/// Write single-level of mapped-grid data in HDF5 format
/**
   Writes a single level in HDF5 format.  Only available if the preprocessor
   macro HDF5 is defined at compilation.  This is blocking.

   Opens and writes to 2 files:
   fileRoot.<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot.map.<DIM>d.hdf5 contains the node-centered mapped-grid information
   (which is the real-space locations of the nodes in mapped space.)

   \b Arguments:
   - \b a_fileRoot: root of files to output to.
   - \b a_grids   : grids
   - \b a_data    : data
   - \b a_names   : names of variables.
   - \b a_CoordSys: coordinate system
*/
void WriteMappedUGHDF5(const string&               a_fileRoot,
                       const DisjointBoxLayout&    a_grids,
                       const LevelData<FArrayBox>& a_data,
                       const FourthOrderCoordSys&  a_CoordSys);




///
/**
   Writes hierarchy of mapped levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot.<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot.map.<DIM>d.hdf5 contains the node-centered mapped-grid information
   (which is the real-space locations of the nodes in mapped space.)

   {\bf Arguments:}\\
   fileRoot  :  root of files to output to.\\
   a_vectGrids : grids at each level.\\
   a_vectData :  data at each level.\\
   a_vectNames:  names of variables.\\
   a_vectCoordSys: coordinate systems for each AMR level. \\
   a_dt     :  time step at coarsest level.\\
   a_time     :  time.\\
   a_vectRatio :  refinement ratio at all levels
   (ith entry is refinement ratio between levels i and i + 1).\\
   a_numLevels :  number of levels to output.\\

   This is blocking.

*/
void
WriteMappedAMRHierarchyHDF5(const string& fileRoot,
                            const Vector<DisjointBoxLayout>& a_vectGrids,
                            const Vector<LevelData<FArrayBox>* > & a_vectData,
                            const Vector<string>& a_vectNames,
                            const Vector<const FourthOrderCoordSys* >& a_vectCoordsys,
                            const Real& a_dt,
                            const Real& a_time,
                            const Vector<int>& a_vectRatio,
                            const int& a_numLevels);


///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot.<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot.map.<DIM>d.hdf5 contains the node-centered mapped-grid information
   (which is the real-space locations of the nodes in mapped space.)

   {\bf Arguments:}\\
   fileRoot  :  root of files to output to.\\
   a_vectGrids : grids at each level.\\
   a_vectData :  data at each level.\\
   a_vectCoordSys: coordinate systems for each AMR level. \\
   a_vectRatio :  refinement ratio at all levels
   (ith entry is refinement ratio between levels i and i + 1).\\
   a_numLevels :  number of levels to output.\\
   This is blocking.
*/
void
WriteMappedAMRHierarchyHDF5(const string& fileRoot,
                            const Vector<DisjointBoxLayout>& a_vectGrids,
                            const Vector<LevelData<FArrayBox>* > & a_vectData,
                            const Vector<const FourthOrderCoordSys* >& a_vectCoordsys,
                            const Vector<int>& a_vectRatio,
                            const int& a_numLevels);


#endif // CH_USE_HDF5


#include "NamespaceFooter.H"

#endif
