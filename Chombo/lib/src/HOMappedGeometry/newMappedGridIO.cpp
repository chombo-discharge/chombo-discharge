#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DFMartin, March 22, 2008

#include "newMappedGridIO.H"
#include "AMRIO.H"
#include "NodeAMRIO.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"


#ifdef CH_USE_HDF5
/// Write single-level of mapped-grid data in HDF5 format
/**
   Writes a single level in HDF5 format.  Only available if the preprocessor
   macro HDF5 is defined at compilation.  This is blocking.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot<DIM>d.map.hdf5 contains the node-centered mapped-grid information
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
                       const NewCoordSys&  a_CoordSys,
                       const Box& a_domainBox,
                       const Real& a_time)
{
  // simplest thing to do here is to just call the AMR version...
  const Vector<DisjointBoxLayout> vectGrids(1, a_grids);
  const Vector<LevelData<FArrayBox>* > vectData(1, const_cast<LevelData<FArrayBox>* >(&a_data) );
  const Vector<const NewCoordSys*> vectCoordSys(1, &a_CoordSys);
  Vector<int> refRatio(1,1);
  int numLevels = 1;

  WriteMappedAMRHierarchyHDF5(a_fileRoot, vectGrids, vectData,
                              vectCoordSys, a_domainBox, refRatio,
                              numLevels, a_time);

}


///
/**
   Writes hierarchy of mapped levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot.<DIM>d.map.hdf5 contains the node-centered mapped-grid information

   Writes hierarchy of mapped levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot.<DIM>d.map.hdf5 contains the node-centered mapped-grid information
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
WriteMappedAMRHierarchyHDF5(const string& a_fileRoot,
                            const Vector<DisjointBoxLayout>& a_vectGrids,
                            const Vector<LevelData<FArrayBox>* > & a_vectData,
                            const Vector<string>& a_vectNames,
                            const Vector<const NewCoordSys*>& a_vectCoordSys,
                            const Box& a_baseDomain,
                            const Real& a_dt,
                            const Real& a_time,
                            const Vector<int>& a_vectRatio,
                            const int& a_numLevels)
{
  // first, write out data to "regular" hdf5 file:
  std::vector<char> iter_str(a_fileRoot.size()+25);
  sprintf(&(iter_str[0]), "%s%dd.hdf5", a_fileRoot.c_str(), SpaceDim);
  string dataFileName(&(iter_str[0]));

  RealVect dx = a_vectCoordSys[0]->dx();
  Vector<IntVect> refRatios;
  for(int i=0; i<a_vectRatio.size(); i++)
    { 
      refRatios.push_back(IntVect::Unit*a_vectRatio[i]);
    }
  WriteAnisotropicAMRHierarchyHDF5(dataFileName,
                        a_vectGrids,
                        a_vectData,
                        a_vectNames,
                        a_baseDomain,
                        dx,
                        a_dt,
                        a_time,
                        refRatios,
                        a_numLevels);


  // now create node-centered data for geometric info
  Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(a_numLevels, NULL);
  int dimension = a_vectCoordSys[0]->dimension();
  for (int level=0; level<a_numLevels; level++)
    {
      const DisjointBoxLayout& levelGrids = a_vectGrids[level];
      // use same ghosting as cell-centered data used
      IntVect ghostVect = a_vectData[level]->ghostVect();
      vectNodeLoc[level] = new LevelData<NodeFArrayBox>(levelGrids,
                                                        dimension,
                                                        ghostVect);

      LevelData<NodeFArrayBox>& levelNodeData = *vectNodeLoc[level];
      const NewCoordSys& levelCS = *(a_vectCoordSys[level]);
      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          NodeFArrayBox& thisNodeFAB = levelNodeData[dit];
          // node-centered FAB
          FArrayBox& thisFAB = thisNodeFAB.getFab();
          levelCS.getNodeRealCoordinates(thisFAB, thisFAB.box());
        } // end loop over grids

    } // end loop over levels

  // create names
  Vector<string> locationNames(dimension);
  for (int i=0; i<dimension; i++)
    {
      if (i ==0)
        {
          locationNames[0] = "x";
        }
      else if (i == 1)
        {
          locationNames[1] = "y";
        }
      else if (i == 2)
        {
          locationNames[2] = "z";
        }
      else if (i == 3)
        {
          locationNames[3] = "u";
        }
      else if (i == 4)
        {
          locationNames[4] = "v";
        }
      else if (i == 5)
        {
          locationNames[5] = "w";
        }
    }

  sprintf(&(iter_str[0]), "%s%dd.map.hdf5", a_fileRoot.c_str(), SpaceDim);
  string gridInfoFileName(&(iter_str[0]));

  // now call nodal WriteAMRHierarchy function...
  WriteAMRHierarchyHDF5(gridInfoFileName,
                        a_vectGrids,
                        vectNodeLoc,
                        locationNames,
                        a_baseDomain,
                        dx[0],
                        a_dt,
                        a_time,
                        a_vectRatio,
                        a_numLevels);

  // clean up after ourselves here
  for (int level=0; level<a_numLevels; level++)
    {
      if (vectNodeLoc[level] != NULL)
        {
          delete vectNodeLoc[level];
          vectNodeLoc[level] = NULL;
        }
    }
}


///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot<DIM>d.map.hdf5 contains the node-centered mapped-grid information
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
WriteMappedAMRHierarchyHDF5(const string& a_fileRoot,
                            const Vector<DisjointBoxLayout>& a_vectGrids,
                            const Vector<LevelData<FArrayBox>* > & a_vectData,
                            const Vector<const NewCoordSys*>& a_vectCoordSys,
                            const Box& a_baseDomainBox,
                            const Vector<int>& a_vectRatio,
                            const int& a_numLevels,
                            const Real& a_time)
{

  // create names for the variables and placeholder values for dt and time,
  // then call the "real" one...
  Real dt = 1.0;
 
  int nComp = a_vectData[0]->nComp();
  Vector<string> compNames(nComp);
  for (int n=0; n<nComp; n++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", n);
      string label(labelChSt);
      compNames[n] = label;
    }

  WriteMappedAMRHierarchyHDF5(a_fileRoot,
                              a_vectGrids,
                              a_vectData,
                              compNames,
                              a_vectCoordSys,
                              a_baseDomainBox,
                              dt,
                              a_time,
                              a_vectRatio,
                              a_numLevels);
}


///
/**
   Writes hierarchy of mapped levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot.<DIM>d.map.hdf5 contains the node-centered mapped-grid information
   (which is the real-space locations of the nodes in mapped space.)

   {\bf Arguments:}\\
   fileRoot  :  root of files to output to.\\
   a_vectGrids : grids at each level.\\
   a_vectData :  data at each level.\\
   a_vectNames:  names of variables.\\
   a_vectCoordSys: coordinate systems for each AMR level. \\
   a_dt     :  time step at coarsest level.\\
   a_time     :  time.\\
   a_vectRatio :  refinement ratio in each direction at all levels
   (ith entry is refinement ratio between levels i and i + 1).\\
   a_numLevels :  number of levels to output.\\

   This is blocking.

*/
void
WriteAnisotropicMappedAMRHierarchyHDF5(
    const string& a_fileRoot,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Vector<const NewCoordSys*>& a_vectCoordSys,
    const Box& a_baseDomain,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_vectRatios,
    const int& a_numLevels)
{
  // first, write out data to "regular" hdf5 file:
  char iter_str[80];

  sprintf(iter_str, "%s%dd.hdf5", a_fileRoot.c_str(), SpaceDim);
  string dataFileName(iter_str);

  RealVect dx = a_vectCoordSys[0]->dx();

  WriteAnisotropicAMRHierarchyHDF5(
      dataFileName,
      a_vectGrids,
      a_vectData,
      a_vectNames,
      a_baseDomain,
      dx,
      a_dt,
      a_time,
      a_vectRatios,
      a_numLevels);


  // now create node-centered data for geometric info
  Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(a_numLevels, NULL);
  int dimension = a_vectCoordSys[0]->dimension();
  for (int level=0; level<a_numLevels; level++)
    {
      const DisjointBoxLayout& levelGrids = a_vectGrids[level];
      // use same ghosting as cell-centered data used
      IntVect ghostVect = a_vectData[level]->ghostVect();
      vectNodeLoc[level] = new LevelData<NodeFArrayBox>(levelGrids,
                                                        dimension,
                                                        ghostVect);

      LevelData<NodeFArrayBox>& levelNodeData = *vectNodeLoc[level];
      const NewCoordSys& levelCS = *(a_vectCoordSys[level]);
      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          NodeFArrayBox& thisNodeFAB = levelNodeData[dit];
          // node-centered FAB
          FArrayBox& thisFAB = thisNodeFAB.getFab();
          levelCS.getNodeRealCoordinates(thisFAB, thisFAB.box());
        } // end loop over grids

    } // end loop over levels

  // create names
  Vector<string> locationNames(dimension);
  for (int i=0; i<dimension; i++)
    {
      if (i ==0)
        {
          locationNames[0] = "x";
        }
      else if (i == 1)
        {
          locationNames[1] = "y";
        }
      else if (i == 2)
        {
          locationNames[2] = "z";
        }
      else if (i == 3)
        {
          locationNames[3] = "u";
        }
      else if (i == 4)
        {
          locationNames[4] = "v";
        }
      else if (i == 5)
        {
          locationNames[5] = "w";
        }
    }

  sprintf(iter_str, "%s%dd.map.hdf5", a_fileRoot.c_str(), SpaceDim);
  string gridInfoFileName(iter_str);

  // now call nodal WriteAMRHierarchy function...
  WriteAnisotropicAMRHierarchyHDF5(
      gridInfoFileName,
      a_vectGrids,
      vectNodeLoc,
      locationNames,
      a_baseDomain,
      dx,
      a_dt,
      a_time,
      a_vectRatios,
      a_numLevels);

  // clean up after ourselves here
  for (int level=0; level<a_numLevels; level++)
    {
      if (vectNodeLoc[level] != NULL)
        {
          delete vectNodeLoc[level];
          vectNodeLoc[level] = NULL;
        }
    }
}


///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot<DIM>d.map.hdf5 contains the node-centered mapped-grid information
   (which is the real-space locations of the nodes in mapped space.)

   {\bf Arguments:}\\
   fileRoot  :  root of files to output to.\\
   a_vectGrids : grids at each level.\\
   a_vectData :  data at each level.\\
   a_vectCoordSys: coordinate systems for each AMR level. \\
   a_vectRatio :  refinement ratio in each direction at all levels
   (ith entry is refinement ratio between levels i and i + 1).\\
   a_numLevels :  number of levels to output.\\
   This is blocking.
*/
void
WriteAnisotropicMappedAMRHierarchyHDF5(
    const string& a_fileRoot,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<const NewCoordSys*>& a_vectCoordSys,
    const Box& a_baseDomainBox,
    const Vector<IntVect>& a_vectRatios,
    const int& a_numLevels,
    const Real& a_time)
{

  // create names for the variables and placeholder values for dt and time,
  // then call the "real" one...
  Real dt = 1.0;
 
  int nComp = a_vectData[0]->nComp();
  Vector<string> compNames(nComp);
  for (int n=0; n<nComp; n++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", n);
      string label(labelChSt);
      compNames[n] = label;
    }

  WriteAnisotropicMappedAMRHierarchyHDF5(
      a_fileRoot,
      a_vectGrids,
      a_vectData,
      compNames,
      a_vectCoordSys,
      a_baseDomainBox,
      dt,
      a_time,
      a_vectRatios,
      a_numLevels);
}


/// Write single-level of multiblock mapped-grid data in HDF5 format
/**
   Writes a single level in HDF5 format.  Only available if the preprocessor
   macro HDF5 is defined at compilation.  This is blocking.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot<DIM>d.map.hdf5 contains the node-centered mapped-grid information
   (which is the real-space locations of the nodes in mapped space.)

   \b Arguments:
   - \b a_fileRoot: root of files to output to.
   - \b a_grids   : grids
   - \b a_data    : data
   - \b a_names   : names of variables.
   - \b a_CoordSys: multiblock coordinate system
*/
void WriteMappedUGHDF5(const string&               a_fileRoot,
                       const DisjointBoxLayout&    a_grids,
                       const LevelData<FArrayBox>& a_data,
                       const MultiBlockCoordSys&  a_CoordSys,
                       const Box& a_domainBox,
                       const Real& a_time)
{
  // simplest thing to do here is to just call the AMR version...
  const Vector<DisjointBoxLayout> vectGrids(1, a_grids);
  const Vector<LevelData<FArrayBox>* > vectData(1, const_cast<LevelData<FArrayBox>* >(&a_data) );
  const Vector<const MultiBlockCoordSys*> vectCoordSys(1, &a_CoordSys);
  Vector<int> refRatio(1,1);
  int numLevels = 1;

  WriteMappedAMRHierarchyHDF5(a_fileRoot, vectGrids, vectData,
                              vectCoordSys, a_domainBox, refRatio,
                              numLevels, a_time);

}


// this one takes a string of output names
void WriteMappedUGHDF5(const string&               a_fileRoot,
                       const DisjointBoxLayout&    a_grids,
                       const LevelData<FArrayBox>& a_data,
                       const MultiBlockCoordSys&  a_CoordSys,
                       const Box& a_domainBox,
                       const Vector<string>& outputNames,
                       const Real& a_time)
{
  // simplest thing to do here is to just call the AMR version...
  const Vector<DisjointBoxLayout> vectGrids(1, a_grids);
  const Vector<LevelData<FArrayBox>* > vectData(1, const_cast<LevelData<FArrayBox>* >(&a_data) );
  const Vector<const MultiBlockCoordSys*> vectCoordSys(1, &a_CoordSys);
  Vector<int> refRatio(1,1);
  int numLevels = 1;
  Real dt=1.0;


  WriteMappedAMRHierarchyHDF5(a_fileRoot, vectGrids, vectData, outputNames,
                              vectCoordSys, a_domainBox, dt, a_time, refRatio,
                              numLevels);

}

///
/**
   Writes hierarchy of multiblock mapped levels in HDF5 format.
   Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot.<DIM>d.map.hdf5 contains the node-centered mapped-grid information
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
WriteMappedAMRHierarchyHDF5(const string& a_fileRoot,
                            const Vector<DisjointBoxLayout>& a_vectGrids,
                            const Vector<LevelData<FArrayBox>* > & a_vectData,
                            const Vector<string>& a_vectNames,
                            const Vector<const MultiBlockCoordSys*>& a_vectCoordSys,
                            const Box& a_baseDomain,
                            const Real& a_dt,
                            const Real& a_time,
                            const Vector<int>& a_vectRatio,
                            const int& a_numLevels)
{
  // first, write out data to "regular" hdf5 file:
  char iter_str[80];

  sprintf(iter_str, "%s%dd.hdf5", a_fileRoot.c_str(), SpaceDim);
  string dataFileName(iter_str);

  const NewCoordSys* baseCoordSysPtr = a_vectCoordSys[0]->getCoordSys(0);
  RealVect dx = baseCoordSysPtr->dx();

  WriteAMRHierarchyHDF5(dataFileName,
                        a_vectGrids,
                        a_vectData,
                        a_vectNames,
                        a_baseDomain,
                        dx[0],
                        a_dt,
                        a_time,
                        a_vectRatio,
                        a_numLevels);


  // now create node-centered data for geometric info
  Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(a_numLevels, NULL);

  int dimension = baseCoordSysPtr->dimension();
  for (int level=0; level<a_numLevels; level++)
    {
      const DisjointBoxLayout& levelGrids = a_vectGrids[level];
      // use same ghosting as cell-centered data used
      IntVect ghostVect = a_vectData[level]->ghostVect();
      vectNodeLoc[level] = new LevelData<NodeFArrayBox>(levelGrids,
                                                        dimension,
                                                        ghostVect);

      LevelData<NodeFArrayBox>& levelNodeData = *vectNodeLoc[level];
      const MultiBlockCoordSys& levelCS = *(a_vectCoordSys[level]);
      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          NodeFArrayBox& thisNodeFAB = levelNodeData[dit];
          // node-centered FAB
          FArrayBox& thisFAB = thisNodeFAB.getFab();
          const NewCoordSys* blockCS = levelCS.getCoordSys(gridBox);

          blockCS->getNodeRealCoordinates(thisFAB, thisFAB.box());
        } // end loop over grids

    } // end loop over levels

  // create names

  // create names
  Vector<string> locationNames(dimension);
  for (int i=0; i<dimension; i++)
    {
      if (i ==0)
        {
          locationNames[0] = "x";
        }
      else if (i == 1)
        {
          locationNames[1] = "y";
        }
      else if (i == 2)
        {
          locationNames[2] = "z";
        }
      else if (i == 3)
        {
          locationNames[3] = "u";
        }
      else if (i == 4)
        {
          locationNames[4] = "v";
        }
      else if (i == 5)
        {
          locationNames[5] = "w";
        }
    }
  sprintf(iter_str, "%s%dd.map.hdf5", a_fileRoot.c_str(), SpaceDim);
  string gridInfoFileName(iter_str);

  // now call nodal WriteAMRHierarchy function...
  WriteAMRHierarchyHDF5(gridInfoFileName,
                        a_vectGrids,
                        vectNodeLoc,
                        locationNames,
                        a_baseDomain,
                        dx[0],
                        a_dt,
                        a_time,
                        a_vectRatio,
                        a_numLevels);

  // clean up after ourselves here
  for (int level=0; level<a_numLevels; level++)
    {
      if (vectNodeLoc[level] != NULL)
        {
          delete vectNodeLoc[level];
          vectNodeLoc[level] = NULL;
        }
    }
}


///
/**
   Writes hierarchy of multiblock mapped levels in HDF5 format.
   Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot.<DIM>d.map.hdf5 contains the node-centered mapped-grid information
   (which is the real-space locations of the nodes in mapped space.)

   {\bf Arguments:}\\
   fileRoot  :  root of files to output to.\\
   a_vectGrids : grids at each level.\\
   a_vectData :  data at each level.\\
   a_vectNames:  names of variables.\\
   a_vectCoordSys: coordinate systems for each AMR level. \\
   a_dt     :  time step at coarsest level.\\
   a_time     :  time.\\
   a_vectRatios:  refinement ratio in each direction at all levels
   (ith entry is refinement ratio between levels i and i + 1).\\
   a_numLevels :  number of levels to output.\\

   This is blocking.

*/
void
WriteAnisotropicMappedAMRHierarchyHDF5(
    const string& a_fileRoot,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Vector<const MultiBlockCoordSys*>& a_vectCoordSys,
    const Box& a_baseDomain,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_vectRatios,
    const int& a_numLevels)
{
  // first, write out data to "regular" hdf5 file:
  char iter_str[80];

  sprintf(iter_str, "%s%dd.hdf5", a_fileRoot.c_str(), SpaceDim);
  string dataFileName(iter_str);

  const NewCoordSys* baseCoordSysPtr = a_vectCoordSys[0]->getCoordSys(0);
  RealVect dx = baseCoordSysPtr->dx();

  WriteAnisotropicAMRHierarchyHDF5(
      dataFileName,
      a_vectGrids,
      a_vectData,
      a_vectNames,
      a_baseDomain,
      dx,
      a_dt,
      a_time,
      a_vectRatios,
      a_numLevels);


  // now create node-centered data for geometric info
  Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(a_numLevels, NULL);

  int dimension = baseCoordSysPtr->dimension();
  for (int level=0; level<a_numLevels; level++)
    {
      const DisjointBoxLayout& levelGrids = a_vectGrids[level];
      // use same ghosting as cell-centered data used
      IntVect ghostVect = a_vectData[level]->ghostVect();
      vectNodeLoc[level] = new LevelData<NodeFArrayBox>(levelGrids,
                                                        dimension,
                                                        ghostVect);

      LevelData<NodeFArrayBox>& levelNodeData = *vectNodeLoc[level];
      const MultiBlockCoordSys& levelCS = *(a_vectCoordSys[level]);
      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          NodeFArrayBox& thisNodeFAB = levelNodeData[dit];
          // node-centered FAB
          FArrayBox& thisFAB = thisNodeFAB.getFab();
          const NewCoordSys* blockCS = levelCS.getCoordSys(gridBox);

          blockCS->getNodeRealCoordinates(thisFAB, thisFAB.box());
        } // end loop over grids

    } // end loop over levels

  // create names

  // create names
  Vector<string> locationNames(dimension);
  for (int i=0; i<dimension; i++)
    {
      if (i ==0)
        {
          locationNames[0] = "x";
        }
      else if (i == 1)
        {
          locationNames[1] = "y";
        }
      else if (i == 2)
        {
          locationNames[2] = "z";
        }
      else if (i == 3)
        {
          locationNames[3] = "u";
        }
      else if (i == 4)
        {
          locationNames[4] = "v";
        }
      else if (i == 5)
        {
          locationNames[5] = "w";
        }
    }
  sprintf(iter_str, "%s%dd.map.hdf5", a_fileRoot.c_str(), SpaceDim);
  string gridInfoFileName(iter_str);

  // now call nodal WriteAMRHierarchy function...
  WriteAnisotropicAMRHierarchyHDF5(
      gridInfoFileName,
      a_vectGrids,
      vectNodeLoc,
      locationNames,
      a_baseDomain,
      dx,
      a_dt,
      a_time,
      a_vectRatios,
      a_numLevels);

  // clean up after ourselves here
  for (int level=0; level<a_numLevels; level++)
    {
      if (vectNodeLoc[level] != NULL)
        {
          delete vectNodeLoc[level];
          vectNodeLoc[level] = NULL;
        }
    }
}


///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot<DIM>d.map.hdf5 contains the node-centered mapped-grid information
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
WriteMappedAMRHierarchyHDF5(const string& a_fileRoot,
                            const Vector<DisjointBoxLayout>& a_vectGrids,
                            const Vector<LevelData<FArrayBox>* > & a_vectData,
                            const Vector<const MultiBlockCoordSys*>& a_vectCoordSys,
                            const Box& a_baseDomainBox,
                            const Vector<int>& a_vectRatio,
                            const int& a_numLevels,
                            const Real& a_time)
{

  // create names for the variables and placeholder values for dt and time,
  // then call the "real" one...
  Real dt = 1.0;
  int nComp = a_vectData[0]->nComp();
  Vector<string> compNames(nComp);
  for (int n=0; n<nComp; n++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", n);
      string label(labelChSt);
      compNames[n] = label;
    }

  WriteMappedAMRHierarchyHDF5(a_fileRoot,
                              a_vectGrids,
                              a_vectData,
                              compNames,
                              a_vectCoordSys,
                              a_baseDomainBox,
                              dt,
                              a_time,
                              a_vectRatio,
                              a_numLevels);
}


///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   Opens and writes to 2 files:
   fileRoot<DIM>d.hdf5 contains the cell-centered dataset, while
   fileRoot<DIM>d.map.hdf5 contains the node-centered mapped-grid information
   (which is the real-space locations of the nodes in mapped space.)

   {\bf Arguments:}\\
   fileRoot  :  root of files to output to.\\
   a_vectGrids : grids at each level.\\
   a_vectData :  data at each level.\\
   a_vectCoordSys: coordinate systems for each AMR level. \\
   a_vectRatio :  refinement ratio in each direction at all levels
   (ith entry is refinement ratio between levels i and i + 1).\\
   a_numLevels :  number of levels to output.\\
   This is blocking.
*/
void
WriteAnisotropicMappedAMRHierarchyHDF5(
    const string& a_fileRoot,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<const MultiBlockCoordSys*>& a_vectCoordSys,
    const Box& a_baseDomainBox,
    const Vector<IntVect>& a_vectRatios,
    const int& a_numLevels,
    const Real& a_time)
{

  // create names for the variables and placeholder values for dt and time,
  // then call the "real" one...
  Real dt = 1.0;
 
  int nComp = a_vectData[0]->nComp();
  Vector<string> compNames(nComp);
  for (int n=0; n<nComp; n++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", n);
      string label(labelChSt);
      compNames[n] = label;
    }

  WriteAnisotropicMappedAMRHierarchyHDF5(
      a_fileRoot,
      a_vectGrids,
      a_vectData,
      compNames,
      a_vectCoordSys,
      a_baseDomainBox,
      dt,
      a_time,
      a_vectRatios,
      a_numLevels);
}


#endif // CH_USE_HDF5

#include "NamespaceFooter.H"


