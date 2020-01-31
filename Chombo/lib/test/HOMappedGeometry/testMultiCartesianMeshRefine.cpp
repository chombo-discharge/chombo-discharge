#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::ifstream;
using std::ios;
#include <string>

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"

#include "MultiBlockCoordSys.H"
#include "MultiBlockMeshRefine.H"
#include "MultiCartesianCS.H"
#include "SPMD.H"
#include "newMappedGridIO.H"

#include "CH_Attach.H"

#include "UsingNamespace.H"


void testMCMultiBlockMeshRefine();

void setFunc(const Box&, int m, FArrayBox& F)
{
  F.setVal(procID());
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
  setChomboMPIErrorHandler();
#endif

  // Check for an input file
  char* inFile = NULL;


  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage: " << a_argv[0] << " <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  // Run the tests
  testMCMultiBlockMeshRefine();

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
 MPI_Finalize();
#endif


}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void tagCells(IntVectSet& a_tags, const IntVect& a_center, int a_radius)
{
  IntVect lo = (a_center - a_radius)*IntVect::Unit;
  IntVect hi = (a_center + a_radius)*IntVect::Unit;
  Box boxTags(lo,hi);
  a_tags |= boxTags;
}

void testMCMultiBlockMeshRefine()
{

  // Read ParmParse inputs
  ParmParse pp;

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  pp.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  pp.query("domain_length",domainLength);

  // Set the resolution of the coarsest level
  int numCells = 32;
  pp.query("num_cells", numCells);
  CH_assert(numCells > 0);

  // Determine how many blocks in each dimension.
  // Default is 2 in every dimension except the last, where it is 1.
  IntVect count = 2*IntVect::Unit - BASISV(SpaceDim-1);
  {
    Vector<int> counta(SpaceDim);
    if (pp.queryarr("count", counta, 0, SpaceDim))
      {
        count = IntVect(counta);
      }
  }

  // Determine which spatial directions wrap around.
  // Default is all dimensions wrap around except the highest.
  IntVect wraparound = IntVect::Unit - BASISV(SpaceDim-1);
  {
    Vector<int> wraparounda(SpaceDim);
    if (pp.queryarr("wraparound", wraparounda, 0, SpaceDim))
      {
        wraparound = IntVect(wraparounda);
      }
  }

  // Determine which directions are to be refined.
  // Default is all of them.
  IntVect refineDirs = IntVect::Unit;
  {
    Vector<int> refineDirsa(SpaceDim);
    if (pp.queryarr("refine_dirs", refineDirsa, 0, SpaceDim))
      {
        refineDirs = IntVect(refineDirsa);
      }
  }

  // Maximum AMR level limit
  int maxLevel = 0;
  pp.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  pp.queryarr("ref_ratio", refRatios, 0, numReadLevels+1);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  pp.query("tag_buffer_size",tagBufferSize);

  // Minimum dimension of a grid
  int blockFactor = 1;
  pp.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  pp.query("max_grid_size",maxGridSize);

  Vector<int> sh(CH_SPACEDIM,1);
  pp.queryarr("shift", sh, 0, CH_SPACEDIM);
  IntVect shift(&(sh[0]));

  Real fillRatio = 0.75;
  pp.query("fill_ratio",fillRatio);

  int gridBufferSize = 1;
  pp.query("grid_buffer_size",gridBufferSize);

  // Print the parameters
  pout() << "count = " << count << endl;
  pout() << "wraparound = " << wraparound << endl;
  pout() << "refine_dirs = " << refineDirs << endl;
  pout() << "number of cells = " << numCells << endl;
  pout() << "maximum level = " << maxLevel << endl;
  pout() << "refinement ratio = ";
  for (int i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
  pout() << endl;
  pout() << "blocking factor = " << blockFactor << endl;
  pout() << "max grid size = " << maxGridSize << endl;
  pout() << "shift = " <<shift <<endl;
  pout() << "fill ratio = " << fillRatio << endl;
  pout() << "grid_buffer_size = "<<gridBufferSize << endl;

  int radius = 1;
  pp.query("patch_radius",radius);

  // Create coordinate system factory and initialize
  // Set up the coordinate system... factory
  Real dx = domainLength / (Real(numCells));
  RealVect dxVect = dx * RealVect::Unit;
  Vector<int> ratios(refRatios);

  IntVect levelDomainLength = numCells * (2*count - IntVect::Unit);

  RefCountedPtr<MultiCartesianCSFactory> coordSysFactPtr;
  coordSysFactPtr = RefCountedPtr<MultiBlockCoordSysFactory>(new MultiCartesianCSFactory);
  coordSysFactPtr->setCount(count);
  coordSysFactPtr->setWraparound(wraparound);

  Box levelDomainBox(IntVect::Zero, levelDomainLength-IntVect::Unit);
  ProblemDomain levelDomain(levelDomainBox);
  Vector<const MultiBlockCoordSys*> coords(maxLevel+1);
  ProblemDomain pd(levelDomain);
  RealVect dv(dxVect);
  for(int r=0; r<=maxLevel; r++)
    {
      coords[r] = coordSysFactPtr->getCoordSys(pd, dv);
      dv /= (refRatios[r]);
      pd.refine(refRatios[r]);
    }

  MultiBlockMeshRefine mbmr(levelDomain, dxVect, coordSysFactPtr, refRatios, fillRatio, blockFactor, gridBufferSize, maxGridSize) ;
  mbmr.setRefineDirs(refineDirs);

  const Vector<Box>& blockBoxes = coords[0]->mappingBlocks();
  int nblocks = blockBoxes.size();
  Vector<Box> allBoxes;
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      Box blockBox = blockBoxes[iblock];
      mbmr.restrictUnrefined(blockBox);
      Vector<Box> thisBlockBoxes;
      domainSplit(blockBox, thisBlockBoxes,
                  maxGridSize, blockFactor, refineDirs);
      allBoxes.append(thisBlockBoxes);
    }

  // Test the multiblock grid generation in regrid
  Vector<Vector<Box> > newGrids(maxLevel);
  IntVectSet tags;
  //  // Tag the full bottom layer in block 0.
  //  Box taggedBox(blockBoxes[0]);
  //  taggedBox.setBig(SpaceDim-1, 0);
  IntVect center0 = (blockBoxes[0].bigEnd()-shift) * refineDirs;
  tagCells(tags, center0, radius);
  if (blockBoxes.size() > 1)
    {
      IntVect center1 = (blockBoxes[1].smallEnd()+shift*2) * refineDirs;
      tagCells(tags, center1, radius);
      if (blockBoxes.size() > 3)
        {
          IntVect center3 = (blockBoxes[3].bigEnd()- shift*3) * refineDirs;
          tagCells(tags, center3, radius);
        }
    }
  //  tags |= taggedBox;
  mbmr.restrictUnrefined(tags);

  pout() << "Tags min box: " << tags.minBox() << endl;
  int baseLevel = 0;
  int topLevel = maxLevel-1;

  Vector<Vector<Box> > oldGrids(topLevel+1);
  oldGrids[0] = allBoxes;
  mbmr.regrid(newGrids, tags, baseLevel, topLevel, oldGrids);

  bool fullRange = false;
  if (fullRange && (refineDirs != IntVect::Unit))
    { // Expand newGrids so they take up full range in unrefined dimensions.
      for (int lev=0; lev<=maxLevel; lev++)
        {
          Vector<Box>& newGridsLevel = newGrids[lev];
          const Vector<Box>& blockBoxes = coords[lev]->mappingBlocks();
          // ASSUMES blockBoxes[0] has full range in unrefined dimensions.
          const Box& baseBox = blockBoxes[0];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (refineDirs[idir] == 0)
                {
                  int start = baseBox.smallEnd(idir);
                  int len = baseBox.size(idir);
                  for (int ibox = 0; ibox < newGridsLevel.size(); ibox++)
                    {
                      newGridsLevel[ibox].setRange(idir, start, len);
                    }
                }
            }
        }
    }

  Vector<DisjointBoxLayout> db(newGrids.size());
  Vector<std::string> names(1, "phi");
  {
    CH_TIME("dumpFiles");
    Vector<LevelData<FArrayBox>*> data(newGrids.size());
    for(int l = 0; l<newGrids.size(); ++l)
      {
        pout() << "Boxes at level " << l << ":  " << newGrids[l].size() << endl;
        mortonOrdering(newGrids[l]);
        Vector<int> procs;
        LoadBalance(procs, newGrids[l]);
        DisjointBoxLayout dbl(newGrids[l], procs);
        db[l] = dbl;
        data[l] = new LevelData<FArrayBox>(dbl, 1);
        data[l]->apply(setFunc);
      }
    // data is the processor number.
    WriteMappedAMRHierarchyHDF5(std::string("meshes"), db, data, names, coords, levelDomainBox, 1, 1, refRatios, db.size());
    
    for(int l = 0; l<newGrids.size(); ++l)
      {
        delete data[l];
      }
  }
  
  for(int l = 0; l<newGrids.size(); ++l)
    {
      delete coords[l];
    }
}
