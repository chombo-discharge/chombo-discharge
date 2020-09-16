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
#include "DoubleCartesianCS.H"
#include "CubedSphere2DCS.H"
#include "CylinderEquiangularCS.H"
#include "SPMD.H"
#include "newMappedGridIO.H"

#include "CH_Attach.H"

#include "UsingNamespace.H"


void testMBMeshRefine();

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
  testMBMeshRefine();

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

void testMBMeshRefine()
{

  // Read ParmParse inputs
  ParmParse pp;

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  pp.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Stop after this number of steps
  int nstop = 0;
  pp.query("max_step",nstop);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  pp.query("domain_length",domainLength);

  // Set the location of the lower left corner
  std::vector<Real> x0a(SpaceDim,0.0);
  pp.queryarr("x0",x0a,0,SpaceDim);
  RealVect x0;
  for ( int d=0 ; d<SpaceDim ; ++d ) x0[d] = x0a[d] ;

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  pp.queryarr("num_cells",numCells,0,SpaceDim);
  // Check that every component of numCells is positive
  for (int i = 0; i < SpaceDim; ++i)
    {
      CH_assert(numCells[i] > 0);
    }
  // IGNORE dimensions other than the first one.
  int lengthCells = numCells[0];

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  pp.queryarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      isPeriodic[dim] = (isPeriodica[dim] == 1);
      if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
        pout() << "Using Periodic BCs in direction: " << dim << endl;
    }

  // Maximum AMR level limit
  int maxLevel = 0;
  pp.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  pp.queryarr("ref_ratio",refRatios,0,numReadLevels+1);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  pp.query("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  pp.query ("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  pp.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  pp.query("max_grid_size",maxGridSize);
  int maxBaseGridSize = 0;
  pp.query("max_base_grid_size",maxBaseGridSize);

  Vector<int> sh(CH_SPACEDIM,1);
  pp.queryarr("shift", sh, 0, CH_SPACEDIM);
  IntVect shift(&(sh[0]));

  Real fillRatio = 0.75;
  pp.query("fill_ratio",fillRatio);

  int gridBufferSize = 1;
  pp.query("grid_buffer_size",gridBufferSize);

  // Print the parameters
  pout() << "maximum step = " << nstop << endl;
  pout() << "number of cells =";
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      pout() << " " << numCells[idir];
    }
  pout() << endl;
  pout() << "maximum level = " << maxLevel << endl;
  pout() << "refinement ratio = ";
  for (int i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
  pout() << endl;
  pout() << "refinement threshold = " << refineThresh << endl;
  pout() << "blocking factor = " << blockFactor << endl;
  pout() << "max grid size = " << maxGridSize << endl;
  pout() << "shift = " <<shift <<endl;
  pout() << "max base grid size = " << maxBaseGridSize << endl ;
  pout() << "fill ratio = " << fillRatio << endl;
  pout() << "grid_buffer_size = "<<gridBufferSize << endl;
  int radius = 1;
  pp.query("patch_radius",radius);

  // Create coordinate system factory and initialize
  // Set up the coordinate system... factory
  Real dx = domainLength / (Real(numCells[0]));
  RealVect dxVect = dx * RealVect::Unit;
  Vector<int> ratios(refRatios);

  IntVect levelDomainLo, levelDomainHi;

  levelDomainLo = IntVect::Zero;
  levelDomainHi = (3*lengthCells-1) * IntVect::Unit; // fine for DoubleCartesianCS
  
  RefCountedPtr<MultiBlockCoordSysFactory> coordSysFactPtr;
#if CH_SPACEDIM == 2
  //coordSysFactPtr = RefCountedPtr<MultiBlockCoordSysFactory>(new DoubleCartesianCSFactory);
  coordSysFactPtr = RefCountedPtr<MultiBlockCoordSysFactory>(new CubedSphere2DCSFactory);
  levelDomainHi = IntVect(11*lengthCells-1, lengthCells-1);
#endif    
#if CH_SPACEDIM == 3
      levelDomainLo = IntVect(-2*lengthCells, -2*lengthCells, 0);
      levelDomainHi = IntVect(3*lengthCells-1, 3*lengthCells-1, 3*lengthCells-1);
      CylinderEquiangularCSFactory* f = new CylinderEquiangularCSFactory;
      f->setCenterPoint(RealVect::Zero);
      RealVect sizeV(2,2,3);
      f->setCentralRectSize(sizeV);
      f->setOuterRadius(3);
      coordSysFactPtr =  RefCountedPtr<MultiBlockCoordSysFactory>(f);
#endif    


  Box levelDomainBox(levelDomainLo, levelDomainHi);
  ProblemDomain levelDomain(levelDomainBox);
  Vector<const MultiBlockCoordSys*> coords(maxLevel+1);
  ProblemDomain d(levelDomain);
  RealVect dv(dxVect);
  for(int r=0; r<=maxLevel; r++)
    {
      coords[r] = coordSysFactPtr->getCoordSys(d, dv);
      dv/=(refRatios[r]);
      d.refine(refRatios[r]);
    }

  const Vector<Box>& blockBoxes = coords[0]->mappingBlocks();
  int nblocks = blockBoxes.size();
  Vector<Box> allBoxes;
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      Vector<Box> thisBlockBoxes;
      domainSplit(blockBoxes[iblock], thisBlockBoxes,
                  maxGridSize, blockFactor);
      allBoxes.append(thisBlockBoxes);
    }


  MultiBlockMeshRefine mbmr(levelDomain, dxVect ,coordSysFactPtr, refRatios, fillRatio, blockFactor, gridBufferSize, maxGridSize) ;


  // Test the multiblock grid generation in regrid
  Vector<Vector<Box> > newGrids(maxLevel);
  IntVectSet tags;
  tagCells(tags, blockBoxes[0].bigEnd()-shift, radius);
  tagCells(tags, blockBoxes[1].smallEnd()+shift*2, radius);
  tagCells(tags, blockBoxes[3].bigEnd()- shift*3, radius);
  pout() << "Tags min box: " << tags.minBox() << endl;
  int baseLevel = 0;
  int topLevel = maxLevel-1;
  Vector<Vector<Box> > oldGrids(topLevel+1);
  oldGrids[0] = allBoxes;
  mbmr.regrid(newGrids, tags, baseLevel, topLevel, oldGrids);

 
  Vector<DisjointBoxLayout> db(newGrids.size());
  Vector<std::string> names(1, "phi");
  {
    CH_TIME("dumpFiles1");  
  Vector<LevelData<FArrayBox>*> data(newGrids.size());
  for(int l = 0; l<newGrids.size(); ++l)
    {
      mortonOrdering(newGrids[l]);
      Vector<int> procs;
      LoadBalance(procs, newGrids[l]);
      DisjointBoxLayout dbl(newGrids[l], procs);
      db[l] = dbl;
      data[l] = new LevelData<FArrayBox>(dbl, 1);
      data[l]->apply(setFunc);
    }
  WriteMappedAMRHierarchyHDF5(std::string("meshes"), db, data, names, coords, levelDomainBox, 1, 1, refRatios, db.size());

  
  for(int l = 0; l<newGrids.size(); ++l)
    {
      delete data[l];
    }
  }


  oldGrids = newGrids;
  newGrids.resize(0);
  
  Vector<IntVectSet> tagv(maxLevel);
  
  for(DataIterator dit = db[1].dataIterator(); dit.ok(); ++dit)
    {
      tagv[1] |= db[1][dit];
    }

  tagv[1].shift(shift);
  tagv[1].grow(1);

  // base level 1, top level 1
  mbmr.regrid(newGrids, tagv, 1, topLevel, oldGrids);

  {
    CH_TIME("dumpFiles2"); 
  Vector<DisjointBoxLayout> dbnew(newGrids.size());
  Vector<LevelData<FArrayBox>*> data(newGrids.size());
  for(int l = 0; l<newGrids.size(); ++l)
    {
      mortonOrdering(newGrids[l]);
      Vector<int> procs;
      LoadBalance(procs, newGrids[l]);
      DisjointBoxLayout dbl(newGrids[l], procs);
      dbnew[l] = dbl;
      data[l] = new LevelData<FArrayBox>(dbl, 1, IntVect::Unit*gridBufferSize);
      data[l]->apply(setFunc);
    }
  WriteMappedAMRHierarchyHDF5("remeshes", dbnew, data, names, coords, levelDomainBox, 1, 1, refRatios, dbnew.size());

  
  for(int l = 0; l<newGrids.size(); ++l)
    {
      delete data[l];
      delete coords[l];
    }
  }
}
