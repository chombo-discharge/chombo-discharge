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


void setProc(const Box&, int m, FArrayBox& F)
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

Real speed = 1.0/20;
Real radius = 1.0;
Real subradius = radius/2;
Real dx = 1;

int STEP = -1;
Real gradeFunction(const RealVect& a_x, int a_step, Real a_tolerance)
{
  static Real c[3];
  if(a_step != STEP)
    {
      STEP = a_step;
      Real theta = 2*M_PI*STEP*speed;
      c[0] = radius*cos(theta);
      c[1] = radius*sin(theta);
      c[2] = 2*dx*speed*STEP;
    }
  Real x=c[0]-a_x[0];
  Real y=c[1]-a_x[1];
  Real d=x*x+y*y;
#if CH_SPACEDIM == 3
  Real z=c[2]-a_x[2];
  d+= z*z;
#endif
  d=sqrt(d);
  Real s = d/subradius;
  Real g = a_tolerance*(dx*s*s -2.3*dx*s + 1.6*dx);
  if(g>dx) g = dx;
  return g;
}

void tagCells(IntVectSet& a_tags, const Box& a_box, const MultiBlockCoordSys* a_coords, 
              Real a_tolerance, int a_step)
{
  int block = a_coords->whichBlock(a_box);
  const NewCoordSys* map = a_coords->getCoordSys(block);
  FArrayBox x(a_box, 3);
  map->getCenterMappedCoordinates(x, a_box);
  const RealVect& dx = map->dx();
  Real m = dx[0];
  for(int i=1; i<CH_SPACEDIM; i++)
    {
      if(dx[i]>m) m=dx[1];
    }
  RealVect p;
  for(BoxIterator bit=BoxIterator(a_box); bit.ok(); ++bit)
    {
      const IntVect& i=bit();
      //D_EXPR(p[0]=x(i,0), p[1]=x(i,1), p[2]=x(i,2));
      p = map->realCoord(map->centerMappedCoordinates(i));
      Real g = gradeFunction(p, a_step, a_tolerance);
      if(m>g) 
        {
          a_tags|=i;
        }
    }
}
 
void setValues(FArrayBox& f,  int a_step, int a_tolerance, const MultiBlockCoordSys* a_coords)
{
  const Box& b = f.box();
  f.setVal(procID(), 1);
  int block = a_coords->whichBlock(b);
  const NewCoordSys* map = a_coords->getCoordSys(block);
  FArrayBox x(b, 3);
  map->getCenterMappedCoordinates(x, b);
  const RealVect& dx = map->dx();
  Real m = dx[0];
  for(int i=1; i<CH_SPACEDIM; i++)
    {
      if(dx[i]>m) m=dx[1];
    }
  RealVect p;
  for(BoxIterator bit=BoxIterator(b); bit.ok(); ++bit)
    {
      const IntVect& i=bit();
      //D_EXPR(p[0]=x(i,0), p[1]=x(i,1), p[2]=x(i,2));
      p = map->realCoord(map->centerMappedCoordinates(i));
      Real g = gradeFunction(p, a_step, a_tolerance);
      f(i,0)=g;
    }
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
  int nsteps = 0;
  pp.query("max_step",nsteps);

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

  // Maximum AMR level limit
  int maxLevel = 0;
  pp.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  std::string prefix("plot.");
  if(pp.contains("prefix")) pp.get("prefix", prefix);

  if(pp.contains("speed")) pp.get("speed", speed);
  if(pp.contains("radius")) pp.get("radius", radius);
  if(pp.contains("subradius")) pp.get("subradius", subradius);
  
  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  pp.queryarr("ref_ratio",refRatios,0,numReadLevels+1);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  pp.query("tag_buffer_size",tagBufferSize);

  Real dt = 0.1;
  if(pp.contains("dt")) pp.get("dt",dt);

  // tolerance sets the tagging threshold.  a cell is tagged based in "tolerance" < sum(df^2/dx^2) * dx^2
  Real tolerance = 0.01;
  if(pp.contains("tolerance")) pp.get("tolerance",tolerance);
  
  // Minimum dimension of a grid
  int blockFactor = 1;
  pp.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  pp.query("max_grid_size",maxGridSize);
  int maxBaseGridSize = 0;
  pp.query("max_base_grid_size",maxBaseGridSize);

  Real fillRatio = 0.75;
  pp.query("fill_ratio",fillRatio);

  int gridBufferSize = 1;
  pp.query("grid_buffer_size",gridBufferSize);

  // Print the parameters
  pout() << "maximum step = " << nsteps << endl;
  pout() << "number of cells =";
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      pout() << " " << numCells[idir];
    }
  pout() << endl;
  pout() << "maximum level = " << maxLevel << endl;
  pout() << "tolerance = " << tolerance << endl;
  pout() << "refinement ratio = ";
  for (int i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
  pout() << endl;

  pout() << "blocking factor = " << blockFactor << endl;
  pout() << "max grid size = " << maxGridSize << endl;
  pout() << "max base grid size = " << maxBaseGridSize << endl ;
  pout() << "fill ratio = " << fillRatio << endl;
  pout() << "grid_buffer_size = "<<gridBufferSize << endl;
  int radius = 1;
  pp.query("patch_radius",radius);

  // Create coordinate system factory and initialize
  // Set up the coordinate system... factory
  dx = domainLength / (Real(numCells[0]));
  RealVect dxVect = dx * RealVect::Unit;
  Vector<int> ratios(refRatios);

  IntVect levelDomainLo, levelDomainHi;

  levelDomainLo = IntVect::Zero;
  levelDomainHi = (3*lengthCells-1) * IntVect::Unit; // fine for DoubleCartesianCS
  
  RefCountedPtr<MultiBlockCoordSysFactory> coordSysFactPtr;
  levelDomainLo = -2*lengthCells*IntVect::Unit;


  //levelDomainLo = IntVect(D_TERM(-2*lengthCells, -2*lengthCells, 0));
  levelDomainHi = (3*lengthCells-1)*IntVect::Unit;
  CylinderEquiangularCSFactory* f = new CylinderEquiangularCSFactory;
  f->setCenterPoint(RealVect::Zero);
  RealVect sizeV= 2*RealVect::Unit;
  f->setCentralRectSize(sizeV);
  f->setOuterRadius(3);
  coordSysFactPtr =  RefCountedPtr<MultiBlockCoordSysFactory>(f);
   


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

  //  OK, base grid is set up now.  Build the MeshGenerator

  MultiBlockMeshRefine mbmr(levelDomain, dxVect ,coordSysFactPtr, refRatios, fillRatio, blockFactor, gridBufferSize, maxGridSize) ;


  // Set up initial grids
  Vector<Vector<Box> > newGrids(maxLevel);
  IntVectSet tags;
 
  int baseLevel = 0;
  int topLevel = maxLevel-1;
  Vector<Vector<Box> > oldGrids(topLevel+1);
  oldGrids[0] = allBoxes;
  Vector<int> procs;
  LoadBalance(procs, allBoxes);
  DisjointBoxLayout dbl(allBoxes, procs);
  
  for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      tagCells(tags, dbl[dit], coords[0], tolerance, 0); // create initial distributed tag set, just on base level
    }
  
  mbmr.regrid(newGrids, tags, baseLevel, topLevel, oldGrids);

 

  Vector<DisjointBoxLayout> db(newGrids.size());
  Vector<LevelData<FArrayBox>*> data(newGrids.size());
  Vector<std::string> names(2);
  names[0] = std::string("grading");
  names[1] = std::string("processor");

 
  for(int l = 0; l<newGrids.size(); ++l)
    {
      mortonOrdering(newGrids[l]);
      Vector<int> procs;
      LoadBalance(procs, newGrids[l]);
      DisjointBoxLayout dbl(newGrids[l], procs);
      db[l] = dbl;
      data[l] = new LevelData<FArrayBox>(dbl, 2);
      LevelData<FArrayBox>& ld = *(data[l]);
      for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) setValues(ld[dit], 0, tolerance, coords[l]);
      if(verbosity >=1)
        pout()<<" level "<<l<<" elements:"<<dbl.numCells();
    }
  if(verbosity >=1)
    pout()<<"\n";

  char pcount[5]="000.";
  std::string fname=prefix; fname.append(pcount);
  
  WriteMappedAMRHierarchyHDF5(fname, db, data, names, coords, levelDomainBox, 1, 0, refRatios, db.size());
  
  
  for(int l = 0; l<newGrids.size(); ++l)
    {
      delete data[l];
    }
  Vector<IntVectSet> tagv(topLevel+1);
  for(int step = 1; step<nsteps; ++step)
    {
     if(verbosity >=1)
       pout()<<"step "<<step<<" ";

      for(int i=0 ; i<=topLevel; i++)
        {
          const DisjointBoxLayout& dbl = db[i];
          IntVectSet& tags = tagv[i];
          tags.makeEmpty();
          for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
            {
              tagCells(tags, dbl[dit], coords[i], tolerance, step);
            }
        }
      
      {
        CH_TIME("Gridding");
        mbmr.regrid(newGrids, tagv, 0, topLevel, oldGrids);
      }
      db.resize(newGrids.size());
      for(int l=0; l<newGrids.size(); l++)
        {
          mortonOrdering(newGrids[l]);
          Vector<int> procs;
          LoadBalance(procs, newGrids[l]);
          DisjointBoxLayout dbl(newGrids[l], procs);
          db[l] = dbl;
          data[l] = new LevelData<FArrayBox>(dbl, 2);
          LevelData<FArrayBox>& ld = *(data[l]);
          for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) setValues(ld[dit], step, tolerance, coords[l]);
          if(verbosity>=1)
            pout()<<" level "<<l<<" elements:"<<dbl.numCells();
        }
      if(verbosity >=1)
        pout()<<"\n";
      sprintf(pcount, "%03i.", step);
      fname=prefix; fname.append(pcount);
      WriteMappedAMRHierarchyHDF5(fname, db, data, names, coords, levelDomainBox, 1, step, refRatios, db.size());
      for(int l=0; l<newGrids.size(); l++)
        {
          delete data[l];
        }
    }
  
  for(int i=0; i<coords.size(); i++) delete coords[i];
  
}
